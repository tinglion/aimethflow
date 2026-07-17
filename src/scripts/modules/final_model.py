#!/usr/bin/env python3
"""
final_model.py — Step 15c: Train production models on the full cohort.

Step 15 (`classify`) evaluates models with cross-validation. This step re-fits
the same 11-model lineup on **all** retained samples (no CV holdout) and
persists the resulting artefacts to result/<wd>/model/, ready to drive a
deployment-time inference path (e.g., the WebUI diagnosis page).

Output layout
-------------
    result/<wd>/model/
    ├── manifest.json            metadata: model list, panel size, hyperparams,
    │                            training-cohort counts, source classifier_summary
    ├── preprocessing.joblib     {feature_names, col_medians, scaler,
    │                             ca125_median, ca125_scaler}
    ├── LR.joblib                cfDNA-only models
    ├── SVM.joblib
    ├── RF.joblib
    ├── XGB.joblib
    ├── LGBM.joblib
    ├── LR+CA125.joblib          fusion models (only when CA-125 data available)
    ├── SVM+CA125.joblib
    ├── RF+CA125.joblib
    ├── XGB+CA125.joblib
    └── LGBM+CA125.joblib

The fusion-model artefacts are skipped when no CA-125 data is available;
manifest.json records which models were trained.

Inference contract
------------------
A consumer (e.g. WebUI diagnosis page) loads:

    pre = joblib.load("model/preprocessing.joblib")
    clf = joblib.load("model/RF.joblib")
    # ...build feature vector x in the order pre['feature_names']
    x_imp  = np.where(np.isnan(x), pre['col_medians'], x)
    x_sc   = pre['scaler'].transform(x_imp.reshape(1, -1))
    p_mal  = clf.predict_proba(x_sc)[0, 1]

For fusion models, the CA-125 column is appended before scaling — see the
`apply_preprocessing()` helper below for the canonical recipe.
"""

import json
import os
import warnings
from typing import Optional

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

from utils.logger import logger
from scripts.modules.classifier_utils import (
    _load_feature_matrix,
    _select_features_in_fold,
    _impute_fit,
    _norm_id,
    _build_sid_to_case_map,
)


def _build_classifier_factories(ml_cfg: dict, spw: float) -> dict:
    """Fresh classifier instances with the same hyperparameters as Step 15."""
    return {
        "LR": LogisticRegression(
            C=ml_cfg.get("lr_C", 1.0), class_weight="balanced",
            max_iter=2000, random_state=42,
        ),
        "SVM": SVC(
            kernel="rbf", C=ml_cfg.get("svm_C", 1.0), probability=True,
            class_weight="balanced", random_state=42,
        ),
        "RF": RandomForestClassifier(
            n_estimators=ml_cfg.get("rf_n_estimators", 500),
            class_weight="balanced", random_state=42, n_jobs=-1,
        ),
        "XGB": XGBClassifier(
            n_estimators=ml_cfg.get("xgb_n_estimators", 200),
            max_depth=ml_cfg.get("xgb_max_depth", 4),
            learning_rate=ml_cfg.get("xgb_learning_rate", 0.1),
            scale_pos_weight=spw, use_label_encoder=False,
            eval_metric="logloss", random_state=42, n_jobs=-1, verbosity=0,
        ),
        "LGBM": LGBMClassifier(
            n_estimators=ml_cfg.get("lgbm_n_estimators", 200),
            max_depth=ml_cfg.get("lgbm_max_depth", 4),
            learning_rate=ml_cfg.get("lgbm_learning_rate", 0.1),
            scale_pos_weight=spw, random_state=42, n_jobs=-1, verbosity=-1,
        ),
    }


def run_train_final_models(workdir, result_dir, config,
                           benign_samples, malignant_samples):
    """
    Train production models on the full cohort (no CV holdout).

    The output layout is documented at the top of this module.
    """
    try:
        import joblib
    except ImportError as exc:
        logger.error(f"[final_model] joblib not installed; cannot persist models ({exc})")
        return None

    logger.info("━━ Step 15c: Train production models on full cohort ━━")

    # 1. Load feature matrix and apply outlier exclusion
    X, y, sample_ids, feature_names, region_in_88 = _load_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )
    if X.shape[0] == 0 or X.shape[1] == 0:
        logger.error("[final_model] Empty feature matrix — aborting")
        return None

    ml_cfg = (config or {}).get("ml", {})
    exclude_set = set(ml_cfg.get("exclude_samples", []))
    if exclude_set:
        keep_mask = np.array([sid not in exclude_set for sid in sample_ids])
        n_excluded = int((~keep_mask).sum())
        logger.info(f"  Excluding {n_excluded} sample(s) per exclude_samples")
        X = X[keep_mask]
        y = y[keep_mask]
        sample_ids = [s for s, k in zip(sample_ids, keep_mask) if k]

    n_neg = int((y == 0).sum())
    n_pos = int((y == 1).sum())
    spw   = n_neg / n_pos if n_pos > 0 else 1.0
    logger.info(f"  Training cohort: {len(y)} samples ({n_neg} benign, {n_pos} malignant), "
                f"{X.shape[1]} features → feature selection ...")

    # 2. Feature selection on the FULL cohort (no fold split)
    sel_idx = _select_features_in_fold(X, y, region_in_88, config)
    if len(sel_idx) == 0:
        logger.warning("  Feature selection returned 0 features; falling back to all features")
        sel_idx = np.arange(X.shape[1])
    selected_feature_names = [feature_names[i] for i in sel_idx]
    X_sel = X[:, sel_idx]
    logger.info(f"  Selected {len(sel_idx)} features for training")

    # 3. Impute (full-cohort medians) + standardize (fit on full cohort)
    X_imp, col_medians = _impute_fit(X_sel)
    scaler = StandardScaler()
    X_sc = scaler.fit_transform(X_imp)

    # 4. Optionally load CA-125 → fusion-model preprocessing
    include_ca125 = bool(ml_cfg.get("include_ca125_baseline", True))
    ca125_vals    = None
    ca125_median  = None
    ca_scaler     = None
    if include_ca125:
        try:
            from methylation.sample_grouping import load_ca125_map
            mode_cfg    = (config or {}).get("plasma", {})
            excel_sheet = mode_cfg.get("excel_sheet", "plasma")
            ca125_map   = load_ca125_map(category=excel_sheet)
            if ca125_map:
                sid_to_case = _build_sid_to_case_map(config)
                def _ca125_for(sid):
                    cid = sid_to_case.get(sid) or sid
                    return ca125_map.get(_norm_id(cid), np.nan)
                ca125_vals = np.array([_ca125_for(s) for s in sample_ids], dtype=float)
                n_with_ca = int((~np.isnan(ca125_vals)).sum())
                logger.info(f"  CA-125: {n_with_ca}/{len(ca125_vals)} samples have measurements")
                if n_with_ca == 0:
                    ca125_vals = None
        except Exception as exc:
            logger.warning(f"  CA-125 loading failed ({exc}); skipping fusion models")
            ca125_vals = None

    if ca125_vals is not None:
        ca125_median = float(np.nanmedian(ca125_vals)) if not np.all(np.isnan(ca125_vals)) else 0.0
        ca125_imp    = np.where(np.isnan(ca125_vals), ca125_median, ca125_vals)
        ca_scaler    = StandardScaler()
        ca125_sc     = ca_scaler.fit_transform(ca125_imp.reshape(-1, 1))
        X_aug_sc     = np.hstack([X_sc, ca125_sc])
    else:
        X_aug_sc = None

    # 5. Train classifiers on full cohort
    out_dir = os.path.join(result_dir, "model")
    os.makedirs(out_dir, exist_ok=True)

    # Save preprocessing first so a partial run still leaves a usable artefact set.
    preprocessing = {
        "feature_names":  selected_feature_names,
        "col_medians":    col_medians.tolist() if hasattr(col_medians, "tolist") else list(col_medians),
        "scaler":         scaler,
        "ca125_median":   ca125_median,
        "ca125_scaler":   ca_scaler,
    }
    joblib.dump(preprocessing, os.path.join(out_dir, "preprocessing.joblib"))
    logger.info(f"  Saved preprocessing → {out_dir}/preprocessing.joblib")

    classifier_factories = _build_classifier_factories(ml_cfg, spw)
    trained_models = []

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                message="X does not have valid feature names",
                                category=UserWarning)
        # cfDNA-only
        for name, clf in classifier_factories.items():
            clf.fit(X_sc, y)
            joblib.dump(clf, os.path.join(out_dir, f"{name}.joblib"))
            trained_models.append(name)
            logger.info(f"    ✓  {name}")

        # cfDNA + CA-125 fusion
        if X_aug_sc is not None:
            fusion_factories = _build_classifier_factories(ml_cfg, spw)
            for base_name, clf in fusion_factories.items():
                clf.fit(X_aug_sc, y)
                fname = f"{base_name}+CA125"
                joblib.dump(clf, os.path.join(out_dir, f"{fname}.joblib"))
                trained_models.append(fname)
                logger.info(f"    ✓  {fname}")

    # 6. Manifest — describes what was trained, for downstream consumers
    manifest = {
        "training_cohort": {
            "n_total":     len(y),
            "n_benign":    n_neg,
            "n_malignant": n_pos,
            "sample_ids":  sample_ids,
            "excluded":    sorted(exclude_set),
        },
        "panel": {
            "n_features": int(len(sel_idx)),
            "features":   selected_feature_names,
        },
        "ca125_available":  ca125_vals is not None,
        "ca125_threshold":  float(ml_cfg.get("ca125_threshold", 35.0)),
        "ca125_median":     ca125_median,
        "trained_models":   trained_models,
        "hyperparameters": {
            "lr_C":              ml_cfg.get("lr_C", 1.0),
            "svm_C":             ml_cfg.get("svm_C", 1.0),
            "rf_n_estimators":   ml_cfg.get("rf_n_estimators", 500),
            "xgb_n_estimators":  ml_cfg.get("xgb_n_estimators", 200),
            "xgb_max_depth":     ml_cfg.get("xgb_max_depth", 4),
            "xgb_learning_rate": ml_cfg.get("xgb_learning_rate", 0.1),
            "lgbm_n_estimators": ml_cfg.get("lgbm_n_estimators", 200),
            "lgbm_max_depth":    ml_cfg.get("lgbm_max_depth", 4),
            "lgbm_learning_rate":ml_cfg.get("lgbm_learning_rate", 0.1),
        },
        "feature_selection_mode": ml_cfg.get("feature_selection_mode", "union"),
        "max_features_clf":       ml_cfg.get("max_features_clf", 20),
    }
    manifest_path = os.path.join(out_dir, "manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False, default=str)
    logger.info(f"  Manifest → {manifest_path}")
    logger.info(f"━━ Trained {len(trained_models)} model(s) for inference: {trained_models}")
    return manifest


def apply_preprocessing(features_dict: dict, ca125: Optional[float],
                        preprocessing: dict, fusion: bool = False) -> np.ndarray:
    """
    Convenience helper for inference-time preprocessing.

    Parameters
    ----------
    features_dict : dict
        {feature_name: methylation_beta_value}. Missing keys are imputed with
        the training-fold median.
    ca125 : float or None
        Patient's serum CA-125 (U/mL). Required if `fusion=True`; ignored otherwise.
        Pass None to use the training-set median (treats unknown as average).
    preprocessing : dict
        Dictionary loaded from `model/preprocessing.joblib`.
    fusion : bool
        If True, append the (imputed, scaled) CA-125 value as the last column.
        Pair with a `*+CA125` classifier.

    Returns
    -------
    np.ndarray of shape (1, n_features) ready for `clf.predict_proba(...)`.
    """
    feat_names  = preprocessing["feature_names"]
    col_medians = np.asarray(preprocessing["col_medians"], dtype=float)

    x = np.array([
        features_dict.get(name, np.nan) for name in feat_names
    ], dtype=float)
    x_imp  = np.where(np.isnan(x), col_medians, x)
    x_sc   = preprocessing["scaler"].transform(x_imp.reshape(1, -1))

    if fusion:
        if preprocessing.get("ca125_scaler") is None:
            raise ValueError("fusion=True but no CA-125 scaler in preprocessing")
        ca = ca125 if (ca125 is not None and not np.isnan(ca125)) \
            else preprocessing["ca125_median"]
        ca_sc = preprocessing["ca125_scaler"].transform(np.array([[ca]], dtype=float))
        x_sc = np.hstack([x_sc, ca_sc])

    return x_sc
