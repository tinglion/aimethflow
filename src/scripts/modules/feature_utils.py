"""feature_utils — stub shim.

The real implementation lives in the private companion package
``aimethflow-core``, which contains patent-pending core algorithms and is
released only to journal-assigned peer reviewers during the pre-publication
phase.

Behaviour
---------
When ``aimethflow-core`` is installed in the same Python environment, importing
this module transparently substitutes ``aimethflow_core.feature_utils`` for
this module object in ``sys.modules``.  All existing call sites of the form::

    from scripts.modules.feature_utils import (
        build_feature_matrix, run_lasso_selection, ...,
    )

continue to work unchanged, because Python resolves the attribute on the
substituted module after the swap.

If ``aimethflow-core`` is not installed, this import raises ``ImportError``
with instructions for obtaining access.
"""
import importlib
import sys

try:
    _real = importlib.import_module("aimethflow_core.feature_utils")
except ImportError as exc:  # pragma: no cover - environment-specific
    raise ImportError(
        "AIMethFlow's feature_utils requires the protected core package "
        "`aimethflow-core`, which is private during the pre-publication phase.\n"
        "  - Peer reviewers: contact the corresponding author for access.\n"
        "  - See: https://github.com/tinglion/aimethflow#status"
    ) from exc

sys.modules[__name__] = _real
