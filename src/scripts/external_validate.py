"""external_validate — stub shim.

The real implementation lives in the private companion package
``aimethflow-core``, which contains patent-pending core algorithms and is
released only to journal-assigned peer reviewers during the pre-publication
phase.

Behaviour
---------
When ``aimethflow-core`` is installed in the same Python environment, importing
this module transparently substitutes ``aimethflow_core.external_validate``
for this module object in ``sys.modules``.  This shim is reachable both as a
plain module import and as a CLI script: invoking::

    python -m scripts.external_validate ...

executes the swapped-in implementation's ``main()`` after the substitution.

If ``aimethflow-core`` is not installed, this import raises ``ImportError``
with instructions for obtaining access.
"""
import importlib
import sys

try:
    _real = importlib.import_module("aimethflow_core.external_validate")
except ImportError as exc:  # pragma: no cover - environment-specific
    raise ImportError(
        "AIMethFlow's external_validate requires the protected core package "
        "`aimethflow-core`, which is private during the pre-publication phase.\n"
        "  - Peer reviewers: contact the corresponding author for access.\n"
        "  - See: https://github.com/tinglion/aimethflow#status"
    ) from exc

sys.modules[__name__] = _real

if __name__ == "__main__":
    _real.main()
