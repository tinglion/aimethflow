#!/usr/bin/env bash
# docker-entrypoint.sh — auto-install protected core if a wheel is mounted.
#
# This container ships with stub modules for Steps 14, 15, 19. To run the full
# pipeline, drop the `aimethflow_core-<version>-py310-none-any.whl` file
# (received from the corresponding author during the iMeta peer-review phase)
# into the `./docs/` directory on the host. docker-compose mounts it at
# `/docs/` inside the container (read-only), and this entrypoint installs it
# on startup.
#
# Without the wheel, Steps 1-13, 16-18, and 20 still run end-to-end; attempting
# Steps 14 / 15 / 19 raises an ImportError with the access path.

set -e

PY=/opt/conda/envs/meth/bin/python
PIP=/opt/conda/envs/meth/bin/pip

if compgen -G "/docs/aimethflow_core-*.whl" > /dev/null 2>&1; then
    WHEEL=$(ls /docs/aimethflow_core-*.whl | head -1)
    if "$PY" -c "import aimethflow_core" 2>/dev/null; then
        echo "[entrypoint] aimethflow_core already installed; skipping wheel install"
    else
        echo "[entrypoint] Installing protected core from mounted wheel: $WHEEL"
        "$PIP" install --quiet --no-deps "$WHEEL"
        echo "[entrypoint] OK: Steps 14, 15, 19 now available"
    fi
else
    echo "[entrypoint] No /docs/aimethflow_core-*.whl found."
    echo "[entrypoint] Steps 1-13, 16-18, 20 will run. Steps 14, 15, 19 will report ImportError."
    echo "[entrypoint] To enable the full pipeline: drop the wheel into ./docs/ on the host and restart."
fi

exec "$@"
