#!/usr/bin/env bash
# Set up RetroSynth: Python venv + packages, claude CLI, and an ORCA check.
# Usage: ./install.sh [venv_dir]   (default venv_dir: retro)
set -euo pipefail

cd "$(dirname "$0")"
VENV="${1:-retro}"
PYTHON="${PYTHON:-python3}"

echo "==> Python virtualenv: $VENV"
if [ ! -x "$VENV/bin/python" ]; then
    "$PYTHON" -m venv "$VENV"
fi
PY="$VENV/bin/python"
"$PY" -m pip install --upgrade pip

echo "==> Python packages"
if ! "$PY" -m pip install -r requirements.txt; then
    echo "!! Full install failed, retrying without tblite (it is often only on conda-forge)"
    grep -v '^tblite' requirements.txt > .requirements.no-tblite.txt
    "$PY" -m pip install -r .requirements.no-tblite.txt
    rm -f .requirements.no-tblite.txt
    if command -v conda >/dev/null 2>&1; then
        echo "   conda found; install tblite into a conda env with:"
        echo "   conda install -c conda-forge tblite-python"
    else
        echo "   tblite missing: the gfn2 backend will be unavailable."
    fi
fi

echo "==> dxtb (without dependency resolution, see requirements-dxtb.txt)"
"$PY" -m pip install --no-deps -r requirements-dxtb.txt

echo "==> Import check"
"$PY" - <<'EOF'
import importlib
missing = []
for module in ["dspy", "rdkit", "numpy", "scipy", "pandas", "PIL", "rich", "dotenv", "tblite", "torch", "dxtb", "pytest"]:
    try:
        importlib.import_module(module)
    except Exception as error:
        missing.append(f"{module} ({error})")
print("   all modules import" if not missing else "!! failed imports: " + ", ".join(missing))
EOF

echo "==> claude CLI (needed for model_type = subscription)"
if command -v claude >/dev/null 2>&1; then
    echo "   found: $(command -v claude)"
elif command -v npm >/dev/null 2>&1; then
    npm install -g @anthropic-ai/claude-code
    echo "   installed; run 'claude' once to log in."
else
    echo "!! npm not found; install Node.js, then: npm install -g @anthropic-ai/claude-code"
fi

echo "==> ORCA (optional QC backend)"
if command -v orca >/dev/null 2>&1; then
    echo "   found: $(command -v orca)"
else
    echo "   not found; download it manually from https://orcaforum.kofo.mpg.de (registration required)"
    echo "   and set orca_path in config.ini if it is not on PATH."
fi

echo "==> Done. Activate with: source $VENV/bin/activate"
