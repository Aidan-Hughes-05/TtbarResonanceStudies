#!/usr/bin/env bash
set -euo pipefail

MG5_DIR="$PWD/MG5_aMC_v3_5_11"

echo "Checking Conda..."
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda is not available in PATH."
    exit 1
fi

if [[ ! -x "$MG5_DIR/bin/mg5_aMC" ]]; then
    echo "ERROR: MadGraph not found at $MG5_DIR. Run this from TtbarResonanceStudies."
    exit 1
fi

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"

echo "Creating or updating the mc_gen environment with ROOT and Pythia8..."
if conda env list | awk '$1 == "mc_gen" { found=1 } END { exit !found }'; then
    conda install -y -n mc_gen -c conda-forge root pythia8
else
    conda create -y -n mc_gen -c conda-forge root pythia8
fi
conda activate mc_gen

echo "Installing the Pythia8 interface and Delphes into MadGraph..."
commands_file="$(mktemp)"
trap 'rm -f "$commands_file"' EXIT
printf 'set pythia8_path %s \ninstall mg5amc_py8_interface\ninstall Delphes\n' \
    "$CONDA_PREFIX" > "$commands_file"
"$MG5_DIR/bin/mg5_aMC" "$commands_file"

echo "Checking ROOT..."
python -c 'import ROOT; print("ROOT", ROOT.gROOT.GetVersion())'

echo "Checking Pythia8..."
pythia8-config --version

echo "Done. Conda Pythia8, its MadGraph interface, and Delphes are ready."