# Ttbar resonance
This repository contains the code needed to generate events for the ttbar/heavy-higgs interference.
These events are used to train algorithms that regress the parameters governing the BSM models.

## Repository structure

This repository is organized into event generation, conversion, analysis, and ML workflows:

- `Generate/` – MadGraph job options and HTCondor submission scripts for `SIGNAL`, `INTERFERENCE`, and `BG` samples.
- `FileConverters/` – converters from LHE to ROOT, and ROOT to H5.
- `MadGraph_Delphes_Analysis/` – C++ analysis framework to analyse ttbar samples and generate NTuples for training.
- `MLNotebooks/` – notebooks and trained PyTorch checkpoints for parameter/mass estimation.
- `PeakDeepMaster/` – Git submodule containing the PyTorch Lightning training/inference framework.
- `Custom_HELAS/` – custom HELAS Fortran routines used by matrix-element generation.
- `Outputs/` – produced ROOT/HDF5 outputs.
- `MG5_aMC_v3_5_11/` – local MadGraph installation used for event production.
- `install/` and `build/` – compiled analysis binaries, headers, and build artifacts.

---
> [!IMPORTANT]
> **Note:** The steps below are **not sequential**. Each section describes a different part of the workflow and can be followed independently.

## Clone with submodules

Clone the repository together with `PeakDeepMaster`:

```bash
git clone --recurse-submodules <repo-url>
```

If you already cloned the repository without submodules, run:

```bash
git submodule update --init --recursive
```


## 1) Compile, install, and use the C++ analysis framework

The C++ framework is in `MadGraph_Delphes_Analysis/` and builds an executable called `analysis_app`.

### Step 1.1: Prerequisites

Just run in `lxplus` and setup `StatAnalysis` release. This was successfully compiled with:
```bash
setupATLAS --quiet && lsetup git && asetup StatAnalysis,0.7.2
```

### Step 1.2: Configure + build + install

From repository root:

```bash
cmake -S MadGraph_Delphes_Analysis -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j10 --target install
```

### Step 1.3: Set runtime environment

```bash
source install/setup.sh
```

This adds `install/bin` to `PATH` and `install/lib` to `LD_LIBRARY_PATH`.

### Step 1.4: Run the analysis executable

`analysis_app` usage:

```bash
analysis_app <input_file.root OR input_dir> <output_file.root> <mode> [comma_separated_filter]
```

- `mode=h` -> histogram mode
- `mode=n` -> ntuple mode (used for downstream ML conversion)
- If `input_dir` is passed, optional `comma_separated_filter` selects matching `.root` files by name.

Examples:

```bash
# Build ntuples from all ROOT files in a directory
analysis_app /path/to/root_inputs /path/to/TtbarNtuples.root n

# Build ntuples from files that match "Signal" and "m800"
analysis_app /path/to/root_inputs /path/to/Signal_m800.root n Signal,m800

# Build histograms from a single ROOT file
analysis_app /path/to/sample.root /path/to/hists.root h
```

---

## 2) Generate samples with HTCondor

Generation scripts are in `Generate/`.

### Step 2.1: Go to submission directory

```bash
cd Generate
```

### Step 2.2: Submit one job manually

```bash
condor_submit -append "ARGS=Signal_c1p5_m600 SIGNAL 1.5 600" Submit.sub
```

`ARGS` format:

```text
<tag> <PROCESS> <COUPLING> <MASS>
```

- `PROCESS`: `SIGNAL` or `INTERFERENCE`
- `tag` convention: `Signal_cXpY_mZZZ` (example: `Signal_c1p5_m600`)

### Step 2.3: Submit a grid of jobs (recommended)

Edit arrays in `Generate/submit_condor_jobs.sh`:

- `PROCESSES`
- `COUPLINGS`
- `MASSES`

Then submit all combinations:

```bash
bash submit_condor_jobs.sh
```

### Step 2.4: Check outputs

- Condor logs are configured in `Submit.sub` (`log`, `output`, `error`).
- Generated LHE files are written by `executable.sh` under your configured MadGraph output path, grouped by tag.

---

## 3) Convert files: LHE -> ROOT -> H5

Conversion utilities are in `FileConverters/`.

### Step 3.1: LHE -> ROOT (single class/directory)

```bash
cd FileConverters
python3 LHEtoROOT.py \
	--input-dir /path/to/lhe_dir \
	--file-class Signal_c1p5_m600 \
	--output-dir /path/to/root_out \
	--n-processes 8 \
	--no-interactive
```

This produces ROOT files named like `Signal_c1p5_m600_<index>.root`.

### Step 3.2: LHE -> ROOT (batch over many subdirectories)

```bash
python3 LHEtoROOT_Batch.py \
	--input-dir /path/to/lhe_root_dir \
	--output-dir /path/to/root_out
```

Each subdirectory inside `--input-dir` is treated as one class/tag.

### Step 3.3: ROOT -> analysis ntuple ROOT (optional but typical)

If your LHE-derived ROOT files contain `mytree`, convert them to the analysis ntuple format expected by the H5 converter (`Ttbar` tree):

```bash
analysis_app /path/to/root_out /path/to/TtbarNtuples.root n
```

### Step 3.4: ROOT -> H5

```bash
python3 H5_Converter/Ttbar_Events.py \
	-i /path/to/TtbarNtuples.root,/path/to/AnotherNtuple.root \
	-o /path/to/Interference_vs_other.h5
```

Notes:

- Input files must contain a `Ttbar` tree.
- Labels are assigned automatically from filename:
	- files containing `Signal` or `Interference` -> class `1`
	- other files -> class `0`

---

## Legacy generation note

Original generation notes are also available here:
[MadGraph generation note](https://www.notion.so/MadGraph-generation-286570f8119b80c4a6a8d9cec861f13c?source=copy_link)

