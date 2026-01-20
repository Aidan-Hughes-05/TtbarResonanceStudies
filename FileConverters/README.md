# File Conversion (LHE → ROOT)

Convert MadGraph LHE files to ROOT.

- Inputs: LHE files in a directory (or tree of directories).
- Outputs: ROOT files named `<file-class>_<index>.root` in the chosen output directory.
- Requirements: Python 3, environment capable of running ROOT conversion (see LHE_Converter).

## Single Directory
Use the converter on a directory of LHE files that belong to the same class.

```bash
python3 LHEtoROOT.py \
  --input-dir /path/to/lhe_dir \
  --file-class Signal_c1p5_m600 \
  --output-dir /path/to/root_out \
  --n-processes 8 \
  --no-interactive
```

- `--input-dir`: directory containing `.lhe` files (searched recursively).
- `--file-class`: prefix added to output ROOT file names.
- `--output-dir`: where ROOT files will be written (created if missing).
- `--n-processes`: parallel workers (default 4).
- `--no-interactive`: skip confirmation prompt.

## Batch (multiple subdirectories)
Run over each subdirectory inside an input directory. Each subdirectory name becomes the `file-class`.

```bash
python3 LHEtoROOT_Batch.py \
  --input-dir /path/to/lhe_root_dir \
  --output-dir /path/to/root_out
```

- Uses `LHEtoROOT.py` under the hood with `--n-processes 10` and `--no-interactive`.
- Example: if `/path/to/lhe_root_dir` has `Signal_c1p5_m600/`, `Interference_c0p7_m900/`, outputs will be named accordingly in `/path/to/root_out`.

## Files
- [LHEtoROOT.py](LHEtoROOT.py): convert one class/directory of LHE files.
- [LHEtoROOT_Batch.py](LHEtoROOT_Batch.py): iterate over subdirectories and call the single-directory converter.
- [LHE_Converter/](LHE_Converter/): underlying LHE→ROOT conversion utilities.

## H5 Conversion (ROOT → HDF5)
Combine one or more ROOT files into a single HDF5 dataset.

```bash
python3 H5_Converter/Ttbar_Events.py \
  -i /path/to/Signal.root,/path/to/BG.root \
  -o /path/to/out.h5
```

- Inputs: comma-separated list of ROOT files containing a `Ttbar` tree.
- Output: a single `.h5` file with groups `INPUTS`, `LABELS`, `METADATA`.
- Labels: files with names containing `Signal` or `Interference` are labeled class `1`, others `0`.
- Optional: `-p/--print-branches` to print ROOT branches and exit.

Tips:
- Ensure ROOT files have branches required by the converter (e.g., `Top_1_*`, `Top_2_*`, `lumiWeight`, `mttbar`, `Coupling`).
- The converter streams in chunks; no extra flags needed for large files.
