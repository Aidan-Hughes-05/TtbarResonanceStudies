# Generation Instructions

- Prerequisites: HTCondor (`condor_submit`) available on PATH.
- Work from this folder so `Submit.sub` is found.

## Single Job
Submit one job by specifying tag, process, coupling, and mass:

```bash
condor_submit -append "ARGS=Signal_c1p5_m600 SIGNAL 1.5 600" Submit.sub
```

- Tag format: `<Signal|Interference>_c<coupling-with-p>_m<mass>` (e.g., `Signal_c0p7_m900`).
- Coupling: non-negative (e.g., `0.1, 0.3, 0.5, 0.7, 0.8, 1.0, 1.5`).
- Mass: one of `600, 700, 800, 900, 1000`.

## Batch
Recommended: use the central submit script:

```bash
# Edit lists in ../submit_condor_jobs.sh if needed, then run:
../submit_condor_jobs.sh
```
