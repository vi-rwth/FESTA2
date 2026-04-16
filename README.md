# FESTA2

> **Note:** This version is Work-in-Progress and subject to change at any point.

* **Published repository:** [vi-rwth/FESTA](https://github.com/vi-rwth/FESTA.git)
* **Publication:** [ACS JCIM (10.1021/acs.jcim.4c01022)](https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01022)

### Core Changes & Usage
* **No PLUMED dependency:** Histogram creation and reweighting are possible directly from a COLVAR file if no FES file is provided.
  > **CAUTION:** COLVAR and FES files are expected to have no headers.
* `--column <str>`: Sets the column index (1-indexed, not 0) for COLVAR or COLVAR+FES.
  > **CAUTION:** PLUMED usually prints the timestep in the first column of the COLVAR file, adjust `--column` accordingly.
* `--kBT <float>`: Sets the $k_B T$ value for reweighting. Requires specifying `BIAS` and `RCT` columns via `--column`.

### Supported Formats
* **NumPy (`.npy`):** Supported for COLVAR files. Shape must be `(N_frames, N_var)` where `N_var` is 2 (no reweighting) or 4 (reweighting).
* **CFG-format:** Read/write support added for trajectories.
* **Batch Input:** Multiple trajectory and COLVAR files (same format) are accepted via bash globbing.

### New Features & Flags
* `--png only`: Generates a preview of Polygons+FES without frame separation.
* `--stride 0`: Outputs only the single most representative frame of each minimum.
* `--thresh_low <float>`: Omits areas with energies lower than `<float>` (useful for transition state analysis).
* `--pbc`: Standard flag indicating `pbc=True`.

### Accuracy & Performance
* **Exact Polygons:** Removed convex hull approximation. Uses `MultiPolygon` fallback if single Polygon creation fails. `FES-png` now displays true outlines.
* **Multi-core Processing:** Polygon distance calculations and frame extraction/printing (`.cfg`, `.pdb`, `.lammpstrj`) execute in chunks utilizing all CPU cores.
* **Spatial Hashing:** Optimized algorithm skips unnecessary frame evaluations for faster processing.
* **Robustness:** FES files now accept non-finite energy values (e.g., NaNs), and automatic threshold determination works consistently.
