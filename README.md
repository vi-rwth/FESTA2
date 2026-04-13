To the published version: https://github.com/vi-rwth/FESTA.git \
To the publication: https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01022

Main differences:
- no PLUMED dependency anymore:
  - if no FES-file provided: histogram creation from COLVAR-file
  - "--column": sets the column index (first index is 1, not 0) for COLVAR, FES or both. 
  - CAUTION: This means that COLVAR- and FES-files are not expected to have any headers!
  - CAUTION: Usually PLUMED prints the timestep in the first column in COLVAR, adjust --column accordingly!
- new formats supported:
  - NumPy NPY-format for COLVAR-files (shape(N_frames, 2))
  - CFG-format for trajectories (https://doi.org/10.1063/5.0155887) read+write
- new features:
  - "--png only": only generate preview of Polygons+FES (no frame separation).
  - "--stride 0": only outputs the single most representative frame of each minimum.
  - "--thresh_low X": omit all areas with energies lower than X. Useful for transition state analysis.
- increased accuracy:
  - no convex hull approximation anymore, if creation of single Polygon fails -> MultiPolygon
- increased performance:
  - spatial hashing algorithm -> not all frames are evaluated
  - polygon distance calculation in chunks using all cores
  - frame extraction/printing for CustomWriter in chunks using all cores (cfg, pdb, lammpstrj)
  - various other minor performance improvements
- quality of life improvements:
  - multiple trajectory- and COLVAR-files (same format) are now accepted (bash glob supported)
  - FES-png now shows the true (Multi)Polygon outlines instead of raw selected frames
  - "--pbc" now denotes pbc=True, no flag pbc=False
  - FES-file is now allowed to have non-finite energy values (e.g. NaNs)

This is still work-in-progress and subject to change at any point \
If you are using this in your work, please cite the published version above
