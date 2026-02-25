To the published version: https://github.com/vi-rwth/FESTA.git \
To the publication: https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01022

Main differences:
- no PLUMED dependency anymore:
  - if no FES-file provided: histogram creation from COLVAR-file
  - "--column": sets the column index (starting at 1) for COLVAR, FES or both. 
  - CAUTION: This means that COLVAR is not expected to have any headers!
  - CAUTION: Usually PLUMED prints the timestep in the first column in COLVAR, adjust --column accordingly!
- new formats supported:
  - binary NPY-format for COLVAR-files (shape(N_frames, 2))
  - MLIPs CFG-format for trajectories (https://doi.org/10.1063/5.0155887) read+write
- increased accuracy:
  - no convex hull approximation anymore, if creation of single Polygon fails -> MultiPolygon
- increased performance:
  - spatial hashing algorithm -> not all frames are evaluated
  - polygon distance calculation in chunks using all cores
  - frame extraction/printing for CustomWriter in chunks using all cores
  - various other minor performance improvements
- quality of life improvements:
  - multiple trajectory- and COLVAR-files (same format) are now accepted (bash glob supported)
  - FES-png now shows the true (Multi)Polygon outlines instead of raw selected frames
  - "--png only": only generate preview of Polygons+FES (no frame separation)
  - "--stride 0": only output the single most representative frame of each minimum

This is still highly work-in-progress and subject to change at any point \
If you are using this in your work, cite the published version above
