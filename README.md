To the finished and published version: https://github.com/vi-rwth/FESTA.git \
To the publication: https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01022

Main differences:
- no PLUMED dependancy anymore:
  - if no FES-file provided: Histogram creation from COLVAR-file
  - columns read can be set manually for COLVAR- and FES-files
- new formats supported:
  - binary NPY-format for COLVAR-files (shape(N_frames, 2))
  - MLIPs CFG-format for trajectories (https://doi.org/10.1063/5.0155887) read+write
- increased accuracy:
  - no convex hull approximation anymore, if creation of single Polygon fails -> MultiPolygon
- increased speed:
  - pre-sorting of trajectory -> not all frames are evaluated
- small quality of life improvements:
  - multiple trajectory- and COLVAR-files can now be concatenated
  - FES-png now shows the true (Multi)Polygon outlines instead of raw selected frames
  - Preview Mode: only generate FES-png (no trajectory files are written) -> "--png only"

This is still highly work-in-progress and subject to change at any point
If you are using this in your work, cite the published version

Possible future additions:
- output of only unique structures (maybe as standalone script)
