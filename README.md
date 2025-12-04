To the finished and published version: https://github.com/vi-rwth/FESTA.git \
To the publication: https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01022

Main differences:
- no PLUMED dependancy anymore:
  - if no FES-file provided: Histogram creation from COLVAR-file
  - read columns manually setable for COLVAR- and FES-files
- Pre-sorting of trajectory and therefore significant speed increase
- Compatible with MLIPs CFG-format (https://doi.org/10.1063/5.0155887) read+write
- Multiple trajectory- and COLVAR-files can be concatenated

This is purely experimental and highly work-in-progress

possible future additions:
- output of only unique structures
