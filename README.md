# Radial-Field-of-Views MATLAB package
Radial sampling patterns in frequency space to support anisotropic field of views (FOVs) in the dual domain


This archive includes Matlab functions to design acquisition schemes
for radial fourier/frequency space sampling that will support the desired field-of-view (FOV), where the primary application is for MRI acquisition trajectories in k-space.
Included are design functions for 2D and 3D radial frequency space sampling, as well as 
PROPELLER MRI trajectories.  C implementations for designing
elliptical FOV trajectories are also included.

A recent addition is a 2D golden angle ordering (with anisotropic FOV) and a 3D golden angle ordering (only isotropic FOV).  There is currently a need for a
- 3D anisotropic FOV golden angle strategies


Within Matlab, use "help Radial-Field-of-Views" for information on the design
functions included in this package.  Documentation on each specific
function is also included.

This package accompanies the journal article and conference abstract:
Peder E.Z. Larson and Paul T. Gurney and Dwight G. Nishimura.
"Anisotropic Field-of-Views in Radial Imaging."
*IEEE Transactions on Medical Imaging.* 27(1): 47-57 (2008).
* http://dx.doi.org/10.1109/TMI.2007.902799
* https://arxiv.org/abs/2101.04660
* https://radiology.ucsf.edu/sites/radiology.ucsf.edu/files/import/filemanager/research/Larson/LarsonRadialFOVs.pdf

Peder E. Z. Larson, M. Lustig, and Dwight G. Nishimura. “Anisotropic field-of-view shapes for improved PROPELLER imaging.” *Magn Reson Imaging*, 27(4):470–479, May
2009.
* http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2882965/
* https://radiology.ucsf.edu/sites/radiology.ucsf.edu/files/import/filemanager/research/Larson/LarsonPropellerFOV.pdf

Contributors: Peder Larson, Paul Gurney, Wenwen Jiang
(c) 2007-2015, Board of Trustees, Leland Stanford Junior University and The Regents of the University of California. 

