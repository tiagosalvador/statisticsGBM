% This script compiles the mex files needed
% to run the 3D algorithm

addpath(genpath(pwd))
cd 'gbm3d - c_files'/
mex dilation.c
mex get_nhd_grains3d.c
mex get_nhd_grains3doriginal.c
mex ls2vf3d.c
mex updatelevelsetdata3d.c
mex updatelevelsetdata3doriginal.c
mex updatelevelsetdata3dconstant.c
cd ..