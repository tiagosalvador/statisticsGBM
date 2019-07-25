% This script compiles the mex files needed
% to run the 2D algorithm

addpath(genpath(pwd))
cd 'gbm2d - c_files'/
mex get_nhd_grains2d.c
mex get_nhd_grains2doriginal.c
mex updatelevelsetdata2d.c
mex updatelevelsetdata2doriginal.c
mex updatelevelsetdata2dconstant.c
mex pgrow3.c
mex loc_levset_to_volfluid.c
cd ..