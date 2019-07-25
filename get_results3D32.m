% Script to run the 3D simulations for the models (i) and (ii)
% on a 32x32x32 3D grid.

clear
n = 32;
N = 100;
th_ang = 5;
seeds = 1:3;
angBrandon = 30;
list_t = [1 2];
list_t_original = [1 2];
choice = 1;
Rgrains = 5;
Rfamilies = 30;
run_simulations3D(n,N,th_ang,seeds,angBrandon,list_t,list_t_original,choice,Rgrains,Rfamilies)