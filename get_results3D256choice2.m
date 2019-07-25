% Script to run the 3D simulations for the models (i) and (ii)
% with a perturbed initial MDF (choice = 2)

clear
n = 256;
N = 10000;
th_ang = 5;
seeds = 1:10;
angBrandon = 30;
list_t = [30 55];
list_t_original = [12 30];
choice = 2;
Rgrains = 5;
Rfamilies = 30;
run_simulations3D(n,N,th_ang,seeds,angBrandon,list_t,list_t_original,choice,Rgrains,Rfamilies)