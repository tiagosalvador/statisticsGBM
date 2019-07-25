% Script to run 3D simulations for model (iii)

clear
n = 256;
N = 10000;
seeds = 1:10;
list_t = [12 30];
Rgrains = 5;
Rfamilies = 30;
run_simulations3Dconstant(n,N,seeds,list_t,Rgrains,Rfamilies)