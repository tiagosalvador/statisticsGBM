% Script to run the 2D simulations for the models (i) and (ii)

clear
n = 4096;
N = 100000;
th_ang = 1;
seeds = 1:3;
angBrandon = 30;
list_t = [32 140];
list_t_original = [26 105];
run_simulations2D(n,N,th_ang,seeds,angBrandon,list_t,list_t_original)