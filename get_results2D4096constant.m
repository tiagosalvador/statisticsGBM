% Script to run the 2D simulations for model (iii)

clear
n = 4096;
N = 100000;
seeds = 1:3;
list_t = [30 118];
run_simulations2Dconstant(n,N,seeds,list_t)