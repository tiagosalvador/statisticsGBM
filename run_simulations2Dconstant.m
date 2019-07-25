function run_simulations2Dconstant(n,N,seeds,list_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_simulations2Dconstant(n,N,seeds,list_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the simulations for 2D grain boundary motion with many grains
% using Algorithms 1 as described in the paper
%   Salvador, T. & Esedoglu, S. The Role of Surface Tension and
%   Mobility Model in Simulations of Grain Growth.
% Isotropic, equal to 1 surface tensions and mobilities are used.
% INPUT:
%   n = Grid size.
%   N = Number of grains. (roughly we want floor(100*4^(log2(n)-7)))
%   seeds = list of seeds to use when generating the random initial data
%   list_t = vector with times at which grain is saved.

addpath(genpath(pwd)) % add all subfolders to the path

%% Setup
dims = [n n]; % Dimensions of the grid: dims=[n n] for nxn grid.
dt = 10/n^2; % Time step size.
Rgrains = 10;
Rfamilies = 30;

mkdir('Simulations 2D')

fprintf('Grid resolution %d x %d\n',dims)
for l = seeds
    clearvars -except n N seeds dims dt l Rgrains Rfamilies angBrandon list_t
    fprintf('Simulation #%d started\n',l)
    filename = ['./Simulations 2D/','simulation2dconstant-',num2str(dims(1)),'x',num2str(dims(1)),'-',num2str(l)];
	save(filename,'-v7.3');
    rng(l); % Control random number to be able to replicate results if needed
    grains = initialvoronoidata2dconstant(N,dims,Rgrains);
    start_progress(' - Saving grains_t0_constant')
    save_variable(filename,grains,'grains_t0_constant')
    display_progress(1,1,1);
    tic
    gbm2dconstant_sim(dt,grains,dims,Rgrains,Rfamilies,list_t,filename);
    fprintf('Simulation #%d finished in %f seconds.\n',l,toc)
end % (for l).
end