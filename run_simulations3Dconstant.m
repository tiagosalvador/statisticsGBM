function run_simulations3Dconstant(n,N,seeds,list_t,Rgrains,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_simulations3Dconstant(n,N,seeds,list_t,Rgrains,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the simulations for 3D grain boundary motion with many grains
% using Algorithms 1 as described in the paper
%   Salvador, T. & Esedoglu, S. The Role of Surface Tension and
%   Mobility Model in Simulations of Grain Growth.
% Isotropic, equal to 1 surface tensions and mobilities are used.
% INPUT:
%   n = Grid size.
%   N = Number of grains. (roughly we want floor(100*4^(log2(n)-7)))
%   seeds = list of seeds to use when generating the random initial data
%   list_t = vector with times at which grain is saved.
%   Rgrains & Rfamilies are hyperparameters of the Algorithms. See the
%   functions 'gbm3d_sim' or 'gbm3doriginal_sim' for more details.


addpath(genpath(pwd)) % add all subfolders to the path

%% Setup

dims = [n n n]; % Dimensions of the grid: dims=[n n n] for nxnxn grid.
dt = 10/n^2; % Time step size.

mkdir('Simulations 3D')

fprintf('Grid resolution %d x %d x %d\n',dims)
for l = seeds
    clearvars -except n N th_ang seeds option dims dt l Rgrains Rfamilies list_t
    fprintf('Simulation #%d started\n',l)
    filename = ['./Simulations 3D/','simulation3dconstant-',num2str(dims(1)),'x',num2str(dims(2)),'x',num2str(dims(3)),'-',num2str(l)];
    save(filename,'-v7.3');
    rng(l); % Control random number to be able to replicate results if needed
    grains = initialvoronoidata3dconstant(N,dims,Rgrains);
    start_progress(' - Saving grains_t0_constant')
    save_variable(filename,grains,'grains_t0_constant')
    display_progress(1,1,1);
    tic
    gbm3dconstant_sim(dt,grains,dims,Rgrains,Rfamilies,list_t,filename);
    fprintf('Simulation #%d finished in %f seconds.\n',l,toc)
end % (for l).

end