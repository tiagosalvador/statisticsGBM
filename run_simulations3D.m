function run_simulations3D(n,N,th_ang,seeds,angBrandon,list_t,list_t_original,choice,Rgrains,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_simulations3D(n,N,th_ang,seeds,angBrandon,list_t,list_t_original,choice,Rgrains,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the simulations for 3D grain boundary motion with many grains
% using Algorithms 1 and 2 as described in the paper
%   Salvador, T. & Esedoglu, S. The Role of Surface Tension and
%   Mobility Model in Simulations of Grain Growth.
% Isotropic, unequal surface tensions by the Read-Shockley model are used
% by both algorithms. Algorithm 1 uses the reciprocal mobilities, while for
% Algorithm 2 all mobilities are equal to 1.
% INPUT:
%   n = Grid size.
%   N = Number of grains. (roughly we want floor(100*4^(log2(n)-7)))
%   th_ang = Threshold angle to merge grains.
%   seeds = list of seeds to use when generating the random initial data
%   angBrandon = Brandon angle in Read-Shockley model (degrees)
%   list_t = vector with times at which grain is saved.
%   list_t_original = vector with times at which grain is saved.
%   choice = variable to determine the type of perturbation to the initial
%   data. See the function 'initialvoronoidata3d' for implementation
%   details
%   Rgrains & Rfamilies are hyperparameters of the Algorithms. See the
%   functions 'gbm3d_sim' or 'gbm3doriginal_sim' for more details.
addpath(genpath(pwd)) % add all subfolders to the path

%% Setup

dims = [n n n]; % Dimensions of the grid: dims=[n n n] for nxnxn grid.
dt = 10/n^2; % Time step size.

mkdir('Simulations 3D')

fprintf('Grid resolution %d x %d x %d\n',dims)
for l = seeds
    clearvars -except n N th_ang seeds option choice dims dt l Rgrains Rfamilies angBrandon list_t list_t_original
    fprintf('Simulation #%d started\n',l)
    filename = ['./Simulations 3D/','simulation3d-',num2str(choice),'-',num2str(dims(1)),'x',num2str(dims(2)),'x',num2str(dims(3)),'-',num2str(l)];
    save(filename,'-v7.3');
    rng(l); % Control random number to be able to replicate results if needed
    [grains,ori] = initialvoronoidata3d(N,dims,Rgrains,Rfamilies,th_ang,choice);
    start_progress(' - Saving grains_t0')
    save_variable(filename,grains,'grains_t0')
    display_progress(1,1,1);
    start_progress(' - Saving ori_t0')
    save_variable(filename,ori,'ori_t0')
    save_variable(filename,length(list_t),'max_epoch')
    save_variable(filename,length(list_t_original),'max_epoch_original')
    display_progress(1,1,1);
    tic
    gbm3d_sim(dt,grains,dims,ori,Rgrains,Rfamilies,angBrandon,list_t,filename,1); 
    fprintf('Simulation #%d finished in %f seconds.\n',l,toc)
    load(filename,'grains_t0','ori_t0')
    tic
    gbm3doriginal_sim(dt,grains_t0,dims,ori_t0,Rgrains,Rfamilies,angBrandon,list_t_original,filename);
    fprintf('Simulation (original) #%d finished in %f seconds.\n',l,toc)
end % (for l).

end