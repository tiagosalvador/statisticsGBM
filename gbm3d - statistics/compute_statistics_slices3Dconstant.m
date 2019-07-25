function [areas_grains_slices, areas_convhull_slices, nNeighbors_slices,...
    isopratio_slices] = ...
    compute_statistics_slices3Dconstant(grains,dims,Rfamilies,dt,spacing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [areas_grains_slices, areas_convhull_slices, nNeighbors_slices,...
%    isopratio_slices] = ...
%    compute_statistics_slices3Dconstant(grains,dims,Rfamilies,dt,spacing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes different grain statistics for the 2d slices of 3D grains
% for model (iii).

N = size(grains,1);

%% Setting up global variables for slices

dims_slices = dims([1 2]); % Assuming that dims = [n n n]

% Auxiliary vars.
global WORKSPACE KERNEL Z;

m = dims_slices(1); % Size of computational grid.
n = dims_slices(2); % Size of computational grid.

WORKSPACE = int32(-ones(dims_slices)); % Needed for pgrow.c.
Z = -ones(m,n); % Another Workspace var.
                % Needed in loc_levset_to_volfluid.c.

global work_x work_y work_bx work_by work_candidate_x work_candidate_y; % Needed for pgrow3.c.
work_x = int32(zeros(prod(dims_slices),1));
work_y = int32(zeros(prod(dims_slices),1));
work_bx = int32(zeros(prod(dims_slices),1));
work_by = int32(zeros(prod(dims_slices),1));
work_candidate_x = int32(zeros(prod(dims_slices),1));
work_candidate_y = int32(zeros(prod(dims_slices),1));

% Prepare the convolution kernel KERNEL in fourier space:
I = sqrt(-1); % Imaginary number i.
w=exp(I*2*pi/n); % nth root of unity.
[x,y] = meshgrid(1:n,1:n); x = x'; y = y';
KERNEL = exp(-dt*n*n*(4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y)));


slices = 1:spacing:dims(1);

%% Computing statistics slices

areas_slices = cell(length(dims),length(slices));
nNeighbors_slices = cell(length(dims),length(slices));
isopratio_conv_slices = cell(length(dims),length(slices));
isopratio_pixels_slices = cell(length(dims),length(slices));
lengths_conv_slices = cell(length(dims),length(slices));
lengths_pixels_slices = cell(length(dims),length(slices));

for i = 1:length(dims)
    switch i
        case 1
            dims_slices = dims([2 3]);
        case 2
            dims_slices = dims([1 3]);
        case 3
            dims_slices = dims([1 2]);
    end
    start_progress([' - Computing section statistics in dim = ',num2str(i)])
    for s = 1:length(slices)
        [grains_slices, ID_slices] = compute_grains_slices(grains,dims,dims_slices,slices(s),i,Rfamilies);
        [areas_grains_slices{i,s},areas_convhull_slices{i,s},...
            nNeighbors_slices{i,s},...
            isopratio_slices{i,s}] = compute_statistics_slices3D_auxconstant(grains_slices,ID_slices,dims_slices,dt);
        display_progress(s,length(slices),1);
    end
end
end
