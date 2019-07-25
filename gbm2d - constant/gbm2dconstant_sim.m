function gbm2dconstant_sim(dt,grains,dims,Rgrains,Rfamilies,list_t,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gbm2dconstant_sim(dt,grains,dims,Rgrains,Rfamilies,list_t,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D grain boundary motion with many grains. Isotropic, equal
% surface tensions. All mobilities = 1.
%
% INPUT:
%   dt = time step size
%   grains = Data structure containing initial configuration of
%            grains. Cell array, with grains{k,1} = indicies of
%            pix near grain k, and grains{k,2} = level set vals
%            at those locations. Eventually, grains{k,3} and
%            grains{k,4} contains convolution values at same
%            locations.
%   dims = Dimensions of the grid: dims=[n n] for nxn grid.
%   Rgrains = radius by each a grain is enlarged.
%   Rfamilies = Minimum distance to maintain between grains in the same 
%               family.
%   angBrandon = Brandon angle in Read-Shockley model (degrees)
%   list_t = vector with times at which grain is saved.
%   filename = name of the file where intermidiate grains are stored.
%              Required because the C routines modify matlab data in place.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAUTION:
%   Some C routines modify matlab data structures in place.
%   This is done in the interest of minimal memory utilization.
%   It has the undesirable side effect that input arrays cannot
%   be guaranteed to remain unaltered, especially 'grains', even
%   when the 1st output var ~= 3rd input var.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE:
%   This program is based on Algorithm 1 in:
%
%   Salvador, T. & Esedoglu, S. The Role of Surface Tension and
%   Mobility Model in Simulations of Grain Growth.
%
%   which was originally proposed in:
%
%   Esedoglu, S.; Otto, F. Threshold dynamics for networks with
%   arbitrary surface tensions. Communications on Pure and
%   Applied Mathematics. 68:5 (2015), pp. 808-864.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%   1. The order (and number) of grains listed in data structure
%      "grains" may change from one iteration to the next. But,
%      orientations in "ori" follow the grains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Auxiliary vars.
global WORKSPACE KERNEL Z ID;

m = dims(1); % Size of computational grid.
n = dims(2); % Size of computational grid.

WORKSPACE = int32(-ones(dims)); % Needed for pgrow.c.
Z = -ones(m,n); % Another Workspace var.
                % Needed in loc_levset_to_volfluid.c.

ID = (1:1:size(grains,1))';

global work_x work_y work_bx work_by work_candidate_x work_candidate_y; % Needed for pgrow3.c.
work_x = int32(zeros(prod(dims),1));
work_y = int32(zeros(prod(dims),1));
work_bx = int32(zeros(prod(dims),1));
work_by = int32(zeros(prod(dims),1));
work_candidate_x = int32(zeros(prod(dims),1));
work_candidate_y = int32(zeros(prod(dims),1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the convolution kernel KERNEL in fourier space:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = sqrt(-1); % Imaginary number i.
w=exp(I*2*pi/n); % nth root of unity.
[x,y] = meshgrid(1:n,1:n); x = x'; y = y';
KERNEL = exp( -dt*n*n*( 4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y) ) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN TIME LOOP STARTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(grains,1);
epoch = 0;
fprintf('# of initial grains = %d\n',N);
for t = 1:max(list_t)% Main time iteration starts.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONVOLUTION STEP:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For fast convolution, far away grains are grouped into families.
    % Then, the entire family (union of grains) is convolved at once.
    [families,famgrains] = grains2families2doriginal(grains,Rfamilies,dims); % Group grains into families.
    N = size(families,1); % Number of families.
    
	start_progress(' - Convolving families')
    for k=1:N % Loop over families.
        cval = convolvefamily2doriginal(families{k,1},families{k,2},dims);
        % Distribute convolution values to the grains contained in this family:
        numberofgrains = size(famgrains{k},1); % Number of grains in this family.
        listofgrains = famgrains{k};           % Column vector of grain indices.
        for ell = 1:numberofgrains % Loop over grains contained in this family.
            label = listofgrains(ell);
            ind = grains{label,1};
            grains{label,3} = cval(ind); % Read off and record the convolution vals.
        end % (for ell) Loop over grains ends.
        display_progress(k,N,1);
    end % (for k) Loop over families ends.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REDISTRIBUTION STEP:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_progress(' - Determining neighborhood grains')
    % Find grains present in a nhd. of each grid pt.:
    presence = get_nhd_grains2doriginal(grains,dims(1)*dims(2));
    display_progress(1,1,1)
    
    start_progress(' - Distribution step')
    % Redistribution according to Esedoglu-Otto paper:
    updatelevelsetdata2dconstant(presence,grains,ID);
    display_progress(1,1,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REMOVE EMPTY GRAINS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Not necessary, but speeds up the algorithm.
    [grains,ID] = removeemptygrainsoriginal(grains,ID);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REFRESH GRAIN BUFFERS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level set data for each grain extends beyond the interface,
    % to a tubular neighborhood. That neighborhood has to be
    % updated once the interface has (potentially) moved.
    start_progress(' - Refreshing grain buffers')
    N = size(grains,1); % Number of grains.
    for k=1:N % Loop over grains.
        ind = grains{k,1}; % Pixels within a nhd. of grain.
        val = grains{k,2}; % Lev. set. vals. at those pixels.
        cval = grains{k,3}; % Convolution vals. at those pixels.
        Z(ind) = val;      % Lev. set. representation on grid.
        posind = ind(val>0); % Pixels in the interior of grain.
        [x,y] = ind2sub(dims,posind);
        [x2,y2] = pgrow3(int32(x),int32(y),Rgrains,WORKSPACE,...
            work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y); % Dilation.
        ind2 = sub2ind(dims,x2,y2);
        val2 = Z(ind2); % Level set vals.
        Z(ind2) = -1; % Reset Z back to all -1's.
        Z(ind) = cval - 1; % Convolution values - 1.
        cval2 = Z(ind2); % Convolution vals - 1.
        Z(ind2) = -1;
        grains{k,1} = ind2;   % Refresh grain's data structure.
        grains{k,2} = val2;   % Ditto.
        grains{k,3} = cval2 + 1; % Ditto.
        display_progress(k,N,1)
    end % (for k). Loop over grains ends.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save intermediate grains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = size(grains,1);
    fprintf('# of grains at iteration %d = %d\n',t,N);
    
    if ismember(t,list_t)
        epoch = epoch+1;
        N = size(grains,1); % Number of grains.
        start_progress([' - Saving data epoch ',num2str(epoch)])
        save_variable(filename,grains,['grains_t',num2str(epoch),'_constant'])
        save_variable(filename,t,['t',num2str(epoch),'_constant'])
        display_progress(1,1,1);
    end
end % (for t). Main time iteration ends.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN TIME LOOP ENDS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end