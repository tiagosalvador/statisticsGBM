function gbm3dconstant_sim(dt,grains,dims,Rgrains,Rfamilies,list_t,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gbm3dconstant_sim(dt,grains,dims,Rgrains,Rfamilies,list_t,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D grain boundary motion with many grains. Isotropic, equal
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
%   max_epoch = maximum number of epochs. Everytime the total number of
%               grains decreases by a factor of 4 is consider an epoch.
%   Rgrains = radius by each a grain is enlarged.
%   Rfamilies = Minimum distance to maintain between grains in the same 
%               family.
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
global WORKSPACE KERNEL Z LONGARRAY;

n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.
n1n2n3 = n1*n2*n3;

WORKSPACE = int32(-ones(dims)); % Needed for dilation.c.
LONGARRAY = int32(zeros(120000000,1)); % Needed for dilation.c.
Z = -ones(dims); % Another Workspace var.
% Needed in ls2vf3D.c.

ID = [1:1:size(grains,1)]'; % Used for keeping track of grains
% as some of them disappear.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the convolution kernel KERNEL in fourier space:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = sqrt(-1);      % Imaginary number i.
wx=exp(I*2*pi/n1); % nth root of unity.
wy=exp(I*2*pi/n2); % nth root of unity.
wz=exp(I*2*pi/n3); % nth root of unity.
[x,y,z]=meshgrid(1:n1,1:n2,1:n3);
x=permute(x,[2 1 3]);
y=permute(y,[2 1 3]);
z=permute(z,[2 1 3]);
KERNEL = n1*n1*(2-wx.^(x-1)-wx.^(1-x)) + n2*n2*(2-wy.^(y-1)-wy.^(1-y)) + ...
    n3*n3*(2-wz.^(z-1)-wz.^(1-z));
KERNEL = exp(-dt*KERNEL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN TIME LOOP STARTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N0 = size(grains,1);
epoch = 0;
fprintf('# of initial grains = %d\n',N0);
for t = 1:max(list_t)% Main time iteration starts.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONVOLUTION STEP:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For fast convolution, far away grains are grouped into families.
    % Then, the entire family (union of grains) is convolved at once.
    % Finally, the families are broken up into individual grains again.
    [families,famgrains] = grains2families3doriginal(grains,Rfamilies,dims); % Group grains into families.
    N = size(families,1); % Number of families.
    
    start_progress(' - Convolving families')
    for k=1:N % Loop over families.
        cval = convolvefamily3doriginal(families{k,1},families{k,2},dims); % Calulate convolution.
        % Distribute convolution values to the grains contained in this family:
        numberofgrains = size(famgrains{k},1); % Number of grains in this family.
        listofgrains = famgrains{k}; % Column vector of grain indices.
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
    presence = get_nhd_grains3doriginal(grains,dims(1)*dims(2)*dims(3));
    display_progress(1,1,1)
    
    start_progress(' - Distribution step')
    % Redistribution according to Esedoglu-Otto paper:
    updatelevelsetdata3dconstant(presence,grains,ID);
    display_progress(1,1,1)

    % Clear memory not needed
    clear presence
    
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
        [x,y,z] = ind2sub(dims,posind);
        [x2,y2,z2] = dilation(int32(x),int32(y),int32(z),Rgrains,WORKSPACE,LONGARRAY); % Dilation.
        ind2 = sub2ind(dims,x2,y2,z2);
        val2 = Z(ind2);    % Level set vals.
        Z(ind2) = -1;      % Reset Z back to all -1's.
        Z(ind) = cval - 1; % Convolution values - 1.
        cval2 = Z(ind2);   % Convolution vals - 1.
        Z(ind2) = -1;      % Reset Z back to all -1's.
        grains{k,1} = ind2;   % Refresh grain's data structure.
        grains{k,2} = val2;   % Ditto.
        grains{k,3} = cval2 + 1; % Ditto.
        display_progress(k,N,1)
    end % (for k). Loop over grains ends.

    N = size(grains,1);
    fprintf('# of grains at iteration %d = %d\n',t,N);
    if ismember(t,list_t)
        epoch = epoch+1;
        start_progress([' - Saving data epoch ',num2str(epoch)])
        save_variable(filename,grains,['grains_t',num2str(epoch),'_constant'])
        save_variable(filename,ID,['ID_t',num2str(epoch),'_constant'])
        save_variable(filename,t,['t',num2str(epoch),'_constant'])
        display_progress(1,1,1);
    end
end % (for t). Main time iteration ends.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN TIME LOOP ENDS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end