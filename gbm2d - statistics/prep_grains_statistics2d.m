function grains_statistics = prep_grains_statistics2d(grains,dims,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grains_statistics = prep_grains_statistics2d(grains,dims,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decomposes the grains into his components.

% Auxiliary vars.
global WORKSPACE KERNEL Z;

m = dims(1); % Size of computational grid.
n = dims(2); % Size of computational grid.

WORKSPACE = int32(-ones(dims)); % Needed for pgrow.c.
Z = -ones(m,n); % Another Workspace var.
                % Needed in loc_levset_to_volfluid.c.

global work_x work_y work_bx work_by work_candidate_x work_candidate_y; % Needed for pgrow3.c.
work_x = int32(zeros(prod(dims),1));
work_y = int32(zeros(prod(dims),1));
work_bx = int32(zeros(prod(dims),1));
work_by = int32(zeros(prod(dims),1));
work_candidate_x = int32(zeros(prod(dims),1));
work_candidate_y = int32(zeros(prod(dims),1));

grains_statistics = cell(0);
label = zeros(dims);
N = size(grains,1); % number of grain types
start_progress(' - Enlarging Grain''s')
for k=1:N
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind = ind(val>0); % Pixels in the interior of grain.
    [x,y] = ind2sub(dims,posind);
    [x2,y2] = pgrow3(int32(x),int32(y),Rfamilies,WORKSPACE,...
        work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y);
    ind2 = sub2ind(dims,x2,y2);
    label(ind2) = -1;
    label(posind) = 1;
    grains_statistics{k,1} = ind2;
    grains_statistics{k,2} = label(ind2);
    grains_statistics{k,3} = 0*ind2; % Convolution vals. init to 0.
    display_progress(k,N,1);
    display_progress(k,N,1);
end

% Prepare the convolution kernel KERNEL in fourier space:

dt = 2/n^2;
I = sqrt(-1); % Imaginary number i.
w=exp(I*2*pi/n); % nth root of unity.
[x,y] = meshgrid(1:n,1:n); x = x'; y = y';
KERNEL = exp(-dt*n*n*(4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y)));

% For fast convolution, far away grains are grouped into families.
% Then, the entire family (union of grains) is convolved at once.
[families,famgrains] = grains2families2doriginal(grains_statistics,Rfamilies,dims); % Group grains into families.
N = size(families,1); % Number of families.

start_progress(' - Convolving families')
for k=1:N % Loop over families.
    cval = convolvefamily2doriginal(families{k,1},families{k,2},dims);
    % Distribute convolution values to the grains contained in this family:
    numberofgrains = size(famgrains{k},1); % Number of grains in this family.
    listofgrains = famgrains{k};           % Column vector of grain indices.
    for ell = 1:numberofgrains % Loop over grains contained in this family.
        label = listofgrains(ell);
        ind = grains_statistics{label,1};
        grains_statistics{label,3} = cval(ind); % Read off and record the convolution vals.
    end % (for ell) Loop over grains ends.
    display_progress(k,N,1);
end % (for k) Loop over families ends.
end