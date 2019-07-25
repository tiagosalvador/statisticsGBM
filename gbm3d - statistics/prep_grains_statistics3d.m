function grains_statistics = prep_grains_statistics3d(grains,dims,Rfamilies,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grains_statistics = prep_grains_statistics3d(grains,dims,Rfamilies,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decomposes the grains into his components.

% Auxiliary vars.
global WORKSPACE KERNEL Z LONGARRAY;

n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.

WORKSPACE = int32(-ones(dims)); % Needed for dilation.c.
LONGARRAY = int32(zeros(120000000,1)); % Needed for dilation.c.
Z = -ones(dims); % Needed in ls2vf3D.c.

grains_statistics = cell(0);
label = zeros(dims);
N = size(grains,1); % number of grain types
start_progress(' - Enlarging Grain''s')
for k=1:N
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind = ind(val>0); % Pixels in the interior of grain.
    
    [x,y,z] = ind2sub(dims,posind);
    [x2,y2,z2] = dilation(int32(x),int32(y),int32(z),Rfamilies,WORKSPACE,LONGARRAY); % Dilation.
    ind2 = sub2ind(dims,x2,y2,z2);
    label(ind2) = -1;
    label(posind) = 1;
    grains_statistics{k,1} = ind2;
    grains_statistics{k,2} = label(ind2);
    grains_statistics{k,3} = 0*ind2; % Convolution vals. init to 0.
    display_progress(k,N,1);
end

% Prepare the convolution kernel KERNEL in fourier space:

I = sqrt(-1);      % Imaginary number i.
wx=exp(I*2*pi/n1); % nth root of unity.
wy=exp(I*2*pi/n2); % nth root of unity.
wz=exp(I*2*pi/n3); % nth root of unity.
[x,y,z]=meshgrid(1:n1,1:n2,1:n3);
x=permute(x,[2 1 3]);
y=permute(y,[2 1 3]);
z=permute(z,[2 1 3]);
aux = n1*n1*(2-wx.^(x-1)-wx.^(1-x)) + n2*n2*(2-wy.^(y-1)-wy.^(1-y)) + ...
    n3*n3*(2-wz.^(z-1)-wz.^(1-z));
KERNEL =  exp(-dt*aux);
% For fast convolution, far away grains are grouped into families.
% Then, the entire family (union of grains) is convolved at once.
[families,famgrains] = grains2families3doriginal(grains_statistics,Rfamilies,dims); % Group grains into families.
N = size(families,1); % Number of families.

start_progress(' - Computing convolution values')
for k=1:N % Loop over families.
    cval = convolvefamily3doriginal(families{k,1},families{k,2},dims);
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