function ori_new = gradient_descent_ori(dims,grains,ori,parameters,density_target,dist_funct,Rfamilies)

global WORKSPACE LONGARRAY

N = size(grains,1);
start_progress(' - Prep grain for gradient descent');
% Dilate the grain neighborhood, and define a level set function
% on it (1 inside, -1 outside):
for k=1:N % Loop over the grains.
    ind = grains{k,1};
    [x,y,z] = ind2sub(dims,ind);
    [x2,y2,z2] = dilation(int32(x),int32(y),int32(z),Rfamilies,WORKSPACE,LONGARRAY); % Dilation.
    ind2 = sub2ind(dims,x2,y2,z2);
    dist_funct(ind2) = -1;
    dist_funct(ind) = 1;
    grains{k,1} = ind2;
    grains{k,2} = dist_funct(ind2);
    grains{k,3} = 0*ind2;       % Convolution vals. init to 0.
    display_progress(k,N,1);
end % (for k). Loop over grains ends.

global KERNEL
% Prepare the convolution kernel KERNEL in fourier space:
n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.
dt = 2/n1^2;
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
KERNEL = exp( -dt*KERNEL );

% For fast convolution, far away grains are grouped into families.
% Then, the entire family (union of grains) is convolved at once.
[families,famgrains] = grains2families3doriginal(grains,Rfamilies,dims); % Group grains into families.
N = size(families,1); % Number of families.

start_progress(' - Convolving families')
for k=1:N % Loop over families.
    cval = convolvefamily3doriginal(families{k,1},families{k,2},dims);
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

grainlabels = zeros(dims); % Grid describing which grain each pixel belongs to.
N = size(grains,1);
start_progress(' - Prep to compute surface areas')
posind = cell(1,N);
for k = 1:N % Loop over grains.
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind{k} = ind(val>0); % Pixels in the interior of grain.
    grainlabels(posind{k}) = k; % Populate grid of grain labels.
    display_progress(k,N,1);
end % (for k). Loop over grains ends.


start_progress(' - Computing surface areas')
areas = zeros(N,N);
nNeighbors = zeros(N,1);
for k = 1:N
    indk = posind{k};
    [I,J,K] = ind2sub(dims,indk);
    [~,~,~,~,~,~,~,~,indup,...
        ~,~,~,~,~,~,~,~,inddown,...
        indN,~,indE,~,indS,~,indW,~] = find_ind_neighbors3d(I,J,K,dims);
    ind_neighbors = [indN;indE;indS;indW;indup;inddown]';
    neighbors = setdiff(unique(grainlabels(ind_neighbors)),k);
    cval = grains{k,3};
    nNeighbors(k) = length(neighbors);
    [ind,select] = setdiff(grains{k,1},indk);
    cval = cval(select);
    for l = neighbors
        [~,~,select] = intersect(posind{l},ind);
        % aux(k,l) = 1;
        areas(k,l) = sqrt(pi)*sum(cval(select))/prod(dims)/sqrt(dt);
    end % (for l).
    display_progress(k,N,1);
end % (for k).

start_progress(' - Vectorizing surface areas')
Nsurface_areas = N*(N+1)/2-N;
areas_vector = zeros(1,Nsurface_areas);
for i=1:N
    for j = (i+1):N
        areas_vector(j-1/2*i.*(1+i-2*N)-N) = areas(i,j);
        display_progress(j-1/2*i.*(1+i-2*N)-N,Nsurface_areas,1);
    end
end
max_iter = 50;
tolerance = 10^-4;
start_progress(' - Gradient descent')
ori_matrix = vertcat(ori{:});
[~,ori_matrix] = gradient_descent_aux(density_target,parameters,ori_matrix,areas_vector/sum(areas_vector),max_iter,tolerance);
ori_new = cell(N,1);
for k=1:N
    ori_new{k} = ori_matrix((3*k-2):3*k,:);
end


end

