function [grains,ori] = initialvoronoidata3d(N,dims,Rgrains,Rfamilies,th_ang,choice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [grains,ori] = initialvoronoidata3d(N,dims,Rgrains,Rfamilies,th_ang,choice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates initial condition with "N" grains on a grid of size "dims".
% Grains with misorientation angle smaller than th_ang degrees are merged.
% Each grain will have several disconnected components.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    choice = 1;
end

disp('Generating initial grain')

% Auxiliary vars.
global WORKSPACE Z LONGARRAY;

n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.

WORKSPACE = int32(-ones(dims)); % Needed for dilation.c.
Z = -ones(dims); % Another Workspace var.
% Needed in ls2vf3D.c.

LONGARRAY = int32(zeros(120000000,1));

label = -1e10 * ones(dims);

r = 0.1;
go = 1;
while go
    Naux = floor(N*(1+r));
    x = 1 + floor(rand(1,Naux)*(n1-1));
    y = 1 + floor(rand(1,Naux)*(n2-1));
    z = 1 + floor(rand(1,Naux)*(n3-1));
    
    % We ensure that there are not repeated seeds for the voronoi region
    aux = unique([x;y;z]','rows');
    x = aux(:,1)';
    y = aux(:,2)';
    z = aux(:,3)';
    Naux = length(x);
    
    if Naux>=N
        go = 0;
        select = randperm(Naux,N);
        x = x(select);
        y = y(select);
        z = z(select);
    else
        r = r + 0.1;
    end
end

grains = cell(N,1);

% Parameter d controls the extent to which the distance funciton
% to each point in the dataset is constructed:
d = 1 + ceil( (n1*n2*n3/N)^(1/3) );
d = 2*d;


start_progress(' - Constructing distance function');
% Construct distance function, label, to the union of all the points:
for k=1:N % Loop over the random points.

  [x2,y2,z2] = meshgrid((x(k)-d):(x(k)+d),(y(k)-d):(y(k)+d),(z(k)-d):(z(k)+d));
  dist = -sqrt( (x2-x(k)).^2 + (y2-y(k)).^2 + (z2-z(k)).^2 );
  x2 = 1 + mod(x2-1,n1);
  y2 = 1 + mod(y2-1,n2);
  z2 = 1 + mod(z2-1,n3);
  ind = sub2ind( dims , x2(:) , y2(:) , z2(:) );
  label(ind) = max(dist(:),label(ind));
  display_progress(k,N,1)
end % (for k). Loop over random points ends.

% If the union of d neighborhoods of the random points do not cover
% the entire computational domain, we cannot trust the construction:
if min(label(:)) < -0.5*1e10
  error('Parameter d too small.');
end

start_progress(' - Constructing Voronoi regions');
% Associate each grid point with the random point it is closest to,
% forming the grains:
indgridQ = zeros(dims);
for k=1:N % Loop over the random points again.

  [x2,y2,z2] = meshgrid((x(k)-d):(x(k)+d),(y(k)-d):(y(k)+d),(z(k)-d):(z(k)+d));
  dist = -sqrt( (x2-x(k)).^2 + (y2-y(k)).^2 + (z2-z(k)).^2 );
  x2 = 1 + mod(x2-1,n1);
  y2 = 1 + mod(y2-1,n2);
  z2 = 1 + mod(z2-1,n3);
  
  ind = sub2ind( dims , x2(:) , y2(:) , z2(:));
  ind2 = ind( dist(:) >= label(ind) );
  ind2 = ind2(not(indgridQ(ind2)));
  indgridQ(ind2) = 1;
  grains{k,1} = ind2;
  display_progress(k,N,1)

end % (for k). Loop over random points ends.

ori = cell(N,1);
parameters = zeros(N,3);
start_progress(' - Constructing orientations');
for k=1:N % Loop over the grains.
    parameters(k,:) = [rand rand rand];
    ori{k} = SO3(parameters(k,1),parameters(k,2),parameters(k,3));
    switch choice
        % Introduces texture by the method described in Gruber et. al.
        % with alpha = 4.2
        case 4
        while true
            minangle = misorientation3d_simple(ori{k});
            p = rand;
            if p < exp(-4.2*minangle)
                break
            else
                parameters(k,:) = [rand rand rand];
                ori{k} = SO3(parameters(k,1),parameters(k,2),parameters(k,3));
            end
        end
        % Introduces texture by the method described in Gruber et. al.
        % with alpha = 1
        case 5
            while true
                minangle = misorientation3d_simple(ori{k});
                p = rand;
                if p < exp(-minangle)
                    break
                else
                    parameters(k,:) = [rand rand rand];
                    ori{k} = SO3(parameters(k,1),parameters(k,2),parameters(k,3));
                end
            end
        % Attempt to perturb the initial MDF by imposing that 10% of the
        % grains have a certain misorientation angle.
        case 6
            if k/N < 0.1
                while true
                    minangle = misorientation3d_simple(ori{k});
                    if minangle/(2*pi)*360 < th_ang
                        break
                    else
                        parameters(k,:) = [rand rand rand];
                        ori{k} = SO3(parameters(k,1),parameters(k,2),parameters(k,3));
                    end
                end
            end
    end
    display_progress(k,N,1)
end % (for k). Loop over grains ends.

switch choice
    case 2
        density_target = @(x) mackenzie_function_perturbed(x);
        ori = gradient_descent_ori(dims,grains,ori,parameters,density_target,label,Rfamilies);
    case 3
        density_target = @(x) mackenzie_function_perturbed_2(x);
        ori = gradient_descent_ori(dims,grains,ori,parameters,density_target,label,Rfamilies);
end


minang = misorientation3d(ori,N);
minang(1:N+1:N^2) = Inf;

start_progress(' - Merging grains');
[angle_min,~] = min(minang(:));
angle_min = angle_min/(2*pi)*360;
while angle_min < th_ang
    select = find(minang(:)/(2*pi)*360 < th_ang);
    [iaux,jaux] = ind2sub(size(minang),select);
    i = min([iaux';jaux'])';
    j = max([iaux';jaux'])';
    
    [i,select,~] = unique(i);
    j = j(select);
    
    [j,select] = setdiff(j,i);
    i = i(select);
    
    ori(j) = [];
    minang(:,j) = [];
    minang(j,:) = [];
    grains_new = grains;
    for k = 1:length(i)
        grains_new{i(k),1} = [grains{i(k),1};grains{j(k),1}];
    end
    grains_new(j,:) = [];
    grains = grains_new;
    [angle_min,~] = min(minang(:));
    angle_min = angle_min/(2*pi)*360;
    display_progress(1000-floor(abs(angle_min-th_ang)/th_ang*1000),1000,1);
end
display_progress(100,100,1)

N = size(grains,1);
start_progress(' - Constructing grains');
% Dilate the grain neighborhood, and define a level set function
% on it (1 inside, -1 outside):
for k=1:N % Loop over the grains.    
    ind = grains{k,1};
    [x,y,z] = ind2sub(dims,ind);
    [x2,y2,z2] = dilation(int32(x),int32(y),int32(z),Rgrains,WORKSPACE,LONGARRAY); % Dilation.
    ind2 = sub2ind(dims,x2,y2,z2);
    label(ind2) = -1;
    label(ind) = 1;
    grains{k,1} = ind2;
    grains{k,2} = label(ind2);
    grains{k,3} = 0*ind2;       % Convolution vals. init to 0.
    grains{k,4} = 0*ind2;       % Convolution vals. init to 0.
    display_progress(k,N,1);
end % (for k). Loop over grains ends.


end