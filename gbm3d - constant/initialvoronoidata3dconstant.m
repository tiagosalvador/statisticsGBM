function grains = initialvoronoidata3dconstant(N,dims,Rgrains)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grains = initialvoronoidata3d(N,dims,Rgrains,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates initial condition with "N" grains on a grid of size "dims".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    display_progress(k,N,1);
end % (for k). Loop over grains ends.

end