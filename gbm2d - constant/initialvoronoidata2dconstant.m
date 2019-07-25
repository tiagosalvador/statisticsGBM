function grains = initialvoronoidata2dconstant(N,dims,Rgrains)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grains = initialvoronoidata2dconstant(N,dims,Rgrains)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates initial condition with "N" grains on a grid of size "dims"
% to be used in the function 'gbm2dconstant_sim'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


label = -1e10 * ones(dims);
m = dims(1); n = dims(2);
W = int32(-ones(dims));
work_x = int32(zeros(prod(dims),1));
work_y = int32(zeros(prod(dims),1));
work_bx = int32(zeros(prod(dims),1));
work_by = int32(zeros(prod(dims),1));
work_candidate_x = int32(zeros(prod(dims),1));
work_candidate_y = int32(zeros(prod(dims),1));

r = 0.1;
go = 1;
while go
    Naux = floor(N*(1+r));
    x = 1 + floor(rand(1,Naux)*(m-1));
    y = 1 + floor(rand(1,Naux)*(n-1));
    
    % We ensure that there are not repeated seeds for the voronoi region
    aux = unique([x;y]','rows','stable');
    x = aux(:,1)';
    y = aux(:,2)';
    Naux = length(x);
    
    if Naux>=N
        go = 0;
        select = randperm(Naux,N);
        x = x(select);
        y = y(select);
    else
        r = r + 0.1;
    end
end

grains = cell(N,1);

% Parameter d controls the extent to which the distance funciton
% to each point in the dataset is constructed:
d = 1 + ceil( sqrt(m*n/N) );
d = 2*d;

start_progress(' - Constructing distance function');
% Construct distance function, label, to the union of all the points:
for k=1:N % Loop over the random points.
    [x2,y2] = meshgrid((x(k)-d):(x(k)+d),(y(k)-d):(y(k)+d));
    dist = -sqrt( (x2-x(k)).^2 + (y2-y(k)).^2 );
    x2 = 1 + mod(x2-1,m);
    y2 = 1 + mod(y2-1,n);
    ind = sub2ind(dims,x2(:),y2(:));
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
    
    [x2,y2] = meshgrid((x(k)-d):(x(k)+d),(y(k)-d):(y(k)+d));
    dist = -sqrt( (x2-x(k)).^2 + (y2-y(k)).^2 );
    x2 = 1 + mod(x2-1,m);
    y2 = 1 + mod(y2-1,n);
    
    ind = sub2ind(dims,x2(:),y2(:));
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
    [x,y] = ind2sub(dims,ind);
    [x2,y2] = pgrow3(int32(x),int32(y),Rgrains,W,...
        work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y);
    ind2 = sub2ind(dims,x2,y2);
    label(ind2) = -1;
    label(ind) = 1;
    grains{k,1} = ind2;
    grains{k,2} = label(ind2);
    grains{k,3} = 0*ind2;       % Convolution vals. of alpha kernel init to 0.
    display_progress(k,N,1);
end % (for k). Loop over grains ends.