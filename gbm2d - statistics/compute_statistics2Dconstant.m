function [areas_grains,areas_convhull, nNeighbors, isopratio] = compute_statistics2Dconstant(grains,dims,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [areas_grains,areas_convhull, nNeighbors, isopratio] = compute_statistics2Dconstant(grains,dims,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes different grain statistics for 2D grains with model (iii)

N = size(grains,1);
aux = zeros(N,N);
grainlabels = zeros(dims); % Grid describing which grain each pixel belongs to.

start_progress(' - Computing areas')
areas_grains = zeros(N,1);
areas_convhull = zeros(N,1);
posind = cell(1,N);
for k = 1:N % Loop over grains.
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind{k} = ind(val>0); % Pixels in the interior of grain.
    grainlabels(posind{k}) = k; % Populate grid of grain labels.
    [areas_grains(k),areas_convhull(k)] = compute_areas(posind{k},dims);
    display_progress(k,N,1);
end % (for k). Loop over grains ends.

start_progress(' - Computing lengths')
lengths = zeros(N,N);
nNeighbors = zeros(N,1);
for k = 1:N
    indk = posind{k};
    [I,J] = ind2sub(dims,indk);
    [indN,~,indE,~,indS,~,indW,~] = find_ind_neighbors2d(I,J,dims);
    ind_neighbors = [indN;indE;indS;indW]';
    neighbors = setdiff(unique(grainlabels(ind_neighbors)),k);
    cval = grains{k,3};
    nNeighbors(k) = length(neighbors);
    [ind,select] = setdiff(grains{k,1},indk);
    cval = cval(select);
    for l = neighbors
        [~,~,select] = intersect(posind{l},ind);
        aux(k,l) = 1;
        lengths(k,l) = sqrt(pi)*sum(cval(select))/prod(dims)/sqrt(dt);
        % The lengths can alternatively be computed by simply counting the 
        % number of pixels with the following code:
        % length(intersect(ind_neighbors,posind{l}))*1/dims(1);
    end % (for l).
    display_progress(k,N,1,1);
end % (for k).

start_progress(' - Computing isopratio')
isopratio = zeros(N,1);
for k = 1:N
    d = length(dims);
    A = areas_grains(k);
    L = (sum(lengths(k,:))+sum(lengths(:,k)))/2;
    isopratio(k) = 4*pi*A^(d-1)/L^d;
    display_progress(k,N,1,1);
end % (for k).

end
