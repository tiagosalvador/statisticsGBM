function [areas_grains,areas_convhull, nNeighbors,...
	isopratio] = compute_statistics_slices3D_auxconstant(grains_section,ID_section,dims_section,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [areas_grains,areas_convhull, nNeighbors,...
%	isopratio] = compute_statistics_slices3D_auxconstant(grains_section,ID_section,dims_section,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(grains_section,1);
aux = zeros(N,N);
grainlabels = zeros(dims_section); % Grid describing which grain each pixel belongs to.

areas_grains = zeros(N,1);
areas_convhull = zeros(N,1);
posind = cell(1,N);
for k = 1:N % Loop over grains.
    ind = grains_section{k,1}; % Pixels within a nhd. of grain.
    val = grains_section{k,2}; % Lev. set. vals. at those pixels.
    posind{k} = ind(val>0); % Pixels in the interior of grain.
    grainlabels(posind{k}) = k; % Populate grid of grain labels.
    [areas_grains(k),areas_convhull(k)] = compute_areas(posind{k},dims_section);
end % (for k). Loop over grains ends.

lengths = zeros(N,N);
nNeighbors = zeros(N,1);
for k = 1:N
    indk = posind{k};
    [I,J] = ind2sub(dims_section,indk);
    [indN,~,indE,~,indS,~,indW,~] = find_ind_neighbors2d(I,J,dims_section);
    ind_neighbors = [indN;indE;indS;indW]';
    neighbors = setdiff(unique(grainlabels(ind_neighbors)),k);
    cval = grains_section{k,3};
    nNeighbors(k) = length(neighbors);
    [ind,select] = setdiff(grains_section{k,1},indk);
    cval = cval(select);
    for l = neighbors
        [~,~,select] = intersect(posind{l},ind);
        aux(k,l) = 1;
        lengths(k,l) = sqrt(pi)*sum(cval(select))/prod(dims_section)/sqrt(dt);
        %lengths_pixels(k,l) = length(intersect(ind_neighbors,posind{l}))*1/(dims_section(1));
    end % (for l).
end % (for k).


isopratio = zeros(N,1);
for k = 1:N
    d = length(dims_section);
    A = areas_grains(k);
    L = (sum(lengths(k,:))+sum(lengths(:,k)))/2;
    isopratio(k) = 4*pi*A^(d-1)/L^d;
end % (for k).