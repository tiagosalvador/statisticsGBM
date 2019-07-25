function [volumes_grains, volumes_convhull,...
    volumes_grains_alt, volumes_convhull_alt,...
    surface_areas_grains, nNeighbors] = ...
    compute_statistics3Dconstant(grains,dims,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [volumes_grains, volumes_convhull,...
%    volumes_grains_alt, volumes_convhull_alt,...
%    surface_areas_grains, nNeighbors] = ...
%    compute_statistics3Dconstant(grains,dims,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes different grain statistics for 3D grains with model (iii)

N = size(grains,1);
aux = zeros(N,N);
grainlabels = zeros(dims); % Grid describing which grain each pixel belongs to.

start_progress(' - Computing volumes')
volumes_grains = zeros(N,1);
volumes_convhull = zeros(N,1);
volumes_grains_alt = zeros(N,1);
volumes_convhull_alt = zeros(N,1);
posind = cell(1,N);
for k = 1:N % Loop over grains.
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind{k} = ind(val>0); % Pixels in the interior of grain.
    grainlabels(posind{k}) = k; % Populate grid of grain labels.
    [volumes_grains(k),volumes_convhull(k),volumes_grains_alt(k),volumes_convhull_alt(k)] = compute_volumes(posind{k},dims);
    display_progress(k,N,1);
end % (for k). Loop over grains ends.

start_progress(' - Computing surface areas')
areas = zeros(N,N);
surface_areas_grains = zeros(N,1);
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
        aux(k,l) = 1;
        areas(k,l) = sqrt(pi)*sum(cval(select))/prod(dims)/sqrt(dt);
        %areas_pixels(k,l) = length(intersect(ind_neighbors,posind{l}))*1/(dims(1)*dims(2));
    end % (for l).
    display_progress(k,N,1);
end % (for k).

for k = 1:N
    surface_areas_grains(k) = (sum(areas(k,:))+sum(areas(:,k)))/2;
end

end
