function S = surface_tension3d(ori,angBrandon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = surface_tension3d(ori,angBrandon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the surface tension matrix.

N = length(ori);

minang = misorientation3d(ori,N);

start_progress(' - Computing surface tension')
S = ones(N,N);
dispstat(sprintf('Surface tension matrix initialized'),'timestamp');

% Read-Shockley with Bradon angle %
angBrandonRad = angBrandon*pi/180; % 30 degrees 
select = minang<=angBrandonRad;
S(select) = minang(select)/angBrandonRad.*(1-log(minang(select)/angBrandonRad));
S(1:N+1:N^2) = 0;
display_progress(N,N,1);

end