function u = shownetwork2d(grains,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  u = shownetwork(grains,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(grains,1);
u = zeros(dims);

Z = -ones(dims);

for k=1:N
  ind = grains{k,1};
  val = grains{k,2};
  [x,y] = ind2sub(dims,ind);
  vf = loc_levset_to_volfluid(int32(x),int32(y),val,Z);
  u(ind) = u(ind) + vf.*(1-vf);
end

u = -u;

imagesc(u)
axis square
axis off
colormap gray(256)
