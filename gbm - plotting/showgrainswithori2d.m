function showgrainswithori2d(grains,dims,ori,id)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% showgrainswithori2d(grains,dims,ori,id)
% Last input variable, "id", is optional.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
  id = (1:1:size(grains,1))';
end

N = size(grains,1); % Number of grains.
u1 = zeros(dims);
u2 = zeros(dims);

for k=1:N % Loop over grains.
  ind = grains{k,1}; % Pixels in a nhd. of the grain.
  val = grains{k,2}; % Level set values. 
  ind2 = ind(val>0); % Pixels in the interior of grain.
  ang = ori(id(k));
  u1(ind2) = (1+cos(ang))/2;
  u2(ind2) = (1+sin(ang))/2;
end
U(:,:,1) = u1;
U(:,:,2) = u2;
U(:,:,3) = ones(dims);
clf
image(U);
axis square;
colormap jet

end