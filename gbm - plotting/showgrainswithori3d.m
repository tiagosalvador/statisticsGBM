function showgrainswithori3d(grains,dims,ori,index_reference,id)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% showgrainswithori3d(grains,dims,ori,index_reference,id)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
  id = (1:1:size(grains,1))';
end

N = size(grains,1);

xgrid = linspace(0,1,dims(1));
ygrid = xgrid;
zgrid = xgrid;
for k=1:N
    disp(k)
    ind = grains{k,1};
    val = grains{k,2};
    ind2 = ind(val>0);
    [x,y,z] = ind2sub(dims,ind2);
    
    if id(k) == index_reference
        theta = 0;
        v = [0;0;0];
    else
        g1 = ori{id(k)};
        g2 = ori{index_reference};
        [theta,v] = misorientationfull(g1*(g2'));
    end
    c = max(min(colormap432(tan(theta/2)*v'),1),0);
    scatter3(xgrid(x),ygrid(y),zgrid(z),2,c,'fill')
    axis([0 1 0 1 0 1]);
    axis square;
    hold on
end
aux_zeros = zeros(size(xgrid));
aux_ones = ones(size(xgrid));
s = 0.2;
scatter3(xgrid,aux_zeros,aux_zeros,s,'black')
scatter3(xgrid,aux_zeros,aux_ones,s,'black')
scatter3(xgrid,aux_ones,aux_zeros,s,'black')
scatter3(xgrid,aux_ones,aux_ones,s,'black')
scatter3(aux_zeros,ygrid,aux_zeros,s,'black')
scatter3(aux_zeros,ygrid,aux_ones,s,'black')
scatter3(aux_ones,ygrid,aux_zeros,s,'black')
scatter3(aux_ones,ygrid,aux_ones,s,'black')
scatter3(aux_zeros,aux_zeros,zgrid,s,'black')
scatter3(aux_zeros,aux_ones,zgrid,s,'black')
scatter3(aux_ones,aux_zeros,zgrid,s,'black')
scatter3(aux_ones,aux_ones,zgrid,s,'black')

% 
% N = size(grains,1); % Number of grains.
% u1 = zeros(dims);
% u2 = zeros(dims);
% 
% for k=1:N % Loop over grains.
%   ind = grains{k,1}; % Pixels in a nhd. of the grain.
%   val = grains{k,2}; % Level set values. 
%   ind2 = ind(val>0); % Pixels in the interior of grain.
%   ang = ori(id(k));
%   u1(ind2) = (1+cos(ang))/2;
%   u2(ind2) = (1+sin(ang))/2;
% end
% U(:,:,1) = u1;
% U(:,:,2) = u2;
% U(:,:,3) = ones(dims);
% clf
% image(U);
% axis square;
% colormap jet

end