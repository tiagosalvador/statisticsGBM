function sphcoord = polarcoord(v)
%%%%% v = n x 3 array of pts
%%%%% output is n x 3 array of spherical coordinates
%%%%% [r theta phi]
r = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
ind1 = find(r == 0);
ind2 = find(r ~= 0);
if(size(ind1,2) > 0)
    theta(ind1)=0; phi(ind1)=0;
end
if(size(ind2,2) > 0)
    theta(ind2) = acos(v(ind2,3)./r(ind2));
    phi(ind2) = atan2(v(ind2,2),v(ind2,1));
end
sphcoord = [r theta' phi'];
end