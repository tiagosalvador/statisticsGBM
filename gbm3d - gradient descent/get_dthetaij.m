function [dthetaijdx1,dthetaijdx2,dthetaijdx3] = get_dthetaij(i,j,theta_ij,ori,dgidx1,dgidx2,dgidx3,index,rot)

n = length(j);
rot_aux = zeros(3*n,3);
rot_aux(1:3:(3*n),:) = rot(((index(i,j)-1)*3+1),:);
rot_aux(2:3:(3*n),:) = rot(((index(i,j)-1)*3+2),:);
rot_aux(3:3:(3*n),:) = rot(((index(i,j)-1)*3+3),:);

iaux = min(i,j);
jaux = max(i,j);
N = size(ori,1)/3;
cos_theta_ij = cos(theta_ij(jaux-1/2*iaux.*(1+iaux-2*N)-N));

aux = -1./sqrt(1-(cos_theta_ij).^2)*(1/2);

q = zeros(3*n,3);
q(1:3:(3*n),:) = ori(3*j-2,:);
q(2:3:(3*n),:) = ori(3*j-1,:);
q(3:3:(3*n),:) = ori(3*j,:);
q = q';
A = (rot_aux').*(dgidx1*q);
dthetaijdx1 = aux.*sum(reshape(sum(A),3,n));
A = (rot_aux').*(dgidx2*q);
dthetaijdx2 = aux.*sum(reshape(sum(A),3,n));
A = (rot_aux').*(dgidx3*q);
dthetaijdx3 = aux.*sum(reshape(sum(A),3,n));

end