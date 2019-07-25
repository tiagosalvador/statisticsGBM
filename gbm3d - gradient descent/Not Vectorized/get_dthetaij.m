function [dthetaijdx1,dthetaijdx2,dthetaijdx3] = get_dthetaij(i,j,theta_ij,ori,dgidx1,dgidx2,dgidx3,index,rot)

rot_aux = rot(((index(i,j)-1)*3+1):(index(i,j)*3),:);
iaux = min(i,j);
jaux = max(i,j);
N = length(ori);
cos_theta_ij = cos(theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N));

aux = -1/sqrt(1-(cos_theta_ij)^2)*(1/2);
q = (ori{j}')*rot_aux;

aux_matrix = [dgidx1;dgidx2;dgidx3];
aux = aux*(aux_matrix(1:3:9,:)*q(:,1)+aux_matrix(2:3:9,:)*q(:,2)+aux_matrix(3:3:9,:)*q(:,3));
dthetaijdx1 = aux(1);
dthetaijdx2 = aux(2);
dthetaijdx3 = aux(3);

end