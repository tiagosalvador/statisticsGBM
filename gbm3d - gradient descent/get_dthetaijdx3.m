function dthetaijdx3 = get_dthetaijdx3(i,j,ori,parameters,index,rot)


rot_aux = rot(((index(i,j)-1)*3+1):(index(i,j)*3),:);
cos_theta_ij = (trace(rot_aux*ori{i}*ori{j}')-1)/2;
x1 = parameters(i,1);
x2 = parameters(i,2);
x3 = parameters(i,3);
dgidx3 = get_dgidx3(x1,x2,x3);
dthetaijdx3 = -1/sqrt(1-(cos_theta_ij)^2)*(1/2)*trace(rot_aux*dgidx3*(ori{j}'));

end