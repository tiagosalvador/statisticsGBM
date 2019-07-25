function M = SO3(x1,x2,x3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M = SO3(x1,x2,x3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine maps three values (x1,x2,x3) in the range [0,1]
% into a 3x3 rotation matrix, M. Uniformly distributed random variables
% x0, x1, and x2 create uniformly distributed random rotation matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = zeros(3,3);
ct = cos(2*pi*x1);
st = sin(2*pi*x1);
R(1,1) = ct;
R(1,2) = st;
R(2,1) = -st;
R(2,2) = ct;
R(3,3) = 1;

v = zeros(3,1);
v(1) = cos(2*pi*x2)*sqrt(x3);
v(2) = sin(2*pi*x2)*sqrt(x3);
v(3) = sqrt(1-x3);

M = (2*(v*v')-eye(3))*R;

end





















%% 
% Unclear if R will be uniformly distributed when v is chosen uniformly on
% the 2-sphere and theta is uniformly chosen in [0,2*pi]

% function R = SO3(v,theta)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % R = SO3(v,theta)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % S03 matrix for a rotation by an angle of theta about an axis in the
% % direction of v
% %
% % INPUT:
% %   v = axis of the rotation
% %   theta = angle of the rotation
% 
% vx = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
% R = cos(theta)*eye(3,3) + sin(theta)*vx + (1-cos(theta))*v'*v;
% 
% end