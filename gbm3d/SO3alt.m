function R = SO3alt(v,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = SO3(v,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S03 matrix for a rotation by an angle of theta about an axis in the
% direction of v
%
% INPUT:
%   v = axis of the rotation
%   theta = angle of the rotation

vx = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
R = cos(theta)*eye(3,3) + sin(theta)*vx + (1-cos(theta))*v'*v;

end