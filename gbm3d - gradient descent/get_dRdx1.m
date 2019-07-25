function r = get_dRdx1(x1,x2,x3)

r = zeros(3,3);
r(1,1) = -2*pi*sin(2*pi*x1);
r(1,2) = 2*pi*cos(2*pi*x1);
r(2,1) = -2*pi*cos(2*pi*x1);
r(2,2) = -2*pi*sin(2*pi*x1);

end