function r = get_dgidx1(x1,x2,x3)

v = get_v(x1,x2,x3);
dRdx1 = get_dRdx1(x1,x2,x3);
r = (2*(v*v')-eye(3))*dRdx1;

end