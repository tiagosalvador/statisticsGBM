function [dgidx1,dgidx2,dgidx3] = get_dgi(x1,x2,x3)

v = get_v(x1,x2,x3);
R = get_R(x1,x2,x3);

dRdx1 = get_dRdx1(x1,x2,x3);
dgidx1 = (2*(v*v')-eye(3))*dRdx1;

dvdx2 = get_dvdx2(x1,x2,x3);
dgidx2 = 2*(dvdx2*v'+v*dvdx2')*R;

dvdx3 = get_dvdx3(x1,x2,x3);
dgidx3 = 2*(dvdx3*v'+v*dvdx3')*R;

end