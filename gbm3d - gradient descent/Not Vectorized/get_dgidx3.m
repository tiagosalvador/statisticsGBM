function r = get_dgidx3(x1,x2,x3)

v = get_v(x1,x2,x3);
dvdx3 = get_dvdx3(x1,x2,x3);
R = get_R(x1,x2,x3);
r = 2*(dvdx3*v'+v*dvdx3')*R;

end