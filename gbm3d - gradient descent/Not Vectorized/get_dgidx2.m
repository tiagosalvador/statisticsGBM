function r = get_dgidx2(x1,x2,x3)

v = get_v(x1,x2,x3);
dvdx2 = get_dvdx2(x1,x2,x3);
R = get_R(x1,x2,x3);
r = 2*(dvdx2*v'+v*dvdx2')*R;

end