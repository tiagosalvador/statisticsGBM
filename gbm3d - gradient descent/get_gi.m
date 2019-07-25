function gi = get_gi(x1,x2,x3)

v = get_v(x1,x2,x3);
R = get_R(x1,x2,x3);
gi = (2*(v*v')-eye(3))*R;

end