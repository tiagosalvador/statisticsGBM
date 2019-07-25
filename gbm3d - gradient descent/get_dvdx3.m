function r = get_dvdx3(x1,x2,x3)

r = [cos(2*pi*x2)/(2*sqrt(x3));sin(2*pi*x2)/(2*sqrt(x3));-1/(2*sqrt(1-x3))];

end