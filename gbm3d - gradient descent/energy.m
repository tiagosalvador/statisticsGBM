function E = energy(x,f,g,l,theta,delta)
dx = x(2)-x(1);
E = sum(abs((f(x)-g(x,l,theta,delta))).^2)*dx;
end