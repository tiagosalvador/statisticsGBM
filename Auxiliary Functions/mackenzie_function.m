function r = mackenzie_function(w)

r = 24/(2*pi^2)*sin(w/2).^2.*chi(tan(w/2));

end

function r = chi(tau)

r = zeros(size(tau));
tau0 = 0;
tau1 = sqrt(2)-1;
tau2 = sqrt(3)/3;
tau3 = 2-sqrt(2);
tau4 = sqrt(23-16*sqrt(2));

select = and(tau0<=tau,tau<tau1);
r(select) = chi1(tau(select));
select = and(tau1<=tau,tau<tau2);
r(select) = chi1(tau(select))+chi2(tau(select));
select = and(tau2<=tau,tau<tau3);
r(select) = chi1(tau(select))+chi2(tau(select))+chi3(tau(select));
select = and(tau3<=tau,tau<tau4);
r(select) = chi1(tau(select))+chi2(tau(select))+chi3(tau(select))+chi4(tau(select));

end


function r = chi1(tau)

r = 4*pi;

end

function r = chi2(tau)

r = -6*S1(alpha1(tau));

end

function r = chi3(tau)

r = -8*S1(alpha2(tau));

end

function r = chi4(tau)

delta1 = pi/2;
delta2 = acos(sqrt(3)/3);
r = 12*S2(alpha1(tau),alpha1(tau),delta1)+24*S2(alpha1(tau),alpha2(tau),delta2);

end

function r = alpha1(tau)

tau1 = sqrt(2)-1;
r = acos(tau1./tau);

end

function r = alpha2(tau)

tau2 = sqrt(3)/3;
r = acos(tau2./tau);

end


function r = S1(alpha)

r = 2*pi*(1-cos(alpha));

end

function r = S2(alpha,beta,gamma)

r = 2*(pi-C(alpha,beta,gamma)-cos(alpha).*C(gamma,alpha,beta)-cos(beta).*C(beta,gamma,alpha));

end

function r = C(alpha,beta,gamma)

r = acos((cos(gamma)-cos(alpha).*cos(beta))./(sin(alpha).*sin(beta)));

end




