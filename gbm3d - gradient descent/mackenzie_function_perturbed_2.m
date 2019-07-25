function r = mackenzie_function_perturbed_2(w)

wmax = 2*atan(sqrt(23-16*sqrt(2)));
w1 = 0.3;
wM = 45/360*2*pi; % maximum of Mackenzie distribution

wnew = (w<=w1).*w*wM/w1 + (w>w1).*((w-wmax)/(abs(w1-wmax))*(abs(wM-wmax))+wmax);

r = mackenzie_function(wnew);

end




