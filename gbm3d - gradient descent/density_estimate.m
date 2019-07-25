function r = density_estimate(x,l,theta,delta)

% r = (bump(x'-theta,delta)*(l'))';

% needed for Matlab version on screms1
[a,b] = meshgrid(x',theta);
r = (((bump(a-b,delta))')*(l'))';

% r = 0;
% for i=1:length(l)
%     r = r+l(i)*bump(x-theta(i),delta);
% end

end