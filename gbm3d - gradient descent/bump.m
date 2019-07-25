function r = bump(x,delta)

%r = 1/delta*bump_aux(x/delta);

r = 1/delta*1/sqrt(2*pi)*exp(-(x/delta).^2/2);


end

% function r = bump_aux(x)
% 
% r = 1/sqrt(2*pi)*exp(-x.^2/2);
% 
% end