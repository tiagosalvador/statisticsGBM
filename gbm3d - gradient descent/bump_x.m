function r = bump_x(x,delta)

%r = 1/delta^2*bump_aux_x(x/delta);

r = 1/delta*1/sqrt(2*pi)*(-x/delta^2).*exp(-(x/delta).^2/2);

end

% function r = bump_aux_x(x)
% 
% r = 1/sqrt(2*pi)*(-x).*exp(-x.^2/2);
% 
% end