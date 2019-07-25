function S=hsvconereg(pts)
rho = atan2(pts(:,2),pts(:,1));
maxrho = 2*pi;
rho = mod(rho,maxrho)./maxrho;
c = hsv2rgb(rho, sqrt(pts(:,1).^2 + pts(:,2).^2)./(pts(:,3) + 1*(pts(:,3)==0)),...
pts(:,3));
S = reshape(c,[size(pts,1),3]);
S = 1-S;
end