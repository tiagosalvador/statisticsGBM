%%%% Function: colormap432
%%%% Input: Array of Rodrigues Vector parameters
%%%% in 432 fundamental zone
%%%% Output: S = [r, g, b]
function S = colormap432(pts)
%%%% Corresponding to Equation C5(2)
kmax = sqrt(2)-1;
theta = atan2(pts(:,3),pts(:,2));
t2ind2 = find(pts(:,1) > 1/3 & theta > atan((1-2*pts(:,1))./(pts(:,1))));
tempvar2 = pts(t2ind2,2).*(1-pts(t2ind2,1));
tempvar3 = pts(t2ind2,1).*(pts(t2ind2,2)+pts(t2ind2,3))./(tempvar2 + 1*(tempvar2 == 0));
pts(t2ind2,:) = [pts(t2ind2,1) tempvar3.*pts(t2ind2,2) tempvar3.*pts(t2ind2,3)];
%%%% Corresponding to Equation C5(3)
g1 = vrrotvec2mat([1 0 0 -3*pi/8]);
pts = pts*g1;
pts = [pts(:,1) - kmax pts(:,2) pts(:,3)];
%%%% Corresponding to Equation C5(4)
tempvar = (1 + pts(:,2).*tan(pi/8)./(pts(:,3) + 1*(pts(:,3)==0)));
pts = [pts(:,1) pts(:,2).*tempvar pts(:,3).*tempvar];
%%%% Corresponding to Equation C5(5)
pts = [pts(:,1) pts(:,2)*cos(pi/8)/tan(pi/8) (pts(:,3) - pts(:,1)./cos(pi/8))];
%%%% Corresponding to Equation C5(6)
theta = atan2(-pts(:,1),pts(:,2));
pts = [pts(:,1).*(sin(theta) + abs(cos(theta))) ...
pts(:,2).*(sin(theta) + abs(cos(theta))) pts(:,3)];
%%%% Corresponding to Equation C5(7)
%keep pts
clearvars -except pts
theta = atan2(-pts(:,1),pts(:,2));rad = hypot(pts(:,2),pts(:,1));
pts = [-rad.*sin(2*theta) rad.*cos(2*theta) pts(:,3)];
%%%% Rescaling cone
%%%% Corresponding to Equation C5(8)
kmax = (sqrt(2)-1);tempk = cos(pi/8)/tan(pi/8);
pts(:,1) = pts(:,1)/kmax;pts(:,2) = pts(:,2)/kmax;pts(:,3) = pts(:,3)*tempk;
g1 = vrrotvec2mat([0 0 1 -pi/6]);
pts = pts*g1;
S = hsvconereg(pts);
%%%% Transforming cone and [r, g, b] coordinates to get intuitive colors
g1 = vrrotvec2mat([0 1 0 pi]);
S = S*g1;
S=[S(:,2) S(:,3)+1 S(:,1)+1];