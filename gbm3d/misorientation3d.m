function minang = misorientation3d(ori,N,ID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minang = misorientation3d(ori,N,ID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the misorientation angle matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    ID = (1:N)';
end

c = cos(pi/2);
s = sin(pi/2);

r1 = [1 0 0;...
      0 c s;...
      0 -s c];

r2 = [c 0 s;...
      0 1 0;...
      -s 0 c];

r3 = [c s 0;...
      -s c 0;...
      0 0 1];

rot = zeros(3*24,3);

rot(1:3,:) = r1;
rot(4:6,:) = r1*r1;
rot(7:9,:) = r1*r1*r1;
rot(10:12,:) = r2;
rot(13:15,:) = eye(3,3);
rot(16:18,:) = r2*r2*r2;

rot(19:21,:) = r2*r1;
rot(22:24,:) = r3*r1*r1;
rot(25:27,:) = r2*r1*r1*r1;
rot(28:30,:) = r1*r2;
rot(31:33,:) = r3;
rot(34:36,:) = r1*r2*r2*r2;

rot(37:39,:) = r2*r2*r1;
rot(40:42,:) = r3*r3*r1*r1;
rot(43:45,:) = r2*r2*r1*r1*r1;
rot(46:48,:) = r1*r1*r2;
rot(49:51,:) = r3*r3;
rot(52:54,:) = r1*r1*r2*r2*r2;

rot(55:57,:) = r2*r2*r2*r1;
rot(58:60,:) = r3*r3*r3*r1*r1;
rot(61:63,:) = r2*r2*r2*r1*r1*r1;
rot(64:66,:) = r1*r1*r1*r2;
rot(67:69,:) = r3*r3*r3;
rot(70:72,:) = r1*r1*r1*r2*r2*r2;

start_progress(' - Computing misorientation angles')
minang = zeros(N,N);
dispstat(sprintf('minang initialized'),'timestamp');
for k = 1:(N-1)
    g1 = ori{ID(k)};
    g2 = vertcat(ori{ID((k+1):N)});
    minang(k,(k+1):N) = misorientation3d_aux(rot,g1*(g2'));
    minang((k+1):N,k) = minang(k,(k+1):N);
    display_progress(k,N,1);
end
display_progress(N,N,1);

end

function r = misorientation3d_aux(rot,q)

n = (length(q)/3);
r = min(real(acos((rot(1:3:72,:)*q(:,1:3:(3*n))+rot(2:3:72,:)*q(:,2:3:(3*n))+rot(3:3:72,:)*q(:,3:3:(3*n))-1)/2)),[],1);

end