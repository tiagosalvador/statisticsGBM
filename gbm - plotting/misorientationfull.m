function [minangle,axis] = misorientationfull(q)

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

rot = zeros(3,3,24);

rot(:,:,1) = r1;
rot(:,:,2) = r1*r1;
rot(:,:,3) = r1*r1*r1;
rot(:,:,4) = r2;
rot(:,:,5) = eye(3,3);
rot(:,:,6) = r2*r2*r2;

count = 6;
rot(:,:,count+1) = r2*r1;
rot(:,:,count+2) = r3*r1*r1;
rot(:,:,count+3) = r2*r1*r1*r1;
rot(:,:,count+4) = r1*r2;
rot(:,:,count+5) = r3;
rot(:,:,count+6) = r1*r2*r2*r2;

count = 12;
rot(:,:,count+1) = r2*r2*r1;
rot(:,:,count+2) = r3*r3*r1*r1;
rot(:,:,count+3) = r2*r2*r1*r1*r1;
rot(:,:,count+4) = r1*r1*r2;
rot(:,:,count+5) = r3*r3;
rot(:,:,count+6) = r1*r1*r2*r2*r2;

count = 18;
rot(:,:,count+1) = r2*r2*r2*r1;
rot(:,:,count+2) = r3*r3*r3*r1*r1;
rot(:,:,count+3) = r2*r2*r2*r1*r1*r1;
rot(:,:,count+4) = r1*r1*r1*r2;
rot(:,:,count+5) = r3*r3*r3;
rot(:,:,count+6) = r1*r1*r1*r2*r2*r2;

minangle = 1e100;

angles = zeros(1,24);
for i=1:24
    for j=1:24
        rot_aux = rot(:,:,i)*q*rot(:,:,j);
        theta = real(acos((trace(rot_aux)-1)/2));
        v = 1/(2*sin(theta))*[rot_aux(3,2)-rot_aux(2,3);rot_aux(1,3)-rot_aux(3,1);rot_aux(2,1)-rot_aux(1,2)];
        R = tan(theta/2)*v;
        if and(all(R>=0),and(R(1)<=(sqrt(2)-1),and(R(2)<=R(1),and(R(3)<=R(2),sum(R)<=1))))
            if theta < minangle
                minangle = theta;
                axis = v;
                bestrot1 = rot(:,:,i);
                bestrot2 = rot(:,:,j);
                besti = i;
                bestj = j;
            end
        end
    end
end
for i=1:24
    for j=1:24
        rot_aux = inv(rot(:,:,i)*q*rot(:,:,j));
        theta = real(acos((trace(rot_aux)-1)/2));
        v = 1/(2*sin(theta))*[rot_aux(3,2)-rot_aux(2,3);rot_aux(1,3)-rot_aux(3,1);rot_aux(2,1)-rot_aux(1,2)];
        R = tan(theta/2)*v;
        if and(all(R>=0),and(R(1)<=(sqrt(2)-1),and(R(2)<=R(1),and(R(3)<=R(2),sum(R)<=1))))
            if theta < minangle
                minangle = theta;
                axis = v;
                bestrot1 = rot(:,:,i);
                bestrot2 = rot(:,:,j);
                besti = i;
                bestj = j;
            end
        end
    end
end
theta = minangle;