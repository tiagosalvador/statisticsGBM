function minangle = misorientation3d_simple(q)

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
  
  angle = real( acos( (trace(rot(:,:,i)*q)-1 )/2 ) );
  angles(i) = angle;
  if angle < minangle
    minangle = angle;
  end
  
end

end