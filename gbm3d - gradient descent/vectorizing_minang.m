function theta = vectorizing_minang(minang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta = vectorizing_minang(minang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transforms the symmetric misorienation angles matrix into a vector taking
% advantage of its symmetric form.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(minang,1);

Nminang = N*(N+1)/2-N;
theta = zeros(1,Nminang);
for i=1:N
    for j = (i+1):N
        theta(j-1/2*i.*(1+i-2*N)-N) = minang(i,j);
    end
end

end
