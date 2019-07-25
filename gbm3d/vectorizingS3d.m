function Svector = vectorizingS3d(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Svector = vectorizingS3d(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transforms the symmetric surface tension matrix into a vector taking
% advantage of its symmetric form.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(S,1);

start_progress(' - Vectorizing S')
NSvector = N*(N+1)/2-N;
Svector = zeros(NSvector,1);
for i=1:N
    for j = (i+1):N
        Svector(j-1/2*i.*(1+i-2*N)-N,1) = S(i,j);
        display_progress(j-1/2*i.*(1+i-2*N)-N,NSvector,1);
    end
end

end
