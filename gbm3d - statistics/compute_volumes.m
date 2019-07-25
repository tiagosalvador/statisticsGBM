function [volume_grain,volume_convhull,volume_grain_alt,volume_convhull_alt] = compute_volumes(ind,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [volume_grain,volume_convhull,volume_grain_alt,volume_convhull_alt] = compute_volumes(ind,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes both the volume of a grain and the volume of its convex envelope.

[I,J,K] = ind2sub(dims,ind);

if sum(I==dims(1))
    x = 1;
    select = (I==x);
    while and(sum(select),x<=dims(1))
        I(select) = I(select)+dims(1);
        x = x+1;
        select = (I==x);
    end
end

if sum(J==dims(2))
    y = 1;
    select = (J==y);
    while and(sum(select),y<=dims(2))
        J(select) = J(select)+dims(2);
        y = y+1;
        select = (J==y);
    end
end

if sum(K==dims(3))
    z = 1;
    select = (K==z);
    while and(sum(select),z<=dims(3))
        K(select) = K(select)+dims(3);
        z = z+1;
        select = (K==z);
    end
end


% The volumes are computed in the same spirit as areas are
% computed for 2D grains. An alternative way using MATLAB
% native function 'boundary' is also used.
[~,volume_grain_alt] = boundary(I,J,K);
volume_grain_alt = volume_grain_alt/prod(dims);
volume_grain = length(ind)/prod(dims);

try
    grainpoints = [I J K];
    [conv_grain,volume_convhull_alt] = convhull(I,J,K);
    volume_convhull_alt = volume_convhull_alt/prod(dims);
    try
        [x,y,z] = meshgrid(min(I):max(I),min(J):max(J),min(K):max(K));
        tol = 1.e-13*mean(abs(grainpoints(:)));
        in = inhull([x(:) y(:) z(:)],[I J K],conv_grain,tol);
        volume_convhull = sum(in>0)/prod(dims);
    catch
        volume_convhull = volume_grain;
    end

catch
    volume_convhull_alt = volume_grain_alt;
    volume_convhull = volume_grain;
end

end


