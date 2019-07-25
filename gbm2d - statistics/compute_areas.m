function [area_grain,area_convhull] = compute_areas(ind,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [area_grain,area_convhull] = compute_areas(ind,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes both the area of a grain and the area of its convex envelope.

[I,J] = ind2sub(dims,ind);

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

area_grain = length(ind)/prod(dims);

try
    [k,~] = convhull(I,J);
    [x,y] = meshgrid(min(I):max(I),min(J):max(J));
    in = inpolygon(x(:),y(:),I(k),J(k));
    area_convhull = sum(in>0)/prod(dims);
catch
    area_convhull = area_grain;
end
end


