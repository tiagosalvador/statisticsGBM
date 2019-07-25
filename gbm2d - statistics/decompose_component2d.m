function grains = decompose_component2d(dims,ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grains = decompose_component2d(dims,ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the components in the set with linear indices ind.

[I,J] = ind2sub(dims,ind);

still_left = sparse(I,J,1,dims(1),dims(2),length(ind));

grain = sparse(I,J,1,dims(1),dims(2),length(ind));

nleft = length(ind);
ncomponents = 0;
while nleft > 0
    ncomponents = ncomponents+1;
    component = false(dims);
    ind_component = find(still_left==1,1);
    still_left(ind_component) = 0;
    component(ind_component) = true;
    change = 1;
    while true
        [I,J] = ind2sub(dims,ind_component);
        [indN,indNE,indE,indSE,indS,indSW,indW,indNW] = find_ind_neighbors2d(I,J,dims);
        ind_neighbors = [indN;indNE;indE;indSE;indS;indSW;indW;indNW];
        ind_component = ind_neighbors(grain(ind_neighbors)>0);
        ind_component = unique(ind_component(not(component(ind_component))));
        component(ind_component) = true;
        change = change+length(ind_component);
        still_left(ind_component) = 0;
        if change == 0
            grains{ncomponents} = find(component);
            break
        else
            nleft = nleft-change;
            change = 0;
        end
    end
end

end