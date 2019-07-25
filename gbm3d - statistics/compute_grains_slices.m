function [grains_section, ID_section] = compute_grains_slices(grains,dims,dims_section,s,i,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [grains_section, ID_section] = compute_grains_slices(grains,dims,dims_section,s,i,Rfamilies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the 2D slices of a 3D grain into a 2D grain.

N = size(grains,1);
grains_section = cell(size(grains));
for k = 1:N
    ind = grains{k,1};
    [I,J,Z] = ind2sub(dims,ind);
    switch i
        case 1
            select = I == s;
            grains_section{k,1} = sub2ind(dims_section,J(select),Z(select));
        case 2
            select = J == s;
            grains_section{k,1} = sub2ind(dims_section,I(select),Z(select));
        case 3
            select = Z == s;
            grains_section{k,1} = sub2ind(dims_section,I(select),J(select));
    end
    grains_section{k,2} = grains{k,2}(select);
    grains_section{k,3} = grains{k,3}(select);
end
ID_section = (1:1:size(grains,1))';
% Remove empty grains
[grains_section,ID_section] = removeemptygrainsoriginal(grains_section,ID_section);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UPDATE CONVOLUTION VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For fast convolution, far away grains are grouped into families.
% Then, the entire family (union of grains) is convolved at once.
[families,famgrains] = grains2families2doriginal(grains_section,Rfamilies,dims_section,0); % Group grains into families.
N = size(families,1); % Number of families.
for k=1:N % Loop over families.
    cval = convolvefamily2doriginal(families{k,1},families{k,2},dims_section);
    % Distribute convolution values to the grains contained in this family:
    numberofgrains = size(famgrains{k},1); % Number of grains in this family.
    listofgrains = famgrains{k};           % Column vector of grain indices.
    for ell = 1:numberofgrains % Loop over grains contained in this family.
        label = listofgrains(ell);
        ind = grains_section{label,1};
        grains_section{label,3} = cval(ind); % Read off and record the convolution vals.
    end % (for ell) Loop over grains ends.
end % (for k) Loop over families ends.
end