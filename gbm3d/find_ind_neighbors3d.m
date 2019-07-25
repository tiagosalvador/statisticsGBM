function [indNup,indNEup,indEup,indSEup,indSup,indSWup,indWup,indNWup,indup,...
        indNdown,indNEdown,indEdown,indSEdown,indSdown,indSWdown,indWdown,indNWdown,inddown,...
        indN,indNE,indE,indSE,indS,indSW,indW,indNW] = find_ind_neighbors3d(I,J,K,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [indNup,indNEup,indEup,indSEup,indSup,indSWup,indWup,indNWup,indup,...
%        indNdown,indNEdown,indEdown,indSEdown,indSdown,indSWdown,indWdown,indNWdown,inddown,...
%        indN,indNE,indE,indSE,indS,indSW,indW,indNW] = find_ind_neighbors3d(I,J,K,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns the linear indices of the 26 neighbours of (I,J,K).
% I, J and K can be vectors.
m = dims(1); n = dims(2); p = dims(3);

Iplus = I+1;
Iplus(Iplus>m) = 1;
Iminus = I-1;
Iminus(Iminus<1) = m;
Jplus = J+1;
Jplus(Jplus>n) = 1;
Jminus = J-1;
Jminus(Jminus<1) = n;
Kplus = K+1;
Kplus(Kplus>p) = 1;
Kminus = K-1;
Kminus(Kminus<1) = p;

indNup = sub2ind(dims,Iminus,J,Kplus);
indNEup = sub2ind(dims,Iminus,Jplus,Kplus);
indEup = sub2ind(dims,I,Jplus,Kplus);
indSEup = sub2ind(dims,Iplus,Jplus,Kplus);
indSup = sub2ind(dims,Iplus,J,Kplus);
indSWup = sub2ind(dims,Iplus,Jminus,Kplus);
indWup = sub2ind(dims,I,Jminus,Kplus);
indNWup = sub2ind(dims,Iminus,Jminus,Kplus);
indup = sub2ind(dims,I,J,Kplus);

indNdown = sub2ind(dims,Iminus,J,Kminus);
indNEdown = sub2ind(dims,Iminus,Jplus,Kminus);
indEdown = sub2ind(dims,I,Jplus,Kminus);
indSEdown = sub2ind(dims,Iplus,Jplus,Kminus);
indSdown = sub2ind(dims,Iplus,J,Kminus);
indSWdown = sub2ind(dims,Iplus,Jminus,Kminus);
indWdown = sub2ind(dims,I,Jminus,Kminus);
indNWdown = sub2ind(dims,Iminus,Jminus,Kminus);
inddown = sub2ind(dims,I,J,Kminus);


indN = sub2ind(dims,Iminus,J,K);
indNE = sub2ind(dims,Iminus,Jplus,K);
indE = sub2ind(dims,I,Jplus,K);
indSE = sub2ind(dims,Iplus,Jplus,K);
indS = sub2ind(dims,Iplus,J,K);
indSW = sub2ind(dims,Iplus,Jminus,K);
indW = sub2ind(dims,I,Jminus,K);
indNW = sub2ind(dims,Iminus,Jminus,K);
end