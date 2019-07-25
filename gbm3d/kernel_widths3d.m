function [alpha,beta,Svector] = kernel_widths3d(ori,angBrandon,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [alpha,beta,Svector] = kernel_widths3d(ori,angBrandon,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the kernel widths to use in the Algorithm 2 from
%   Salvador, T. & Esedoglu, S. The Role of Surface Tension and
%   Mobility Model in Simulations of Grain Growth.
%
%   as described in:
%
%   Salvador, T. & Esedoglu, S. J Sci Comput (2019) 79: 648.
%   https://doi.org/10.1007/s10915-018-0866-8.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = surface_tension3d(ori,angBrandon);

N = size(S,1);

switch option
    case 1
        % All mobilities = 1
        % M = ones(N)-eye(N); no need to create M since it's never needed
    case 2
        % All mobilities = 1/(surface tension)
        M = 1./S;
        M(1:(N+1):end) = 0;
end

%% Computing J*S*J


start_progress(' - Computing J*S*J')
Se = S*ones(N,1);
s = sum(S(:));
S = S+1/N^2*s;
for i = 1:N
    S(i,:) = S(i,:)-1/N*Se(i);
    S(:,i) = S(:,i)-1/N*Se(i);
    display_progress(i,N,1);
end

start_progress(' - Computing eigenvalues of J*S*J')
eigS = eig(S);
[~,i] = min(abs(eigS));
eigS(i) = [];

mineigS = min(eigS);
maxeigS = max(eigS);

display_progress(1,1,1);

start_progress(' - Reverting J*S*J')
for i = 1:N
    S(i,:) = S(i,:)+1/N*Se(i);
    S(:,i) = S(:,i)+1/N*Se(i);
    display_progress(i,N,1);
end
S = S-1/N^2*s;


% This alternative is faster but it doesn't converge all the time

% mineigS = eigs(S,1); % min(eigS)
% maxeigS = eigs(S,2,'smallestab'); % max(eigS)
% maxeigS = maxeigS(2);

if option == 1
    eigreciprocalM = -ones(N-1,1);
    maxeigreciprocalM = -1;
    mineigreciprocalM = -1;
else
    reciprocalM = M.^-1;
    reciprocalM((1:N+1:end)) = 0;
    eigreciprocalM = eig(J*reciprocalM*J);
    [~,i] = min(abs(eigreciprocalM));
    eigreciprocalM(i) = [];
    maxeigreciprocalM = max(eigreciprocalM);
    mineigreciprocalM = min(eigreciprocalM);
end

alpha = mineigS/maxeigreciprocalM;
beta = maxeigS/mineigreciprocalM;

if option==1
    aux = S;
else
    aux = S.*M;
end
aux(1:N+1:N^2) = -Inf;
alpha = max(alpha,max(aux(:)));
aux(1:N+1:N^2) = Inf;
beta = min(beta,min(aux(:)));
    
if sum(eigS > 0)
%if maxeigS > 0
    error('The matrix \sigma is not conditionally negative semidefinite.')
end
if sum(eigreciprocalM > 0)
%if maxeigreciprocalM > 0
    error('The matrix 1/\mu is not conditionally negative semidefinite.')
end

% In the special case where the surface tensions and the reciprocal
% mobilities are the same, alpha and beta will be the same. In such even,
% using the simplified algorithm with two Gaussian kernels is not
% recommended as one Gaussian kernel is indeed enough. Nonetheless the
% algorithm works provided we chose \alpha > \beta.
if alpha == beta
    alpha = 1.1*alpha;
    beta = 0.9*beta;
end

Svector = vectorizingS3d(S);

end
