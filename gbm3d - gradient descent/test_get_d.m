clear
rng(1)
setup_rot


x1 = rand;
x2 = rand;
x3 = rand;
gi = SO3(x1,x2,x3);

h = 10^(-5);

% Test get_dRdx1

exact = get_dRdx1(x1,x2,x3);
approx = (get_R(x1+h,x2,x3)-get_R(x1,x2,x3))/h;
norm(exact-approx)

% Test get_dvdx2

exact = get_dvdx2(x1,x2,x3);
approx = (get_v(x1,x2+h,x3)-get_v(x1,x2,x3))/h;
norm(exact-approx)

% Test get_dvdx3

h = 10^(-10);

exact = get_dvdx3(x1,x2,x3);
approx = (get_v(x1,x2,x3+h)-get_v(x1,x2,x3))/h;
norm(exact-approx)

% Test get_dgidx1

h = 10^(-8);

exact = get_dgidx1(x1,x2,x3);
approx = (get_gi(x1+h,x2,x3)-get_gi(x1,x2,x3))/h;
norm(exact-approx)

% Test get_dgidx2

h = 10^(-8);

exact = get_dgidx2(x1,x2,x3);
approx = (get_gi(x1,x2+h,x3)-get_gi(x1,x2,x3))/h;
norm(exact-approx)

% Test get_dgidx3

h = 10^(-10);

exact = get_dgidx3(x1,x2,x3);
approx = (get_gi(x1,x2,x3+h)-get_gi(x1,x2,x3))/h;
norm(exact-approx)

% Test get_dthetaijdx1, get_dthetaijdx2, get_dthetaijdx3

h = 10^(-6);

N = 3;
parameters = zeros(N,3);
ori = cell(1,N);
for k=1:N
    parameters(k,:) = [rand rand rand];
    ori{k} = SO3(parameters(k,1),parameters(k,2),parameters(k,3));
end
[theta_ij,~,index] = misorientation3d_grad(ori,N);


exact1 = zeros(N,N);
exact2 = zeros(N,N);
exact3 = zeros(N,N);
approx1 = zeros(N,N);
approx2 = zeros(N,N);
approx3 = zeros(N,N);


exact = zeros(N,N);
approx = zeros(N,N);

% Vectorized
for i=1:N
    for j=1:N
        if not(i==j)
            x1 = parameters(i,1);
            x2 = parameters(i,2);
            x3 = parameters(i,3);
            gi = get_gi(x1,x2,x3);
            iaux = min(i,j);
            jaux = max(i,j);
            [exact1(i,j),exact2(i,j),exact3(i,j)] = get_dthetaij(i,j,ori,parameters,index,rot);
            gi_perturbed = get_gi(x1+h,x2,x3);
            ori_perturbed = ori;
            ori_perturbed{i} = gi_perturbed;
            [theta_ij_perturbed,~,~] = misorientation3d_grad(ori_perturbed,N);
            approx1(i,j) = (theta_ij_perturbed(jaux-1/2*iaux*(1+iaux-2*N)-N)-theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N))/h;
            
            gi_perturbed = get_gi(x1,x2+h,x3);
            ori_perturbed = ori;
            ori_perturbed{i} = gi_perturbed;
            [theta_ij_perturbed,~,~] = misorientation3d_grad(ori_perturbed,N);
            approx2(i,j) = (theta_ij_perturbed(jaux-1/2*iaux*(1+iaux-2*N)-N)-theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N))/h;
            
            gi_perturbed = get_gi(x1,x2,x3+h);
            ori_perturbed = ori;
            ori_perturbed{i} = gi_perturbed;
            [theta_ij_perturbed,~,~] = misorientation3d_grad(ori_perturbed,N);
            approx3(i,j) = (theta_ij_perturbed(jaux-1/2*iaux*(1+iaux-2*N)-N)-theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N))/h;
        end
    end
end
norm(exact1-approx1)
norm(exact2-approx2)
norm(exact3-approx3)

for i=1:N
    for j=1:N
        if not(i==j)
            x1 = parameters(i,1);
            x2 = parameters(i,2);
            x3 = parameters(i,3);
            gi = get_gi(x1,x2,x3);
            iaux = min(i,j);
            jaux = max(i,j);
            exact(i,j) = get_dthetaijdx1(i,j,ori,parameters,index,rot);
            gi_perturbed = get_gi(x1+h,x2,x3);
            ori_perturbed = ori;
            ori_perturbed{i} = gi_perturbed;
            [theta_ij_perturbed,~,~] = misorientation3d_grad(ori_perturbed,N);
            approx(i,j) = (theta_ij_perturbed(jaux-1/2*iaux*(1+iaux-2*N)-N)-theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N))/h;
        end
    end
end
norm(exact-approx)

for i=1:N
    for j=1:N
        if not(i==j)
            x1 = parameters(i,1);
            x2 = parameters(i,2);
            x3 = parameters(i,3);
            gi = get_gi(x1,x2,x3);
            iaux = min(i,j);
            jaux = max(i,j);
            exact(i,j) = get_dthetaijdx2(i,j,ori,parameters,index,rot);
            gi_perturbed = get_gi(x1,x2+h,x3);
            ori_perturbed = ori;
            ori_perturbed{i} = gi_perturbed;
            [theta_ij_perturbed,~,~] = misorientation3d_grad(ori_perturbed,N);
            approx(i,j) = (theta_ij_perturbed(jaux-1/2*iaux*(1+iaux-2*N)-N)-theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N))/h;
        end
    end
end
norm(exact-approx)


for i=1:N
    for j=1:N
        if not(i==j)
            x1 = parameters(i,1);
            x2 = parameters(i,2);
            x3 = parameters(i,3);
            gi = get_gi(x1,x2,x3);
            iaux = min(i,j);
            jaux = max(i,j);
            exact(i,j) = get_dthetaijdx3(i,j,ori,parameters,index,rot);
            gi_perturbed = get_gi(x1,x2,x3+h);
            ori_perturbed = ori;
            ori_perturbed{i} = gi_perturbed;
            [theta_ij_perturbed,~,~] = misorientation3d_grad(ori_perturbed,N);
            approx(i,j) = (theta_ij_perturbed(jaux-1/2*iaux*(1+iaux-2*N)-N)-theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N))/h;
        end
    end
end
norm(exact-approx)
