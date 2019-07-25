function [parameters,ori_matrix] = gradient_descent_aux(density_target,parameters,ori_matrix,lengths,max_iter,tolerance)

setup_rot

% Setup grid to approximate the L2 energy function
Ngrid = 100;
wmax = 2*atan(sqrt(23-16*sqrt(2)));
x = linspace(0,wmax,Ngrid);
dx = x(2)-x(1);
% Width of the bump function
delta = 0.02;

% Number of grains
N = size(ori_matrix,1)/3;

% Time step gradient descent
dt = .5;

% Determining misorientation
[theta,~,index] = misorientation3d_grad(ori_matrix,N);

% Since most of the lengths are zero (i.e. most pairs of grains
% do not share an interface), it is very efficient to do the
% computations only on the nonzero ones.
active = lengths>0;

E = energy(x,density_target,@density_estimate,lengths(active),theta(active),delta);

grad_E_x1 = zeros(1,N);
grad_E_x2 = zeros(1,N);
grad_E_x3 = zeros(1,N);

theta_old = theta;
data_x = x;
data_density_estimate_old = density_estimate(x,lengths(active),theta_old(active),delta);

grad_theta_ij = zeros(size(lengths));

t = 0;
while and(t<max_iter,E>tolerance)
    t = t+1;

    % needed for Matlab version on screms1
    [a,b] = meshgrid(x',theta(active));
    aux = (bump_x(a-b,delta))';
    grad_theta_ij(active) = 2*((density_target(x)-density_estimate(x,lengths(active),theta(active),delta))*aux).*lengths(active)*dx;
    
    % grad_theta_ij(active) = 2*((density_target(x)-density_estimate(x,lengths(active),theta(active),delta))*bump_x(x'-theta(active),delta)).*lengths(active)*dx;
    
    for i=1:N
        x1 = parameters(i,1);
        x2 = parameters(i,2);
        x3 = parameters(i,3);
        [dgidx1,dgidx2,dgidx3] = get_dgi(x1,x2,x3);
        j = [1:i-1 i+1:N];
        [dthetaijdx1,dthetaijdx2,dthetaijdx3] = get_dthetaij(i,j,theta,ori_matrix,dgidx1,dgidx2,dgidx3,index,rot);
        iaux = min(i,j);
        jaux = max(i,j);
        grad_E_x1(i) = sum(grad_theta_ij(jaux-1/2*iaux.*(1+iaux-2*N)-N).*dthetaijdx1);
        grad_E_x2(i) = sum(grad_theta_ij(jaux-1/2*iaux.*(1+iaux-2*N)-N).*dthetaijdx2);
        grad_E_x3(i) = sum(grad_theta_ij(jaux-1/2*iaux.*(1+iaux-2*N)-N).*dthetaijdx3);
    end

    for k=1:N
        parameters(k,:) = min(max(parameters(k,:)-dt*[grad_E_x1(k) grad_E_x2(k) grad_E_x3(k)],0),1);
        ori_matrix((3*k-2):3*k,:) = SO3(parameters(k,1),parameters(k,2),parameters(k,3));
    end
    
    [theta,~,index] = misorientation3d_grad(ori_matrix,N);
    E = energy(x,density_target,@density_estimate,lengths(active),theta(active),delta);
    c = clock;
    fprintf('E = %1.5f at iteration %d %02d:%02d:%02d\n',E,t,c(4),c(5),round(c(6)));
    
    %% Control data
    
    data_density_estimate = density_estimate(x,lengths(active),theta(active),delta);
    
    save('data_gradient_descent','data_x','data_density_estimate_old','data_density_estimate')

end
dispstat('Finished.','keepprev');
end
