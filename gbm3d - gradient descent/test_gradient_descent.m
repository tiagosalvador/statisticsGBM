clear
setup_rot

% Control seed
rng(2)
max_angle_disorientation = 2*atan(sqrt(23-16*sqrt(2)));

% Setup grid to approximate the L2 energy function
Ngrid = 100;
x = linspace(0,max_angle_disorientation,Ngrid);
dx = x(2)-x(1);
delta = 0.1;

% Target distribution
density_target = @(x) mackenzie_function_perturbed(x);

% Number of grains
N = 500;

% Generating initial misorientation distribution
parameters = zeros(N,3);
ori = cell(1,N);
for k=1:N
    parameters(k,:) = [rand rand rand];
    ori{k} = SO3(parameters(k,1),parameters(k,2),parameters(k,3));
end

% Determining misorientation
[theta,minang,index] = misorientation3d_grad(ori,N);
l = rand(size(theta));
l(300:end) = 0;
l = l/sum(l);

E = energy(x,density_target,@density_estimate,l,theta,delta);

theta_old = theta;
ori_old = ori;
Eold = E;


% figure(1);
% plot(x,density_estimate(x,l,theta_old,delta),x,density_estimate(x,l,theta,delta))
% figure(2);
% histogram(theta_old)
% hold on
% histogram(theta)
% hold off

dthetaijdx1 = zeros(N,N);
dthetaijdx2 = zeros(N,N);
dthetaijdx3 = zeros(N,N);

tic
profile off
profile on

% tic
% for i=1:N
%     x1 = parameters(i,1);
%     x2 = parameters(i,2);
%     x3 = parameters(i,3);
%     dgidx1 = get_dgidx1(x1,x2,x3);
%     dgidx2 = get_dgidx2(x1,x2,x3);
%     dgidx3 = get_dgidx3(x1,x2,x3);
%     
%     j = setdiff(1:N,i);
%     [dthetaijdx1_vec(i,j),dthetaijdx2_vec(i,j),dthetaijdx3_vec(i,j)] = get_dthetaij_vec(i,j,theta,ori,dgidx1,dgidx2,dgidx3,index,rot);
% 
% end
% toc
% 
% tic
% for i=1:N
%     x1 = parameters(i,1);
%     x2 = parameters(i,2);
%     x3 = parameters(i,3);
%     dgidx1 = get_dgidx1(x1,x2,x3);
%     dgidx2 = get_dgidx2(x1,x2,x3);
%     dgidx3 = get_dgidx3(x1,x2,x3);
%     
%     for j=1:N
%         if not(i==j)
%             [dthetaijdx1(i,j),dthetaijdx2(i,j),dthetaijdx3(i,j)] = get_dthetaij(i,j,theta,ori,dgidx1,dgidx2,dgidx3,index,rot);
%         end
%     end
% end
% toc
% profile viewer

dt = .1;
for t=1:5000
    Eold = E;
    ori_old = ori;
    
    for i=1:N
        x1 = parameters(i,1);
        x2 = parameters(i,2);
        x3 = parameters(i,3);
        dgidx1 = get_dgidx1(x1,x2,x3);
        dgidx2 = get_dgidx2(x1,x2,x3);
        dgidx3 = get_dgidx3(x1,x2,x3);
        
        j = setdiff(1:N,i);
        [dthetaijdx1(i,j),dthetaijdx2(i,j),dthetaijdx3(i,j)] = get_dthetaij_vec(i,j,theta,ori,dgidx1,dgidx2,dgidx3,index,rot);
        
    end
%     for i=1:N
%         x1 = parameters(i,1);
%         x2 = parameters(i,2);
%         x3 = parameters(i,3);
%         dgidx1 = get_dgidx1(x1,x2,x3);
%         dgidx2 = get_dgidx2(x1,x2,x3);
%         dgidx3 = get_dgidx3(x1,x2,x3);
% 
%         for j=1:N
%             if not(i==j)
%                 [dthetaijdx1(i,j),dthetaijdx2(i,j),dthetaijdx3(i,j)] = get_dthetaij(i,j,theta,ori,dgidx1,dgidx2,dgidx3,index,rot);
%             end
%         end
%     end
    
    grad_theta_ij = 2*((density_target(x)-density_estimate(x,l,theta,delta))*bump_x(x'-theta,delta)).*l*dx;
    grad_E_x1 = zeros(size(ori));
    grad_E_x2 = zeros(size(ori));
    grad_E_x3 = zeros(size(ori));
    for i=1:N
        for j=1:N
            if not(i==j)
                iaux = min(i,j);
                jaux = max(i,j);
                grad_E_x1(i) = grad_E_x1(i)+grad_theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N)*dthetaijdx1(i,j);
                grad_E_x2(i) = grad_E_x2(i)+grad_theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N)*dthetaijdx2(i,j);
                grad_E_x3(i) = grad_E_x3(i)+grad_theta_ij(jaux-1/2*iaux*(1+iaux-2*N)-N)*dthetaijdx3(i,j);
            end
        end
    end
    
%     dxx = 0.0001;
%     for i=1:N
%         ori = ori_old;
%         ori{i} = SO3(parameters(i,1)+dxx,parameters(i,2),parameters(i,3));
%         [theta,~,~] = misorientation3d_grad(ori,N);
%         E = energy(x,density_target,@density_estimate,l,theta,delta);
%         grad_E_x1_approx(i) = (E-Eold)/dxx;
%     end
%     max(abs(grad_E_x1-grad_E_x1_approx))
%     for i=1:N
%         ori = ori_old;
%         ori{i} = SO3(parameters(i,1),parameters(i,2)+dxx,parameters(i,3));
%         [theta,~,~] = misorientation3d_grad(ori,N);
%         E = energy(x,density_target,@density_estimate,l,theta,delta);
%         grad_E_x2_approx(i) = (E-Eold)/dxx;
%     end
%     max(abs(grad_E_x2-grad_E_x2_approx))
%     for i=1:N
%         ori = ori_old;
%         ori{i} = SO3(parameters(i,1),parameters(i,2),parameters(i,3)+dxx);
%         [theta,~,~] = misorientation3d_grad(ori,N);
%         E = energy(x,density_target,@density_estimate,l,theta,delta);
%         grad_E_x3_approx(i) = (E-Eold)/dxx;
%     end
%     max(abs(grad_E_x3-grad_E_x3_approx))
    
    for k=1:N
        %parameters{k} = min(max(parameters{k}-dt*[grad_E_x1_approx(k);grad_E_x2_approx(k);grad_E_x3_approx(k)],0),1);
        parameters(k,:) = min(max(parameters(k,:)-dt*[grad_E_x1(k) grad_E_x2(k) grad_E_x3(k)],10^-3),1);
        ori{k} = SO3(parameters(k,1),parameters(k,2),parameters(k,3));
    end
    
    
    [theta,~,index] = misorientation3d_grad(ori,N);
    E = energy(x,density_target,@density_estimate,l,theta,delta);
    
    if mod(t,10)==0
        disp(t)
        E = energy(x,density_target,@density_estimate,l,theta,delta)
        figure(1);
        plot(x,density_estimate(x,l,theta_old,delta),x,density_estimate(x,l,theta,delta),x,density_target(x))
        figure(2);
        histogram(theta_old,50,'Normalization','probability')
        hold on
        histogram(theta,50,'Normalization','probability')
        hold off
        pause(0.001)
    end
    
end
profile viewer
toc
%clf

%plot(x,density_estimate(x,l,theta_old,delta),x,density_estimate(x,l,theta,delta))