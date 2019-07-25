clear
setup_rot

% Control seed
rng(1)

% Setup grid to approximate the L2 energy function
Ngrid = 100;
x = linspace(0,63/360*2*pi,Ngrid);
dx = x(2)-x(1);
delta = 0.1;

% Target distribution
density_target = @(x) 1;

% Number of grains
N = 3;

% Generating initial misorientation distribution
parameters = cell(1,N);
ori = cell(1,N);
for k=1:N
    parameters{k} = [rand;rand;rand];
    ori{k} = SO3(parameters{k}(1),parameters{k}(2),parameters{k}(3));
end

% Determining misorientation
[theta,minang,index] = misorientation3d_grad(ori,N);
l = 1/length(theta)*ones(size(theta));

E = energy(x,density_target,@density_estimate,l,theta,delta);

exact = 2*(density_target(x)-density_estimate(x,l,theta,delta))*bump_x(x'-theta,delta).*l*dx;

h = 10^(-6);

theta_perturbed = theta;
theta_perturbed(1) = theta_perturbed(1)+h;
E_perturbed = energy(x,density_target,@density_estimate,l,theta_perturbed,delta);
approx = (E_perturbed-E)/h;
norm(exact(1)-approx)

theta_perturbed = theta;
theta_perturbed(2) = theta_perturbed(2)+h;
E_perturbed = energy(x,density_target,@density_estimate,l,theta_perturbed,delta);
approx = (E_perturbed-E)/h;
norm(exact(2)-approx)

theta_perturbed = theta;
theta_perturbed(3) = theta_perturbed(3)+h;
E_perturbed = energy(x,density_target,@density_estimate,l,theta_perturbed,delta);
approx = (E_perturbed-E)/h;
norm(exact(3)-approx)

