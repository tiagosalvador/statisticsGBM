clear

%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%

load('mason_3d.mat')
% Normalization for experimental data:
area = 0;
for i=1:22
  add = ( y(i+1)+y(i) ) / 2 * ( x(i+1) - x(i) );
  area = area + add;
end
y = y / area;

load('radii3D-choice1.mat')

dims = [256 256 256];
choice = 1;

folder_name = './Histograms3Dcomparison/';
file_name = ['simulation3dcomparison-',num2str(choice),'-',num2str(dims(1)),'x',num2str(dims(2)),'x',num2str(dims(3))];

title_name = 'Probability Density vs Reduced Effective Radii';
saveto = [folder_name,file_name,'-radii_t'];
axis_bounds = [0 3 0 3];
scale = 1;
legend_text{1} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = 1';
legend_text{2} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = \sigma_{ij}';
legend_text{3} = '\sigma_{ij} = 1, \mu^{-1}_{ij} = 1';
legend_text{4} = 'Mason et. al.';

for j = 2:3
    clf
    hold on
    for i = 1:3
        switch i
            case 1
                marker = '-';
            case 2
                marker = '-';
            case 3
                marker = '-';
        end
        errorbar(xradii3D{i,j},yradii3D{i,j},eradii3D{i,j},marker,'LineWidth',1.5)
    end
    plot(x,y,'-o')
    hold off
    legend(legend_text,'Location','best')
    axis(axis_bounds)
    title(title_name)
    print(gcf, '-depsc2', [saveto,num2str(j-1),'-simulations']);
end
