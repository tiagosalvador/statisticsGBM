clear

set(0,'defaultaxeslinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultpatchlinewidth',1)

for choice = [1]
    load(['data3d-',num2str(choice),'-slices'])
    nSim = size(areas_grains_slices,1);
    dims = [256 256 256];
    max_epoch = 2;
    
    folder_name = './Histograms3Dcomparison/';
    file_name = ['simulation3dcomparison-',num2str(choice),'-',num2str(dims(1)),'x',num2str(dims(2)),'x',num2str(dims(3))];
    
    for k = 1:nSim
        for epoch = 0:max_epoch
            for i = 1:length(dims)
                areas_aux1{k,epoch+1,i} = vertcat(areas_grains_slices{k,epoch+1}{i,:});
                isopratio_aux1{k,epoch+1,i} = vertcat(isopratio_slices{k,epoch+1}{i,:});
            end
        end
    end
    for epoch = 0:max_epoch
        for i = 1:length(dims)
            areas_grains_aux{epoch+1,i} = vertcat(areas_aux1{:,epoch+1,i});
            isopratio_aux{epoch+1,i} = vertcat(isopratio_aux1{:,epoch+1,i});
        end
    end
    for k = 1:nSim
        for epoch = 0:max_epoch
            for i = 1:length(dims)
                areas_aux1{k,epoch+1,i} = vertcat(areas_grains_slices_original{k,epoch+1}{i,:});
                isopratio_aux1{k,epoch+1,i} = vertcat(isopratio_slices_original{k,epoch+1}{i,:});
            end
        end
    end
    for epoch = 0:max_epoch
        for i = 1:length(dims)
            areas_grains_original_aux{epoch+1,i} = vertcat(areas_aux1{:,epoch+1,i});
            isopratio_original_aux{epoch+1,i} = vertcat(isopratio_aux1{:,epoch+1,i});
        end
    end
    
    
    load('data3dconstant-slices')
    
    for k = 1:nSim
        for epoch = 0:max_epoch
            for i = 1:length(dims)
                areas_aux1{k,epoch+1,i} = vertcat(areas_grains_slices{k,epoch+1}{i,:});
                isopratio_aux1{k,epoch+1,i} = vertcat(isopratio_slices{k,epoch+1}{i,:});
            end
        end
    end
    for epoch = 0:max_epoch
        for i = 1:length(dims)
            areas_grains_constant_aux{epoch+1,i} = vertcat(areas_aux1{:,epoch+1,i});
            isopratio_constant_aux{epoch+1,i} = vertcat(isopratio_aux1{:,epoch+1,i});
        end
    end
    
    %% Setup Legend
    
    legend_text{1} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = 1';
    legend_text{2} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = \sigma_{ij}';
    legend_text{3} = '\sigma_{ij} = 1, \mu^{-1}_{ij} = 1';
    
    %% Slices
    
    for i = 3:3
        %% Effective Radii

        data{1} = cellfun(@sqrt,areas_grains_aux(:,i)','UniformOutput',false);
        data{2} = cellfun(@sqrt,areas_grains_original_aux(:,i)','UniformOutput',false);
        data{3} = cellfun(@sqrt,areas_grains_constant_aux(:,i)','UniformOutput',false);
        
        title_name = 'Probability Density vs Reduced Effective Radii';
        saveto = [folder_name,file_name,'-radii_slices_',num2str(i),'_t'];
        axis_bounds = [0 3 0 1.1];
        scale = 1;
        bin_width = [0.1 0.1 0.1];
        [xradii3DX,yradii3DX,eradii3DX] = plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);

        save([folder_name,'radii3DX-choice',num2str(choice),'-slice',num2str(i)],'xradii3DX','yradii3DX','eradii3DX')
        
        %% Isoperimetric ratio
        
        data{1} = isopratio_aux(:,i)';
        data{2} = isopratio_original_aux(:,i)';
        data{3} = isopratio_constant_aux(:,i)';
        title_name = 'Probability Density vs Reduced Isoperimetric Ratio';
        saveto = [folder_name,file_name,'-isopratio_slices_',num2str(i),'_t'];
        axis_bounds = [0 1 0 11];
        scale = 1;
        bin_width = [0.01 0.01 0.01];
        plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);
 
    end
end

%%

clear
close all

%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENTAL DATA:
%%%%%%%%%%%%%%%%%%%%%%
load hull_3dx_experimental.mat

% Normalization for Hull's data:
area = 0;
for i=1:10
  add = ( y(i+1)+y(i) ) / 2 * ( x(i+1) - x(i) );
  area = area + add;
end
y = y / area;

moment = 0;
for i=1:10
  add = ( y(i+1)*x(i+1)+y(i)*x(i) ) / 2 * ( x(i+1) - x(i) );
  moment = moment + add;
end
x = x / moment;

area = 0;
for i=1:10
  add = ( y(i+1)+y(i) ) / 2 * ( x(i+1) - x(i) );
  area = area + add;
end
y = y / area;

xHull = x;
yHull = y;

load groeber_3dx_experimental.mat

% Normalization for Groeber's data:

area = 0;
for i=1:9
  add = ( y(i+1)+y(i) ) / 2 * ( x(i+1) - x(i) );
  area = area + add;
end
y = y / area;

moment = 0;
for i=1:9
  add = ( y(i+1)*x(i+1)+y(i)*x(i) ) / 2 * ( x(i+1) - x(i) );
  moment = moment + add;
end
x = x / moment;

area = 0;
for i=1:9
  add = ( y(i+1)+y(i) ) / 2 * ( x(i+1) - x(i) );
  area = area + add;
end
y = y / area;

xGroeber = x;
yGroeber = y;


load('radii3DX-choice1-slice3.mat')

dims = [256 256 256];
choice = 1;

folder_name = './Histograms3Dcomparison/';
file_name = ['simulation3dcomparison-',num2str(choice),'-',num2str(dims(1)),'x',num2str(dims(2)),'x',num2str(dims(3))];

title_name = 'Probability Density vs Reduced Effective Radii';
saveto = [folder_name,file_name,'-radii_slices_3_t'];
axis_bounds = [0 3 0 1];
scale = 1;
legend_text{1} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = 1';
legend_text{2} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = \sigma_{ij}';
legend_text{3} = '\sigma_{ij} = 1, \mu^{-1}_{ij} = 1';
legend_text{4} = 'Hull et. al.';
legend_text{5} = 'Groeber et. al.';

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
        errorbar(xradii3DX{i,j},yradii3DX{i,j},eradii3DX{i,j},marker,'LineWidth',1.5)
    end
    plot(xHull,yHull,'-o',xGroeber,yGroeber,'-d')
    hold off
    legend(legend_text,'Location','best')
    axis(axis_bounds)
    title(title_name)
    print(gcf, '-depsc2', [saveto,num2str(j-1),'-experiments']);
end

