clear

set(0,'defaultaxeslinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultpatchlinewidth',1)

load('data2d.mat')
nSim = size(areas_grains,1);
max_epoch = 2;
folder_name = './Histograms2Dcomparison/';
file_name = ['simulation2dcomparison-',num2str(dims(1)),'x',num2str(dims(2))];


for epoch = 0:max_epoch
    areas_convhull_aux{epoch+1} = vertcat(areas_convhull{:,epoch+1});
    areas_convhull_original_aux{epoch+1} = vertcat(areas_convhull_original{:,epoch+1});
    areas_grains_aux{epoch+1} = vertcat(areas_grains{:,epoch+1});
    areas_grains_original_aux{epoch+1} = vertcat(areas_grains_original{:,epoch+1});
    isopratio_aux{epoch+1} = vertcat(isopratio{:,epoch+1});
    isopratio_original_aux{epoch+1} = vertcat(isopratio_original{:,epoch+1});
    mdf_aux{epoch+1} = vertcat(mdf{:,epoch+1});
    mdf_original_aux{epoch+1} = vertcat(mdf_original{:,epoch+1});
end

load('data2dconstant')

for epoch = 0:max_epoch
    areas_grains_constant_aux{epoch+1} = vertcat(areas_grains{:,epoch+1});
    areas_convhull_constant_aux{epoch+1} = vertcat(areas_convhull{:,epoch+1});
    isopratio_constant_aux{epoch+1} = vertcat(isopratio{:,epoch+1});
end

%% Setup Legend

legend_text{1} = '\sigma_{ij} = RS, \mu_{ij} = 1';
legend_text{2} = '\sigma_{ij} = RS, \mu_{ij} = \sigma^{-1}_{ij}';
legend_text{3} = '\sigma_{ij} = 1, \mu_{ij} = 1';

%% Volumes

data{1} = areas_grains_aux;
data{2} = areas_grains_original_aux;
data{3} = areas_grains_constant_aux;

title_name = 'Probability Density vs Reduced Area';
saveto = [folder_name,file_name,'-areas_t'];
axis_bounds = [0 5 0 2];
scale = 1;
bin_width = [0.1 0.1 0.1];
[xGSD2D,yGSD2D] = plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);
save([folder_name,'GSD2D'],'xGSD2D','yGSD2D')

%% Effective Radii

data{1} = cellfun(@sqrt,areas_grains_aux,'UniformOutput',false);
data{2} = cellfun(@sqrt,areas_grains_original_aux,'UniformOutput',false);
data{3} = cellfun(@sqrt,areas_grains_constant_aux,'UniformOutput',false);
title_name = 'Probability Density vs Reduced Effective Radii';
saveto = [folder_name,file_name,'-radii_t'];
axis_bounds = [0 3 0 3];
scale = 1;
bin_width = [0.05 0.05 0.05];
[xradii2D,yradii2D] = plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);
save([folder_name,'radii2D'],'xradii2D','yradii2D')

%% Perimeters

data{1} = cellfun(@(A,I)sqrt(4*pi*A./I),areas_grains_aux,isopratio_aux,'UniformOutput',false);
data{2} = cellfun(@(A,I)sqrt(4*pi*A./I),areas_grains_original_aux,isopratio_original_aux,'UniformOutput',false);
data{3} = cellfun(@(A,I)sqrt(4*pi*A./I),areas_grains_constant_aux,isopratio_constant_aux,'UniformOutput',false);
title_name = 'Probability Density vs Reduced Perimeter';
saveto = [folder_name,file_name,'-perimeters_t'];
axis_bounds = [0 4 0 1.6];
scale = 1;
bin_width = [0.1 0.1 0.1];
plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);

%% Isoperimetric ratio

data{1} = isopratio_aux;
data{2} = isopratio_original_aux;
data{3} = isopratio_constant_aux;
title_name = 'Probability Density vs Isoperimetric Ratio';
saveto = [folder_name,file_name,'-isopratio_t'];
axis_bounds = [0.4 1 0 15];
scale = 0;
bin_width = [0.01 0.01 0.01];
[xisopratio2D,yisopratio2D] = plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);
save([folder_name,'isopratio2D'],'xisopratio2D','yisopratio2D')

%% Convexity Ratio

for i=1:length(areas_grains_aux)
    data{1}{i} = (areas_convhull_aux{i}-areas_grains_aux{i})./areas_grains_aux{i};
    data{2}{i} = (areas_convhull_original_aux{i}-areas_grains_original_aux{i})./areas_grains_original_aux{i};
    data{3}{i} = (areas_convhull_constant_aux{i}-areas_grains_constant_aux{i})./areas_grains_constant_aux{i};
end
title_name = 'Probability Density vs Convexity Ratio';
saveto = [folder_name,file_name,'-convexity_t'];
axis_bounds = [0 .1 0 100];
scale = 0;
bin_width = [0.01 0.01 0.01];
plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);

%% MDF

clear legend_text

legend_text{1} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = 1';
legend_text{2} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = \sigma_{ij}';

clear data
data{1} = mdf_aux;
data{2} = mdf_original_aux;

th_ang = 2;
NumBins = 180/th_ang;
axis_bounds = [0 180 0 0.06];
title_name = 'Misorientation Distribution Function';
saveto = [folder_name,file_name,'-MDF_t'];
plot_comparison_mdf(data,NumBins,title_name,saveto,axis_bounds,legend_text);