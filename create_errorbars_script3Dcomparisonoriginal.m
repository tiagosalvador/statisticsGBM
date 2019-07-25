clear

set(0,'defaultaxeslinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultpatchlinewidth',1)

load('data3d-1')
nSim = size(volumes_grains,1);
dims = [256 256 256];
max_epoch = 2;

folder_name = './Histograms3Dcomparisonoriginal/';
file_name = ['simulation3dcomparison-',num2str(choice),'-',num2str(dims(1)),'x',num2str(dims(2)),'x',num2str(dims(3))];

for epoch = 0:max_epoch
    volumes_grains_aux{epoch+1} = vertcat(volumes_grains{:,epoch+1});
    volumes_grains_alt_aux{epoch+1} = vertcat(volumes_grains_alt{:,epoch+1});
    volumes_convhull_aux{epoch+1} = vertcat(volumes_convhull{:,epoch+1});
    surface_areas_grains_aux{epoch+1} = vertcat(surface_areas_grains{:,epoch+1});
    mdf_aux{epoch+1} = vertcat(mdf{:,epoch+1});
end

for epoch = 0:max_epoch
    volumes_grains_original_aux{epoch+1} = vertcat(volumes_grains_original{:,epoch+1});
    volumes_grains_alt_original_aux{epoch+1} = vertcat(volumes_grains_alt_original{:,epoch+1});
    volumes_convhull_original_aux{epoch+1} = vertcat(volumes_convhull_original{:,epoch+1});
    surface_areas_grains_original_aux{epoch+1} = vertcat(surface_areas_grains_original{:,epoch+1});
    mdf_original_aux{epoch+1} = vertcat(mdf_original{:,epoch+1});
end

load('data3doriginal-1')

for epoch = 0:max_epoch
    volumes_grains_original_aux2{epoch+1} = vertcat(volumes_grains_original{:,epoch+1});
    volumes_grains_alt_original_aux2{epoch+1} = vertcat(volumes_grains_alt_original{:,epoch+1});
    volumes_convhull_original_aux2{epoch+1} = vertcat(volumes_convhull_original{:,epoch+1});
    surface_areas_grains_original_aux2{epoch+1} = vertcat(surface_areas_grains_original{:,epoch+1});
    mdf_original_aux2{epoch+1} = vertcat(mdf_original{:,epoch+1});
end

    load('data3dconstant')
    
    for epoch = 0:max_epoch
        volumes_grains_constant_aux{epoch+1} = vertcat(volumes_grains{:,epoch+1});
        volumes_grains_alt_constant_aux{epoch+1} = vertcat(volumes_grains_alt{:,epoch+1});
        volumes_convhull_constant_aux{epoch+1} = vertcat(volumes_convhull{:,epoch+1});
        surface_areas_grains_constant_aux{epoch+1} = vertcat(surface_areas_grains{:,epoch+1});
        mdf_constant_aux{epoch+1} = vertcat(mdf{:,epoch+1});
    end

%% Setup Legend

legend_text{1} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = 1 (merging)';
legend_text{2} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = \sigma_{ij} (merging)';
legend_text{3} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = \sigma_{ij}';
legend_text{4} = '\sigma_{ij} = 1, \mu^{-1}_{ij} = 1';

%% Volumes

data{1} = volumes_grains_aux;
data{2} = volumes_grains_original_aux;
data{3} = volumes_grains_original_aux2;
data{4} = volumes_grains_constant_aux;

title_name = 'Probability Density vs Reduced Volume';
saveto = [folder_name,file_name,'-volumes_t'];
axis_bounds = [0 5 0 2];
scale = 1;
plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale);

%% Effective Radii

data{1} = cellfun(@(x)x.^(1/3),volumes_grains_aux,'UniformOutput',false);
data{2} = cellfun(@(x)x.^(1/3),volumes_grains_original_aux,'UniformOutput',false);
data{3} = cellfun(@(x)x.^(1/3),volumes_grains_original_aux2,'UniformOutput',false);
data{4} = cellfun(@(x)x.^(1/3),volumes_grains_constant_aux,'UniformOutput',false);
title_name = 'Probability Density vs Reduced Effective Radii';
saveto = [folder_name,file_name,'-radii_t'];
axis_bounds = [0 3 0 3];
scale = 1;
bin_width = [0.1 0.1 0.1];
plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);

%% Surface Areas

data{1} = surface_areas_grains_aux;
data{2} = surface_areas_grains_original_aux;
data{3} = surface_areas_grains_original_aux2;
data{4} = surface_areas_grains_constant_aux;
title_name = 'Probability Density vs Reduced Surface Areas';
saveto = [folder_name,file_name,'-surface_areas_t'];
axis_bounds = [0 4 0 1.6];
scale = 1;
plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale);

%% Isoperimetric Ratio

data{1} = cellfun(@(V,S) 36*pi*V.^2./S.^3,volumes_grains_aux,surface_areas_grains_aux,'UniformOutput',false);
data{2} = cellfun(@(V,S) 36*pi*V.^2./S.^3,volumes_grains_original_aux,surface_areas_grains_original_aux,'UniformOutput',false);
data{3} = cellfun(@(V,S) 36*pi*V.^2./S.^3,volumes_grains_original_aux2,surface_areas_grains_original_aux2,'UniformOutput',false);
data{4} = cellfun(@(V,S) 36*pi*V.^2./S.^3,volumes_grains_constant_aux,surface_areas_grains_constant_aux,'UniformOutput',false);
title_name = 'Probability Density vs Reduced Isoperimetric Ratio';
saveto = [folder_name,file_name,'-isopratio_t'];
axis_bounds = [0 1 0 14];
scale = 0;
bin_width = [0.02 0.02 0.02];
plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width);

%% MDF

clear legend_text

legend_text{1} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = 1 (merging)';
legend_text{2} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = \sigma_{ij} (merging)';
legend_text{3} = '\sigma_{ij} = RS, \mu^{-1}_{ij} = \sigma_{ij}';
legend_text{4} = 'Mackenzie Function';

clear data
data{1} = mdf_aux;
data{2} = mdf_original_aux;
data{3} = mdf_original_aux2;


BinWidthDegrees = 2;
NumBins = 360*1/BinWidthDegrees;
axis_bounds = [0 63 0 0.05];
title_name = 'Misorientation Distribution Function';
saveto = [folder_name,file_name,'-MDF_t'];
plot_comparison_mdf(data,NumBins,title_name,saveto,axis_bounds,legend_text);

