clear
%% Load the SMAP data
input_drive = 'Z:';
input_dir = [input_drive, '/IceSheets/data/Greenland Data/rSIR/SMAP/NP/essential_vars/v2/'];
out_dir = 'Z:\IceSheets\Documentation\Manuscripts_n_abstracts\MWA_journal\Figs\';
promice_dir = 'Z:\IceSheets\Weather Stations Data\Greenland\PROMICE\DATA\AWS\V9\combined\';
load([promice_dir,'/locations_promice_aws_revA.mat']')
load([input_drive, '/IceSheets/Weather Stations Data/Greenland/PROMICE/DATA/AWS/V10/SMAP_TB_FS_PROMICE/mean_SMAP_FS_TBVnTBH_at_PROMICE_AWS_2015-2023.mat']);
%% PROMICE Ice Mask
promicedir = [input_drive, '\IceSheets\data\Promice_ice_mask\'];
load([promicedir, 'iceedge_N03_125.mat']);
load("ice_mask_constructed_from_L1_9km_ice_mask.mat")
[lat_area,lon_area] = return_rSIR_lat_lon_subset_v1(input_drive);
%% Mean TB Maps
fs=14;
levels = [1000,1500,2000,2500,3000];
% figure('position', [5.7694e+03 -309.4000 1.7500e+03 850]);
t=tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
greenland('color', colo('grey1'));
% pcolorpsn(lat_area, lon_area, mtbv_refsp2.*ice_mask_totall2);
pcolorpsn(lat_area, lon_area, nanmean(mtbvly,3).*ice_mask2);
% [C,h]= contourpsn(Lat,Lon,tmp_elev,levels,'k--');
plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
set(gca, 'clim', [140 260]);
title(sprintf('%s','V-pol TB (TBV)'),'FontSize', fs);
xlim([-0.7e6 0.9e6]);
ylim([-3.5e6 -0.5e6 ]);
zoom on;
cmap=parula; %
% cmap=get_python_cmap('bwr'); % 
% cmap=get_python_cmap('plasma'); % plasma (perceptually uniform seq)
% cmap=get_python_cmap('gist_ncar');% Very well differentiating

colormap(cmap(:,1:3))
set(gca,'FontWeight','bold')
axis off

nexttile;
greenland('color', colo('grey1'));
% pcolorpsn(lat_area, lon_area, mtbh_refsp2.*ice_mask_totall2);
pcolorpsn(lat_area, lon_area, nanmean(mtbhly,3).*ice_mask2);
% [C,h]= contourpsn(Lat,Lon,tmp_elev,levels,'k--');
plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
set(gca, 'clim', [140 260]);
title(sprintf('%s','H-pol TB (TBH)'),'FontSize', fs);
xlim([-0.7e6 0.9e6]);
ylim([-3.5e6 -0.5e6 ]);
zoom on;
% cmap=get_python_cmap('RdYlBu_r');
colormap(cmap(:,1:3))
c1=colorbar;
% c1.Location = 'southoutside';
c1.Location = 'south';
set(c1, 'Position', [0.18 0.0348 0.3 0.0182]);
set(gca,'FontWeight','bold')
axis off


nexttile;
greenland('color', colo('grey1'));
% pcolorpsn(lat_area, lon_area, ((mtbv_refsp2 - mtbh_refsp2)).*ice_mask_totall2);
% pcolorpsn(lat_area, lon_area, ((mtbv_reffl2 - mtbh_reffl2)).*ice_mask_totall2);
% pcolorpsn(lat_area, lon_area, ((mtbv_refsp2 - mtbh_refsp2)./(mtbv_refsp2 + mtbh_refsp2)).*ice_mask_totall2);
pcolorpsn(lat_area, lon_area, ((nanmean(mtbvly,3) - nanmean(mtbhly,3))./(nanmean(mtbvly,3) + nanmean(mtbhly,3))).*ice_mask2);
% [C,h]= contourpsn(Lat,Lon,tmp_elev,levels,'k--');
plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
% set(gca, 'clim', [10 50]);
set(gca, 'clim', [0 0.15]);
% title(sprintf('%s','TBV - TBH'),'FontSize', fs);
title(sprintf('%s','(TBV - TBH)/(TBV + TBH)'),'FontSize', fs);
xlim([-0.7e6 0.9e6]);
ylim([-3.5e6 -0.5e6 ]);
zoom on;
% cmap=get_python_cmap('RdYlBu_r');
colormap(cmap(:,1:3))
c2=colorbar;
% c1.Location = 'southoutside';
c2.Location = 'south';
c2.FontSize = fs;
c1.FontSize = fs;
set(c2, 'Position', [0.73 0.0348 0.15 0.0182]);
set(gca,'FontWeight','bold')
title(c1,'TB [K]','FontSize', fs,'fontweight','bold');
title(c2,'NPR','FontSize', fs,'fontweight','bold');
axis off
% xlabel(t,'Easting [m]','FontWeight','bold','FontSize', fs);
% ylabel(t,'Norting [m]','FontWeight','bold','FontSize', fs);

filename = sprintf('%sFigure_1_mean_frozen_season_TB.png', out_dir);
% Save the figure with 300 dpi
exportgraphics(gcf, filename, 'Resolution', 300);