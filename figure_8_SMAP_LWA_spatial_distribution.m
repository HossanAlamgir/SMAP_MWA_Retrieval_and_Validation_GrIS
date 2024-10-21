clear
out_dir='Z:\IceSheets\Documentation\Manuscripts_n_abstracts\MWA_journal\Figs\';
%% Load ice edge
driveletter   = 'Z';
datadir1 = [driveletter, ':\IceSheets\data\'];
promicedir = [datadir1, '\Promice_ice_mask\'];
% promicedir = [HOMEDIR, '_Temp\icesheetdata\'];
A_GrIS = 1716555; %km2 -- https://www.promice.dk/RetreatingIce.html;jsessionid=0AB4246BAECD765D889EE6659AB92686?cid=3816
load([promicedir, 'iceedge.mat']);
load("ice_mask_constructed_from_L1_9km_ice_mask.mat")
input_dir = [driveletter, ':/IceSheets/Melt_Detection_Results-AH/Greenland/rSIR/3.125km/v4/Z_eq_10_comp_ref_revB/evening/'];
%% SMAP LWA Maps
fs=16;
input_drive = 'Z:';
pass='Evening';
fpath=[input_drive,'/IceSheets/EM_Model_SnoWR_Algorithm/SnoWR_Results/Greenland-AH/v3/transformed/v1/'];
[lat_area,lon_area] = return_rSIR_lat_lon_subset_v1('Z:'); % lat/lon of corresponding rSIR subset
%
figure('position',  [5.3862e+03 -239.8000 2.1928e+03 1.0856e+03]);
t=tiledlayout(2,5, 'Padding', 'compact', 'TileSpacing', 'compact');
yrs=2015:2023;
% yrs=2018;

for y=1:length(yrs)
    sname=['SnoWR_output_','V_',pass,'_',num2str(yrs(y)),'.mat'];
    load([fpath, sname]);
    mwa=nansum(snowr.mwa,3);
    %%%
    filename = ['SMAP_melt_flags_',num2str(yrs(y)),'.mat'];
    disp([input_dir,filename])
    load([input_dir,filename])

    temp_incr = nan([size(lat_area),366]);
    temp_decr = nan([size(lat_area),366]);

    temp_incr(tbv_anom_refw <= - thresh_tbv )=1;
    temp_incr2=double(nansum(temp_incr,3)  > 0 );
    % I_inc0=find(temp_incr2);

    temp_decr(tbv_anom_refw > thresh_tbv )=1;
    temp_decr2=double(nansum(temp_decr,3)  > 0 );
    % I_dec0=find(temp_decr2);

    temp_incr2(temp_incr2==0)=nan;
    temp_decr2(temp_decr2==0)=nan;

    % temp_incr2(temp_decr2>0)=nan; %don't/ mask the inc/dec pixels
    temp_incr2(temp_incr2>0)=-5;

    lwa_inc_TB = mwa.*temp_decr2.*ice_mask2;
    lwa_dec_TB = mwa.*temp_incr2.*ice_mask2 - 1; % For the mask

    %%%
    % [tmp_elev, Lat, Lon] = load_gimpdem();
    % pcolorpsn(Lat,Lon,tmp_elev);
    % [C,h]= contourfpsn(Lat,Lon,tmp_elev,'k--');
    % [C,h]= contourfpsn(Lat,Lon,tmp_elev,[1000 2000 3000],"ShowText",false,"LabelFormat","%d m", ...
    %     "FaceAlpha",0.25);
    % [C,h]= contourfpsn(Lat,Lon,tmp_elev,[1000 2000 3000],'k--');
    %%%
    h(y)=nexttile;
    % set(gcf, 'position', [64.2000   13.8000  560.0000  849.6000]);
    greenland('color', colo('grey1'));
    pcolorpsn(lat_area, lon_area,lwa_inc_TB );
    pcolorpsn(lat_area, lon_area,lwa_dec_TB);
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
    % set(gca, 'clim', [-1 3000]);
    clim([-5 3000]);
    title(sprintf('%s',num2str(yrs(y))),'FontSize', fs);
    xlim([-0.7e6 0.9e6]);
    ylim([-3.5e6 -0.5e6 ]);

    % ylim([-3.35e6 -2e6 ]);
    % xlim([-0.4e6 0.7e6]);

    lwa = mwa.*temp_decr2.*ice_mask2;
    avg_lwa(y) = nanmean(lwa(:));
    zoom on;
    % colormap(cmap1);
    % cmap=get_python_cmap('Spectral');
    % cmap=get_python_cmap('BrBG');
    % cmap=get_python_cmap('RdYlBu');
    cmap=get_python_cmap('plasma_r', 600);
    % cmap(1,:) = 1;
    % cmap(1,:) = [0.1 .75 0.1 1];
    cmap = [0.7 .7 0.7 0;cmap];
    colormap(cmap(:,1:3))

    % colorbar;
    % pause(1)
    % file_suffix=datestr(dnum_pm_ave(idates(u)), 'mm-dd-yyyy');
    % saveas(gcf, sprintf('%s\\MWA_%s.png',...
    %     [out_dir,'\MWA\'],file_suffix));
    % close

end
xlabel(t,'Easting [m]','FontWeight','bold','FontSize', fs);
ylabel(t,'Norting [m]','FontWeight','bold','FontSize', fs);
cbh = colorbar(h(end));
cbh.Layout.Tile = 'east';
cbh.FontSize=fs;
cbh.Ticks = linspace(0, 3000, 7);
% cb_title=title(cbh, 'Liquid Water Amounts [mm]','Rotation',90);


filename = sprintf('%sAnnual_MWA_LUT2_1.png',out_dir);
exportgraphics(gcf, filename, 'Resolution', 300);

