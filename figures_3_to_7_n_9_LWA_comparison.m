clear
%% Load ice edge
promice_dir = 'Z:/IceSheets/Weather Stations Data/Greenland/PROMICE/DATA/AWS/V9/combined/';
load([promice_dir,'/locations_promice_aws.mat']')
out_dir = 'Z:/IceSheets/Documentation/Manuscripts_n_abstracts/MWA_journal/Figs/';
driveletter   = 'Z:';
datadir1 = [driveletter, '/IceSheets/data/'];
promicedir = [datadir1, '/Promice_ice_mask/'];
A_GrIS = 1716555; %km2
load([promicedir, 'iceedge.mat']);
%% Fig 5: Elevation Map
% stations_sel = [3,4,12,19,32,33];
% figure
% set(gcf, 'position', [5.1126e+03 -44.2000 560 627.6000]);
% greenland('color', colo('grey1'));
% [tmp_elev, Lat, Lon] = load_gimpdem();
% % pcolorpsn(Lat,Lon,tmp_elev);
% % [C,h]= contourfpsn(Lat,Lon,tmp_elev,'k--');
% % [C,h]= contourfpsn(Lat,Lon,tmp_elev,[1000 2000 3000],"ShowText",true,"LabelFormat","%d m", ...
% %     "FaceAlpha",0.25);
% [C,h]= contourfpsn(Lat,Lon,tmp_elev,'k--');
% plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
% colormap parula(6)
% clim([500 3500])
% c=colorbar('Location', 'southoutside');
% c.Position = [0.3538 0.1018 0.3169 0.0340];
% ylabel(c,'Elevation from GIMP DEM [m]', 'FontWeight','bold')
% % plotpsn(lat, lon, 'b.','markersize',8);
% % plotpsn(latest_AWS_locations.lat, latest_AWS_locations.lon, 'r^','markersize',12);
% plotpsn(lat(stations_sel), lon(stations_sel), 'ro','markersize',24,'linewidth',4)
% % textpsn(lat(stations_sel), lon(stations_sel),promice_stations(stations_sel))
% axis off
% xlim([-0.7e6 0.9e6]);
% ylim([-3.5e6 -0.5e6 ]);
% zoom on;
% grid
% title('Selected PROMICE AWS','FontSize',12,'FontWeight','bold')
% saveas(gcf, sprintf('%s\\GrIS_DEM.png',out_dir));
%% Fig 3: SMAP Time Series
stn =2;
yrs =2016;
xl1 = datenum(yrs,3,1);
xl2 = datenum(yrs,11,1);
filename = sprintf('%s/LWA/SMAP_TB_time_series_%s.xls',out_dir,char(stations{stn}));
smap_data = readtable(filename);
smap_data = [datenum(smap_data.time),smap_data.TBV,smap_data.TBH,smap_data.TBV_ref,smap_data.TBV_std];
smap_dnum = smap_data(:,1);
[smap_data]=return_desired_data(smap_data,smap_dnum,yrs);
smap_dnum = smap_data(:,1);
TBV = smap_data(:,2);
TBH = smap_data(:,3);
TBV_ref = smap_data(:,4);
TBV_std = smap_data(:,5);
tbv_thresh = TBV_ref + min(12,10*TBV_std);

figure(Position=[5.9694e+03 -109 1.2328e+03 740.8000]);
t=tiledlayout(1,1,'TileSpacing','compact');

nexttile
plot(smap_dnum,TBV,'k.-','LineWidth',1.5,'MarkerSize',8)
hold
grid
ylim([140 280])
xlim([xl1,xl2]);
datetick('x','mm/dd/yyyy','keeplimits');
set(gca,'FontSize',14,'FontWeight', 'bold');

l1=legend('L-band TBV','L-band Threshold','FontSize',12,'FontWeight','bold');
legend('boxoff')
% l1.Position =  [0.1033    0.8554    0.0923    0.0619];

title(t,sprintf('%s',char(stations{stn})),'FontSize',14,'FontWeight', 'bold');
ylabel(t,'Brightness Temperature [K]','FontSize',14,'FontWeight', 'bold');
% Specify the output file path
filename = sprintf('%sL-band_TB_time_series_%s_%d.png',out_dir,char(stations{stn}),yrs);
% Save the figure with 300 dpi
% exportgraphics(gcf, filename, 'Resolution', 300);
%% Fig 5: 1 Year Comparison
stations=readtable(sprintf('%s/LWA/AWS.xlsx',out_dir), 'ReadVariableNames', false);
stations = table2cell(stations);
% stations_sel = [1,3,2]; yl = [0 80];
stations_sel = [4,5,6]; yl = [0 80];
figure('Position', [5.9902e+03 -261.4000 806.4000 994]);

t=tiledlayout(3,1,'TileSpacing','compact');
yrs=[2023];
y=1;
for stn = stations_sel
    % Load EBM (AWS) LWA time series
    filename1 = sprintf('%sLWA/Samimi_EBM_LWA_time_series_%s.xls',out_dir,char(stations{stn}));
    LWA_ebm = readtable(filename1);
    aws_days = datenum(LWA_ebm.time);
    aws_lwa_daily = LWA_ebm.aws_lwa_daily;  
    [aws_lwa_daily,aws_days]=return_desired_data(aws_lwa_daily,aws_days,yrs(y));

    % Load SMAP LWA time series
    filename = sprintf('%s/LWA/SMAP_LWA_time_series_%s.xls',out_dir,char(stations{stn}));
    LWA_smap = readtable(filename);
    smap_data = [datenum(LWA_smap.time),LWA_smap.LWA];
    [smap_data]=return_desired_data(smap_data,smap_data(:,1),yrs(y));
    % You may average twice daily samples for daily sample
    [smap_days,lwa_smap_daily]=compute_daily_mean(smap_data(:,1),smap_data(:,2));
    [smap_days,lwa_smap_daily] = matchup_time_series(aws_lwa_daily,aws_days,lwa_smap_daily,smap_days);

    if ~isempty(aws_lwa_daily)
        r1 = corr(aws_lwa_daily,lwa_smap_daily,'rows','complete');
        rmsd1 = rmse(aws_lwa_daily,lwa_smap_daily,"omitnan");
        nse1 = nash_sutcliffe_efficiency(aws_lwa_daily,lwa_smap_daily);

    end

    nexttile()
    plot(aws_days,aws_lwa_daily,'b-','LineWidth',3)
    hold   
    plot(smap_days,lwa_smap_daily,'r-.','LineWidth',3,'MarkerSize',16)

    xlim([datenum(yrs(y),6,1),datenum(yrs(y),11,1)]);
    ylim(yl)    
    grid
    % datetick('x','mmm/yyyy','keeplimits');
    datetick('x','dd mmm','keeplimits');
    title(sprintf(char(stations{stn})),'Interpreter', 'none')

    if ~isempty(aws_lwa_daily) && ~isnan(r1)
        ax = gca;
        text(ax,datenum(yrs(y),6,5), 75, sprintf('r = %0.2f', r1), ...%75
            'FontSize', 14, ...
            'FontWeight', 'bold', ...
            'Color', 'black');
        text(ax,datenum(yrs(y),6,5), 63, sprintf('RMSD = %0.0f mm', rmsd1), ...%63
            'FontSize', 14, ...
            'FontWeight', 'bold', ...
            'Color', 'black');
    end
    ax = gca;
    ax.XAxis.FontWeight = 'bold'; % Make X-axis tick labels bold
    ax.YAxis.FontWeight = 'bold'; % Make Y-axis tick labels bold
    set(gca,'FontSize',14,'FontWeight', 'bold');

end

l1=legend('EBM (AWS)','SMAP','FontSize',14,'FontWeight','bold');
legend('boxoff')

filename = sprintf('%sMWA_SE_revLUT2_corrected_data.png',out_dir);
% exportgraphics(gcf, filename, 'Resolution', 300);

%% Figure 6: 3 Way Comparison
stations=readtable(sprintf('%s/LWA/AWS.xlsx',out_dir), 'ReadVariableNames', false);
stations = table2cell(stations);
stations_sel = [3,1,2,6,5,4];
% figure('Position',   [5.0962e+03 -371 2.2872e+03 1.1076e+03]);
figure(Position=[6031 -515.8000 1.2772e+03 1.3544e+03])
t=tiledlayout(6,3,'TileSpacing','tight','Padding','compact');
yrs=[2021:2023];
c=1;
md1 = nan(18, 1);
md2 = nan(18, 1);
md3 = nan(18, 1);
sigma1 = nan(18, 1);
sigma2 = nan(18, 1);
sigma3 = nan(18, 1);
mad1 = nan(18, 1);
mad2 = nan(18, 1);
mad3 = nan(18, 1);
r1 = nan(18, 1);
r2 = nan(18, 1);
r3 = nan(18, 1);
rmsd1 = nan(18, 1);
rmsd2 = nan(18, 1);
rmsd3 = nan(18, 1);
nse1 = nan(18, 1);
nse2 = nan(18, 1);
nse3 = nan(18, 1);
clear mwa_smap_acc mwa_aws_acc mwa_gemb_acc
clear AWS yr smap_monset aws_monset gemb_monset smap_fonset aws_fonset gemb_fonset ...
    smap_mduration' aws_mduration' gemb_mduration' smap_maxsmelt' aws_maxsmelt' gemb_maxsmelt' ...
    mwa_smap_acc' mwa_aws_acc' mwa_gemb_acc
s=1;
for stn = stations_sel

    for y = 1:length(yrs)
        nexttile()
        % Load EBM (AWS) LWA time series
        filename1 = sprintf('%sLWA/Samimi_EBM_LWA_time_series_%s.xls',out_dir,char(stations{stn}));
        LWA_ebm = readtable(filename1);
        aws_days = datenum(LWA_ebm.time);
        aws_lwa_daily = LWA_ebm.aws_lwa_daily;
        [aws_lwa_daily,aws_days]=return_desired_data(aws_lwa_daily,aws_days,yrs(y));

        % Load SMAP LWA time series
        filename = sprintf('%s/LWA/SMAP_LWA_time_series_%s.xls',out_dir,char(stations{stn}));
        LWA_smap = readtable(filename);
        smap_data = [datenum(LWA_smap.time),LWA_smap.LWA];
        [smap_data]=return_desired_data(smap_data,smap_data(:,1),yrs(y));
        % You may average twice daily samples for daily sample
        [smap_days,lwa_smap_daily]=compute_daily_mean(smap_data(:,1),smap_data(:,2));

        % Read GEMB Time Series
        filename = sprintf('%s/LWA/GEMB_LWA_time_series_%s.xls',out_dir,char(stations{stn}));
        LWA_gemb = readtable(filename);
        gemb_data = [datenum(LWA_gemb.time),LWA_gemb.lwa_gemb];
        [gemb_data]=return_desired_data(gemb_data,gemb_data(:,1),yrs(y));

        if ~isempty(aws_lwa_daily) && nansum(aws_lwa_daily)>0
            [smap_days,lwa_smap_daily] = matchup_time_series(aws_lwa_daily,aws_days,lwa_smap_daily,smap_days);
            [gemb_days,gemb_lwa] = matchup_time_series(aws_lwa_daily,aws_days,gemb_data(:,2),gemb_data(:,1));
            md1(c) = nanmean((aws_lwa_daily-lwa_smap_daily));
            md2(c) = nanmean((gemb_lwa-lwa_smap_daily));
            md3(c) = nanmean((aws_lwa_daily-gemb_lwa));
            sigma1(c) = nanstd(aws_lwa_daily-lwa_smap_daily);
            sigma2(c) = nanstd(gemb_lwa-lwa_smap_daily);
            sigma3(c) = nanstd(aws_lwa_daily-gemb_lwa);
            mad1(c) = nanmean(abs(aws_lwa_daily-lwa_smap_daily));
            mad2(c) = nanmean(abs(gemb_lwa-lwa_smap_daily));
            mad3(c) = nanmean(abs(aws_lwa_daily-gemb_lwa));
            r1(c) = corr(aws_lwa_daily,lwa_smap_daily,'rows','complete');
            r2(c) = corr(gemb_lwa,lwa_smap_daily,'rows','complete');
            r3(c) = corr(aws_lwa_daily,gemb_lwa,'rows','complete');
            rmsd1(c) = rmse(aws_lwa_daily,lwa_smap_daily,"omitnan");
            rmsd2(c) = rmse(gemb_lwa,lwa_smap_daily,"omitnan");
            rmsd3(c) = rmse(aws_lwa_daily,gemb_lwa,"omitnan");
            nse1(c) = nash_sutcliffe_efficiency(aws_lwa_daily,lwa_smap_daily);
            nse2(c) = nash_sutcliffe_efficiency(gemb_lwa,lwa_smap_daily);
            nse3(c) = nash_sutcliffe_efficiency(aws_lwa_daily,gemb_lwa);          

        else
            [smap_days,lwa_smap_daily] = matchup_time_series(gemb_data(:,2),gemb_data(:,1),lwa_smap_daily,smap_days);
            [gemb_days,gemb_lwa] = matchup_time_series(gemb_data(:,2),gemb_data(:,1),gemb_data(:,2),gemb_data(:,1));
            md2(c) = nanmean((gemb_lwa - lwa_smap_daily));
            sigma2(c) = nanstd(gemb_lwa - lwa_smap_daily);
            mad2(c) = nanmean(abs(gemb_lwa - lwa_smap_daily));
            r2(c) = corr(gemb_lwa, lwa_smap_daily, 'rows', 'complete');
            rmsd2(c) = rmse(gemb_lwa, lwa_smap_daily, "omitnan");
            nse2(c) = nash_sutcliffe_efficiency(gemb_lwa, lwa_smap_daily);
        end


        [smap_monset(c),smap_fonset(c),smap_mduration(c),smap_maxsmelt(c)] = melt_metrics(smap_days,lwa_smap_daily,yrs(y));
        [aws_monset(c),aws_fonset(c),aws_mduration(c),aws_maxsmelt(c)] = melt_metrics(aws_days,aws_lwa_daily,yrs(y));
        [gemb_monset(c),gemb_fonset(c),gemb_mduration(c),gemb_maxsmelt(c)] = melt_metrics(gemb_days,gemb_lwa,yrs(y));
        AWS(c,:) = string(stn);
        yr(c) = yrs(y);


        plot(gemb_days,gemb_lwa,':','color',[0,0,0],'LineWidth',3,'MarkerSize',6)
        hold
        plot(aws_days,aws_lwa_daily,'b-','LineWidth',3)             
        plot(smap_days,lwa_smap_daily,'r-.','LineWidth',3,'MarkerSize',16)    
   
        text(datenum(yrs(y),6,5), 80, sprintf('(%s)', char(96+c)), ...
            'FontSize', 14, ...
            'FontWeight', 'bold', ...
            'Color', 'black');

        xlim([datenum(yrs(y),6,1),datenum(yrs(y),11,1)]);
        grid on;grid minor
        datetick('x','mm/dd','keeplimits');
        % datetick('x','dd mmm','keeplimits');
        if c<4
            title(sprintf('%d',yrs(y)))
        end

        if mod(c,3) == 0
            ax = gca;
            text(ax,1.05, 0.6, sprintf('%s', char(stations{stn})), ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 14, ...
                'FontWeight', 'bold', ...
                'Color', 'black', ...
                'Interpreter', 'none',...
                'Rotation', 270);
        end

        ylim([0 90])
        ax = gca;
        ax.XAxis.FontWeight = 'bold'; % Make X-axis tick labels bold
        ax.YAxis.FontWeight = 'bold'; % Make Y-axis tick labels bold
        set(gca,'FontSize',14,'FontWeight', 'bold');

        mwa_smap_acc(c) = nansum(lwa_smap_daily(:));
        mwa_aws_acc(c) = nansum(aws_lwa_daily(:));
        mwa_gemb_acc(c) = nansum(gemb_lwa(:));

        c=c+1;

    end
    s = s+1;
end

T1=table(AWS, yr', md1, md2, md3, sigma1, sigma2, sigma3, mad1, mad2, mad3, r1, r2, r3, rmsd1, rmsd2, rmsd3, nse1, nse2, nse3);
% writetable(T1,sprintf('%sPerformace_matrics_promice_t_surf_verify.xls',out_dir),'Sheet',1);

T2=table(AWS,yr',smap_monset',aws_monset',gemb_monset',smap_fonset',aws_fonset',gemb_fonset',...
    smap_mduration',aws_mduration',gemb_mduration',smap_maxsmelt',aws_maxsmelt',gemb_maxsmelt',...
    mwa_smap_acc',mwa_aws_acc',mwa_gemb_acc');
% writetable(T2,sprintf('%sPerformace_matrics_promice_t_surf_verify.xls',out_dir),'Sheet',2);

ylabel(t,'Liquid Water Amount [mm]','FontSize',14,'FontWeight','bold')
xlabel(t,'Month','FontSize',14,'FontWeight','bold')

l1=legend('GEMB (ERA5)','EBM (AWS)','SMAP','FontSize',12,'FontWeight','bold');
legend('boxoff')
l1.Position =  [   0.8184    0.8944    0.0923    0.0619 ];

% Specify the output file path
filename = sprintf('%sLWA_3_way_comparision_promice_corrected_data_revA.png', out_dir);
% Save the figure
% exportgraphics(gcf, filename, 'Resolution', 300);

%% Figure 7: LWA Time series
stations=readtable(sprintf('%s/LWA/AWS.xlsx',out_dir), 'ReadVariableNames', false);
stations = table2cell(stations);
stations_sel = [3,6,1,5,2,4];
figure('Position',  [5.6434e+03 -314.2000 1740 1.0508e+03]);
t=tiledlayout(3,2,'TileSpacing','compact');
yrs=[2021:2023];
fig = [1 4 2 5 3 6];
c=1;
for stn = stations_sel
    % Load SMAP LWA time series
    filename = sprintf('%s/LWA/SMAP_LWA_time_series_%s.xls',out_dir,char(stations{stn}));
    LWA_smap = readtable(filename);
    % You may average twice daily samples for daily sample
    [smap_days,lwa_smap_daily]=compute_daily_mean(datenum(LWA_smap.time),LWA_smap.LWA);
    nexttile()
    b= bar(smap_days,lwa_smap_daily,'FaceColor',[0.8500 0.3250 0.0980]);
    b(1).BarWidth = 15;
    ylim([0 70])
    xlim([datenum(2015,1,1),datenum(2024,1,1)]);

    ax = gca;
    ax.YColor = [0.8500 0.3250 0.0980];

    grid
    datetick('x','mm/yy','keeplimits');
    title(sprintf('%s',char(stations{stn})),'FontSize',14,'FontWeight','bold',Interpreter='none')
    ax = gca;
    ax.YColor = 'k';    
    set(gca, 'FontSize',14,'FontWeight', 'bold');

    text(datenum(2015,6,1), 65, sprintf('(%s)', char(96+fig(c))), ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'Color', 'black');

    c = c+1;
end
% end

ylabel(t,'Liquid Water Amount [mm]','FontSize',14,'FontWeight','bold','Color',[0.0 0.0 0.0]);
xlabel(t,'Year','FontSize',14,'FontWeight','bold')

filename = sprintf('%sTB&MWA_timeSeries_LUT2_revised.png',out_dir);
% exportgraphics(gcf, filename, 'Resolution', 300);

%% Figure 9: Modeled vs Measured Subsurface Temperature
stations=readtable(sprintf('%s/LWA/AWS.xlsx',out_dir), 'ReadVariableNames', false);
stations = table2cell(stations);
stn = 1;
yr=[2021];
fs=14;
filename = sprintf('%s/LWA/Modeled_and_Measured_Temperature_%s.nc',out_dir,char(stations{stn}));
aws_dnum = ncread(filename, 'date_number');
depth_meas= ncread(filename, 'depth_meas');
temp_meas= ncread(filename, 'temp_meas');
temp_mod= ncread(filename, 'temp_mod');
depth_mod = ncread(filename, 'depth_mod');
dnum_meas = repmat(aws_dnum,[1 12]);

depth_mod = repmat(depth_mod',[size(temp_mod,1) 1]);
dnum_mod = repmat(datenum(aws_dnum),[1 size(temp_mod,2)]);
figure(Position=[6.2022e+03 -239.8000 776 852.4000])
t=tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
% plot_profile(aws_dnum,depth_mod,temp_mod')
scatter(dnum_mod(:),depth_mod(:),20,temp_mod(:),'filled')
set(gca,'YDir','reverse')
title('Modeled Temperature Profile','FontWeight','bold', 'FontSize',fs)
% contourf(aws_dnum,z,aws_data1(:,6:48)',40,"ShowText",true,"LabelFormat","%0.1f C");
shading flat;
colorbar("off")
xlim([datenum(yr,6,17),datenum(yr,11,1)]);
ylim([0 10])
clim([-40 0])
datetick('x','mm/dd/yyyy','keeplimits');
yticks([0:10])
yticklabels([0:10])
cmap=get_python_cmap('plasma');
cmap3=cmap(:,1:3);
cmap3(end,:)=[0 1 0];
colormap(ax1,cmap3)
yline(1:9,'LineWidth',0.25,'Color',[0.5,0.5,0.5],LineStyle='--')
xline(xticks,'LineWidth',0.25,'Color',[0.5,0.5,0.5],LineStyle='--')
set(gca,'TickDir','both','TickLength',[0.0052,0.0052],'XMinorTick','off')
set(gca,'FontWeight','bold')
set(gca, 'FontSize',fs)

ax2 = nexttile;
% scatter(temp_aws{stn}(:,4),-1*temp_aws{stn}(:,9),20,temp_aws{stn}(:,8),'filled')
scatter(dnum_meas(:),-1*depth_meas(:),20,temp_meas(:),'filled')
xlim([datenum(yr,6,17),datenum(yr,11,1)]);
ylim([-10 0])
clim([-40 0])
datetick('x','mm/dd/yyyy','keeplimits');
yticks([-10:1:0])
yticklabels([[-10:1:-1]*-1,0])
cmap=get_python_cmap('plasma');
cmap3=cmap(:,1:3);
cmap3(end,:)=[0 1 0];
colormap(ax2,cmap3)
title('Measured Temperature Profile','FontWeight','bold', 'FontSize',fs)
yline(-1*[1:9],'LineWidth',0.25,'Color',[0.5,0.5,0.5],LineStyle='--')
xline(xticks,'LineWidth',0.25,'Color',[0.5,0.5,0.5],LineStyle='--')
set(gca,'TickDir','both','TickLength',[0.0052,0.0052],'XMinorTick','off')
set(gca,'FontWeight','bold')
set(gca, 'FontSize',fs)

cbh = colorbar;
cbh.Layout.Tile = 'east';
cbh.FontSize=fs;

filename = sprintf('%sfig_9_meas_n_mod_temp.png', out_dir);
% Save the figure with 300 dpi
% exportgraphics(gcf, filename, 'Resolution', 300);

