%%% TB variation along a lat transect
% Define the transect of your interest: trans_lat,minlon,maxlon
% Contact: Alamgir Hossan JPL, 1/27/2024 alamgir.hossan@jpl.nasa.gov
%%%

clear
% Direct to the data 
input_drive = 'Z:';
input_dir = [input_drive, '/IceSheets/data/Greenland Data/rSIR/SMAP/NP/essential_vars/v2/'];
input_dir2 = [input_drive, '/IceSheets/data/Greenland Data/rSIR/SSMIS/F18/NP/'];
ouput_dir = 'Z:\IceSheets\EM_Model_SnoWR_Algorithm\SnoWR_Results\Figs\Greenland\time_series_analysis\';

% load lat/lon of rSIR 3.126km data subset 
[lat_area,lon_area] = return_rSIR_lat_lon_subset_v1(input_drive);

promice_dir = 'Z:\IceSheets\Weather Stations Data\Greenland\PROMICE\DATA\AWS\V9\combined\';
load([promice_dir,'/locations_promice_aws.mat']')


%%
yrs = 2015:2023;

for yr=5%:length(yrs)

    filename = ['NSIDC-0738-EASE2_N3.125km-SMAP_LRM-',num2str(yrs(yr)),'.mat'];
    disp([input_dir,filename])
    smap = load([input_dir,filename]);

    % Tbv_smap = cat(3,smap.Tbvm,smap.Tbve);
    % dnum_smap = cat(3,smap.dnumm,smap.dnume);

    % [dnum,ind] = sort(dnum);
    % Tbv = tbv(ind);

    filename2 = ['NSIDC-0630-EASE2_N3.125km-F18_SSMIS-',num2str(yrs(yr)),'.mat'];
    ssmis = load([input_dir2,filename2]);

    % end

%% Provide the lat of the transect you'd like
trans_lat = 67.07; % K transect
% trans_lat = 69.87162; % CP1 transect
minlon = -50; % define min and max lon
maxlon = -36;
step = 0.1;
[r,c,ind]=retun_indexes_lat_transect_gris(trans_lat,minlon,maxlon,step);

    %% TBs
    tbl=smap.Tbve;
    dnuml=smap.dnume;
    tblm=nanmean(tbl(:,:,91:96),3);
    tbls=nanstd(tbl(:,:,91:96),0,3);
    tbl_tr = tblm(ind);
    lon_tr = lon_area(ind);
    Z=4;
    thresh=Z*tbls(ind);


    tbka=ssmis.Tbve;
    dnumka=ssmis.dnume;
    tbkam=nanmean(tbka(:,:,91:96),3);
    tbka_tr=tbkam(ind);

    %
    % doy = doys(1);
    tbl_anom = tbl -tblm;
    tbka_anom = tbka - tbkam;
    %% Select dates, or years below
    dnums = datenum(2019, 7,  31);
    doys=dnums - datenum(yrs(yr),1,1)+1;

    % videoFile = [ouput_dir,'TB_along_CP1_transect_',num2str(yrs(yr)),'.mp4'];
    % writerObj = VideoWriter(videoFile, 'MPEG-4');
    % writerObj.FrameRate = 5; % Adjust the frame rate as needed
    % 
    % open(writerObj);

    for doy=doys(1)%100:250 % Define DOY of interest

        tbls=tbl(:,:,doy);
        tbls_tr = tbls(ind);
        tbls_tr_anom = tbls_tr - tbl_tr;
        meltl = find(abs(tbls_tr_anom)>=thresh);
        frozen = find(abs(tbls_tr_anom)<thresh);

        tbkas=tbka(:,:,doy);
        tbkas_tr = tbkas(ind);
        ka_thresh = 245;

        ind_aws_nearby = find(abs(lat - trans_lat)<0.25);
        lon_aws_nearby=lon(ind_aws_nearby);
        lat_aws_nearby=lat(ind_aws_nearby);
        aws_nearby = string(stations(ind_aws_nearby));

        figure(Position=[4.5342e+03 -228.2000 916.4000 663.2000])
        plot(lon_tr,tbl_tr,'k-','LineWidth',2);
        grid
        hold
        % plot(lon_tr,tbl_tr+thresh,'b--',LineWidth=1)
        % plot(lon_tr,tbl_tr-thresh,'b--',LineWidth=1)

        % x_fill = [lon_tr; flipud(lon_tr)];
        % y_fill_upper = [tbl_tr + thresh; flipud(tbl_tr - thresh)];
        % fill(x_fill, y_fill_upper, 'b', 'FaceAlpha', 0.1);

        title(sprintf('%s',datestr(datenum(yrs(yr),1,1)+doy-1)))
        ylabel('Brightness Temperature (TB) [K]',FontWeight='bold')
        xlabel('Longitude [deg.]',FontWeight='bold')        

        plot(lon_tr,tbls_tr,'k-.','LineWidth',1)
        % plot(lon_tr(frozen),tbls_tr(frozen),'b.','LineWidth',1,MarkerSize=6)
        plot(lon_tr(meltl),tbls_tr(meltl),'ks','LineWidth',1,MarkerSize=4)

        plot(lon_tr,tbka_tr,'b:','LineWidth',2);
        % yline(ka_thresh,':','Color',[0.8,0,0],'LineWidth',1.5)
        plot(lon_tr,tbkas_tr,'b--','LineWidth',1)
        % plot(lon_tr(tbkas_tr<ka_thresh),tbkas_tr(tbkas_tr<ka_thresh),'r.','LineWidth',1,MarkerSize=6)
        plot(lon_tr(tbkas_tr>=ka_thresh),tbkas_tr(tbkas_tr>=ka_thresh),'bo','LineWidth',1,MarkerSize=3)

        text(lon_aws_nearby,140*(ones(size(lon_aws_nearby))),aws_nearby,FontWeight="bold")
        % plot(lon_aws_nearby,145*(ones(size(lon_aws_nearby))),'k*')
        ylim([130 280])

        l=legend('L-band TBV Ref.','L-band TBV','L-band Melt Flag','Ka--band TBV Ref.','Ka-band TBV','Ka-band Melt Flag');
        % l=legend('L-band TBV ref','L-band thresh','L-band TBV','L-band Melt Flag','Ka-band Threshold','Ka-band TBV','Ka-band Melt Flag')
        l.Box = "off";
        l.NumColumns = 2;

        set(gca,"FontWeight","bold")


        % currentFrame = getframe(gcf);
        % 
        % writeVideo(writerObj, currentFrame);
        % close
    end

    % 
    % close(writerObj);
    % 
    % disp(['Video created: ' videoFile]);

end

