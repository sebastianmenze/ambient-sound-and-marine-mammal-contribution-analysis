%%  ambient sound analysis code 

% by Sebastian Menze, 2016
% Institute of Marine Research, Bergen, Norway
% Alfred Wegener Instituite, Bremerhaven, Germany

% Supplementary material to the paper:
% "The influence of sea ice, wind speed and marine mammals ...
% on Southern Ocean ambient sound"
% published in Royal Society Open Science

% Developed and tested using Matlab 2015a with all toolboxes

% before execution download additional m-file and ERA INERIM dataset
% manually, see instructions below:

% shadedErrorBar from Rob Campbell - November 2009
% available on Matlab file exchange 
% Link: http://www.mathworks.com/matlabcentral/fileexchange/26311-shadederrorbar

% get ERA-INTERIM reanalysis wind speed data (u10, v10) from ECMWF website
% http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/
% personal account (free) necessary for downloading data
% insert path to downloaded nc file here and unncomment:
% nc_file='netcdf-atls05-20141001165205-429-0777.nc';

% now the script can be executed and all additional data will be 
% downloaded automatically (requires internet access)
% analysis and plotting will run automatically, reproducing the figures 

% NOTICE: Script takes several hours to complete

%% load ambient sound dataset from netcdf file

nc_filename='AWI_230_6_ambient_sound_data.nc';

info=ncinfo(nc_filename);
info.Variables.Name

a66.data.spec_dB_5min=double(ncread(nc_filename,'psd_in_db_re_1_squared_upa_over_hz'))';
a66.data.frec=double(ncread(nc_filename,'frequency'));
time_strings=ncread(nc_filename,'time');
a66.data.startdate=datenum(time_strings,'yyyy-mm-ddTHH:MM:SS');
a66.data.spl_broadband_5min=double(ncread(nc_filename,'spl_in_db_re_1_upa'))';

clear time_strings info

nc_filename='AWI_232_9_ambient_sound_data.nc';

info=ncinfo(nc_filename);
info.Variables.Name

a69.data.spec_dB_5min=double(ncread(nc_filename,'psd_in_db_re_1_squared_upa_over_hz'))';
a69.data.frec=double(ncread(nc_filename,'frequency'));
time_strings=ncread(nc_filename,'time');
a69.data.startdate=datenum(time_strings,'yyyy-mm-ddTHH:MM:SS');
a69.data.spl_broadband_5min=double(ncread(nc_filename,'spl_in_db_re_1_upa'))';

clear time_strings info

% %% plot longterm spectra with full resolution (requires much RAM!)
% 
% f1=figure(1);
% set(f1,'color',[1 1 1]);
% 
% h = surf(data.startdate,data.frec,data.spec_dB_5min,'EdgeColor','none');
% axis xy; axis tight; colormap(jet); view(0,90)
% datetick('x','yy-mm','keeplimits')
% ylabel('Frequency (Hz)');
% set (gca,'yscale','log')
% grid off
% colormap jet
% set(gca,'clim',[40 100])
% ylim([10 data.samplerate/2])
% t=title([data.ID ' longterm spectrogram'])
% set(t,'interpreter','none')
% clear t
% set(gca,'Tickdir','out')
% colorbar


%% plot long term spectrogram with limited RAM 

frec_bins=logspace(log10(10),log10(16000),1000);
frec_center= frec_bins(1:end-1) + diff(frec_bins)/2 ;

for i_bin=2:numel(frec_bins)
    ix_frec=a66.data.frec>frec_bins(i_bin)-1 & a66.data.frec<frec_bins(i_bin);
    a66_interp(:,i_bin-1)=mean(a66.data.spec_dB_5min(ix_frec,:),1);
    a69_interp(:,i_bin-1)=mean(a69.data.spec_dB_5min(ix_frec,:),1);
end 


 f1=figure(1);
clf
set(f1,'color',[1 1 1]);

subplot(211)
h = surf(a66.data.startdate,frec_center,a66_interp','edgecolor','none');
axis tight; colormap(jet); view(0,90)
box on
    ylabel('Frequency in Hz');
 set (gca,'yscale','log')
    grid off
    colormap jet
    set(gca,'clim',[40 100])
    ylim([10 16000])  
    xlim([datenum(2008,3,10) datenum(2010,12,17)])
set(gca,'Tickdir','out')
set(gca,'XTick',[datenum(2008,1,1),datenum(2008,7,1),datenum(2009,1,1),datenum(2009,7,1),datenum(2010,1,1),datenum(2010,7,1)])

datetick('x','yyyy-mm','keeplimits', 'keepticks')
cb=colorbar
ylabel(cb,'PSD in dB re 1 \mu Pa^2 Hz^-^1')
text(datenum(2008,04,01),25000,200,'a) 66°S','fontweight','bold')

subplot(212)
h = surf(a69.data.startdate,frec_center,a69_interp','edgecolor','none');
axis tight; colormap(jet); view(0,90)
box on
    ylabel('Frequency in Hz');
 set (gca,'yscale','log')
    grid off
    colormap jet
    set(gca,'clim',[40 100])
    ylim([10 16000])  
    xlim([datenum(2008,3,10) datenum(2010,12,17)])
    
set(gca,'Tickdir','out')
set(gca,'XTick',[datenum(2008,1,1),datenum(2008,7,1),datenum(2009,1,1),datenum(2009,7,1),datenum(2010,1,1),datenum(2010,7,1)])

datetick('x','yyyy-mm','keeplimits', 'keepticks')

cb=colorbar
ylabel(cb,'PSD in dB re 1 \mu Pa^2 Hz^-^1')
text(datenum(2008,04,01),25000,200,'b) 69°S','fontweight','bold')
drawnow

%  set(gcf,'PaperPositionMode','auto')
%  print(gcf,'-dpdf','longtermspecs_5min','-r400') 
%   print(gcf,'-dtiff','longtermspecs_5min','-r400') 

%% plot empirical probablity spectra


clear histval i nozero_histval 
for i=1:size(a66.data.spec_dB_5min,1)
a66.histval(i,:)=hist(a66.data.spec_dB_5min(i,:),40:0.1:110);
% disp(i/size(a66.data.spec_dB_5min,1))
end

ix=a66.histval==0;
a66.nozero_histval=a66.histval;
a66.nozero_histval(ix)=NaN;

clear histval i nozero_histval 
for i=1:size(a69.data.spec_dB_5min,1)
a69.histval(i,:)=hist(a69.data.spec_dB_5min(i,:),40:0.1:110);
% disp(i/size(a69.data.spec_dB_5min,1))
end

ix=a69.histval==0;
a69.nozero_histval=a69.histval;
a69.nozero_histval(ix)=NaN;

per66=prctile(a66.data.spec_dB_5min,[5 50 95],2);
per69=prctile(a69.data.spec_dB_5min,[5 50 95],2);

figure(2)
clf
set(gcf,'color',[1 1 1])

subplot(211)
hold on

h = surf(a66.data.frec,40:0.1:110,(a66.nozero_histval./size(a66.data.spec_dB_5min,2) )','EdgeColor','none')
    axis xy; axis tight; colormap(jet); view(0,90)
    xlabel('Frequency (Hz)');
    colormap jet
 set(gca,'xscale','log','xlim',[10 16000],'ylim',[40 110],'clim',[0 0.01],'zlim',[0 500])
grid on
cb=colorbar
ylabel(cb,'Empirical probability density')

plot3(a66.data.frec,filtfilt(ones(3,1)/3,1,per66(:,1)),ones(numel(a66.data.frec),1),'-k','linewidth',1.5)
plot3(a66.data.frec,filtfilt(ones(3,1)/3,1,per66(:,2)),ones(numel(a66.data.frec),1),'-k','linewidth',1.5)
plot3(a66.data.frec,filtfilt(ones(3,1)/3,1,per66(:,3)),ones(numel(a66.data.frec),1),'-k','linewidth',1.5)


% plot noise floor
ylabel('PSD in db re 1 \muPa^2 Hz^{-1}')
plot3([10,100,1000,20000],[54,42,42,42],200*ones(4,1),':k','linewidth',2)

text(1900,100,'a) 66°S','fontweight','bold')

subplot(212)
hold on
h = surf(a69.data.frec,40:0.1:110,(a69.nozero_histval./size(a69.data.spec_dB_5min,2) )','EdgeColor','none')
    axis xy; axis tight; colormap(jet); view(0,90)
    xlabel('Frequency (Hz)');
    colormap jet
 set(gca,'xscale','log','xlim',[10 16000],'ylim',[40 110],'clim',[0 0.01],'zlim',[0 500])
grid on
cb=colorbar
ylabel(cb,'Empirical probability density')

plot3(a69.data.frec,filtfilt(ones(3,1)/3,1,per69(:,1)),ones(numel(a66.data.frec),1),'-k','linewidth',1.5)
plot3(a69.data.frec,filtfilt(ones(3,1)/3,1,per69(:,2)),ones(numel(a66.data.frec),1),'-k','linewidth',1.5)
plot3(a69.data.frec,filtfilt(ones(3,1)/3,1,per69(:,3)),ones(numel(a66.data.frec),1),'-k','linewidth',1.5)

% plot noise floor
ylabel('PSD in db re 1 \muPa^2 Hz^{-1}')
plot3([10,100,1000,20000],[54,42,42,42],200*ones(4,1),':k','linewidth',2)

text(1900,100,'b) 69°S','fontweight','bold')

% 
%  set(gcf,'PaperPositionMode','auto')
%  print(gcf,'-dtiff','empirical_spectral_prob_density_both_recorders','-r400') 
%   print(gcf,'-dpdf','empirical_spectral_prob_density_both_recorders','-r400') 

clearvars -except a66 a69 

%% calculate marine mammal contribution SPL and SNR

for i_time=1:size(a66.data.spec_dB_5min,2)
        
    f=a66.data.frec;
    db_recorded=a66.data.spec_dB_5min(:,i_time);
    pxx=power(10,db_recorded/10);
    
    %%%% blue and fin whale contribution fit
    p=pxx;
    bio_spec= f<10 | f>15 & f<30  | f>50;
    p(bio_spec)=NaN;
    db=10*log10(p);
    
    f_points=logspace(1,4,100);
    db_points=interp1(f,db,f_points);
    [xData, yData] = prepareCurveData( f_points, db_points );
    ft = fittype( 'power2' );
    
    [fitresult, gof] = fit( xData, yData, ft);
    db_fit = feval(fitresult,f);
    p_fit=10.^(db_fit/10);
    
    %%%% blue whale contribution SPL
    con_abs(i_time)=   (sqrt( sum(  pxx(f>=26 & f<=28) ) )) ;
    con_noise(i_time)=   (sqrt( sum(  p_fit(f>=26 & f<=28) ) )) ;
    
    blue_con_snr(i_time)= 20*log10(con_abs(i_time)) - 20*log10(con_noise(i_time));
    if  blue_con_snr(i_time)>3
        blue_whale_con_p(i_time)= con_abs(i_time) - con_noise(i_time) ;
        blue_whale_con_db(i_time)=real(20*log10(blue_whale_con_p(i_time)));
    else
        blue_whale_con_p(i_time)=NaN;
        blue_whale_con_db(i_time)=NaN;
    end
    
    %%%% fin whale contribution  
    p=pxx;
    bio_spec= f<50 | f>97 & f<101  | f>150;
    p(bio_spec)=NaN;
    db=10*log10(p);
    
    f_points=logspace(1,4,100);
    db_points=interp1(f,db,f_points);
    [xData, yData] = prepareCurveData( f_points, db_points );
    ft = fittype( 'poly4' );
    
    [fitresult, gof] = fit( xData, yData, ft);
    db_fit = feval(fitresult,f);
    p_fit=10.^(db_fit/10);
    
    con_abs=   (sqrt( sum(  pxx(f>=96 & f<=99) ) )) ;
    con_noise=   (sqrt( sum(  p_fit(f>=96 & f<=99) ) )) ;
    
    fin_con_snr(i_time)= 20*log10(con_abs) - 20*log10(con_noise);
    if  fin_con_snr(i_time)>1
        fin_whale_con_p(i_time)= con_abs - con_noise ;
        fin_whale_con_db(i_time)=real(20*log10(fin_whale_con_p(i_time)));
    else
        fin_whale_con_p(i_time)=NaN;
        fin_whale_con_db(i_time)=NaN;
    end
    
    %%%% minke whale contribution  
    p=pxx;
    bio_spec= f<30 | f>97 & f<500  | f>1000;
    p(bio_spec)=NaN;
    db=10*log10(p);
    
    f_points=logspace(1,4,100);
    db_points=interp1(f,db,f_points);
    [xData, yData] = prepareCurveData( f_points, db_points );
    ft = fittype( 'power2' );
    
    [fitresult, gof] = fit( xData, yData, ft);
    db_fit = feval(fitresult,f);
    p_fit=10.^(db_fit/10);
    
    con_abs=   (sqrt( sum(  pxx(f>=105 & f<=300) ) )) ;
    con_noise=   (sqrt( sum(  p_fit(f>=105 & f<=300) ) )) ;
    
    minke_con_snr(i_time)= 20*log10(con_abs) - 20*log10(con_noise);
    if  minke_con_snr(i_time)>3
        minke_whale_con_p(i_time)= con_abs - con_noise ;
        minke_whale_con_db(i_time)=real(20*log10(minke_whale_con_p(i_time)));
    else
        minke_whale_con_p(i_time)=NaN;
        minke_whale_con_db(i_time)=NaN;
    end
    
    %%%% leopard seal contribution   
    con_abs=   (sqrt( sum(  pxx(f>=320 & f<=350) ) )) ;
    con_noise=   (sqrt( sum(  p_fit(f>=320 & f<=350) ) )) ;
    
    leo_con_snr(i_time)= 20*log10(con_abs) - 20*log10(con_noise);
    if  leo_con_snr(i_time)>3
        leo_whale_con_p(i_time)= con_abs - con_noise ;
        leo_whale_con_db(i_time)=real(20*log10(leo_whale_con_p(i_time)));
    else
        leo_whale_con_p(i_time)=NaN;
        leo_whale_con_db(i_time)=NaN;
    end
    
end

a66.time=a66.data.startdate;
a66.blue_whale_con_db=blue_whale_con_db;
a66.blue_whale_con_snr=blue_con_snr;
a66.fin_whale_con_db=fin_whale_con_db;
a66.fin_whale_con_snr=fin_con_snr;
a66.minke_whale_con_db=minke_whale_con_db;
a66.minke_whale_con_snr=minke_con_snr;
a66.leo_whale_con_db=leo_whale_con_db;
a66.leo_whale_con_snr=leo_con_snr;

clearvars -except a66 a69 

for i_time=1:size(a69.data.spec_dB_5min,2)
    
    
    f=a69.data.frec;
    db_recorded=a69.data.spec_dB_5min(:,i_time);
    pxx=power(10,db_recorded/10);
    
    %%%% blue and fin whale contribution fit
    p=pxx;
    bio_spec= f<10 | f>15 & f<30  | f>50;
    p(bio_spec)=NaN;
    db=10*log10(p);
    
    f_points=logspace(1,4,100);
    db_points=interp1(f,db,f_points);
    [xData, yData] = prepareCurveData( f_points, db_points );
    ft = fittype( 'power2' );
    
    [fitresult, gof] = fit( xData, yData, ft);
    db_fit = feval(fitresult,f);
    p_fit=10.^(db_fit/10);
    
    %%%% blue whale contribution SPL
    con_abs(i_time)=   (sqrt( sum(  pxx(f>=26 & f<=28) ) )) ;
    con_noise(i_time)=   (sqrt( sum(  p_fit(f>=26 & f<=28) ) )) ;
    
    blue_con_snr(i_time)= 20*log10(con_abs(i_time)) - 20*log10(con_noise(i_time));
    if  blue_con_snr(i_time)>3
        blue_whale_con_p(i_time)= con_abs(i_time) - con_noise(i_time) ;
        blue_whale_con_db(i_time)=real(20*log10(blue_whale_con_p(i_time)));
    else
        blue_whale_con_p(i_time)=NaN;
        blue_whale_con_db(i_time)=NaN;
    end
    
    %%%% fin whale contribution    
    p=pxx;
    bio_spec= f<50 | f>97 & f<101  | f>150;
    p(bio_spec)=NaN;
    db=10*log10(p);
    
    f_points=logspace(1,4,100);
    db_points=interp1(f,db,f_points);
    [xData, yData] = prepareCurveData( f_points, db_points );
    ft = fittype( 'poly4' );
    
    [fitresult, gof] = fit( xData, yData, ft);
    db_fit = feval(fitresult,f);
    p_fit=10.^(db_fit/10);
    
    con_abs=   (sqrt( sum(  pxx(f>=96 & f<=99) ) )) ;
    con_noise=   (sqrt( sum(  p_fit(f>=96 & f<=99) ) )) ;
    
    fin_con_snr(i_time)= 20*log10(con_abs) - 20*log10(con_noise);
    if  fin_con_snr(i_time)>1
        fin_whale_con_p(i_time)= con_abs - con_noise ;
        fin_whale_con_db(i_time)=real(20*log10(fin_whale_con_p(i_time)));
    else
        fin_whale_con_p(i_time)=NaN;
        fin_whale_con_db(i_time)=NaN;
    end
    
    %%%% minke whale contribution    
    p=pxx;
    bio_spec= f<30 | f>97 & f<500  | f>1000;
    p(bio_spec)=NaN;
    db=10*log10(p);
    
    f_points=logspace(1,4,100);
    db_points=interp1(f,db,f_points);
    [xData, yData] = prepareCurveData( f_points, db_points );
    ft = fittype( 'power2' );
    
    [fitresult, gof] = fit( xData, yData, ft);
    db_fit = feval(fitresult,f);
    p_fit=10.^(db_fit/10);
    
    con_abs=   (sqrt( sum(  pxx(f>=105 & f<=300) ) )) ;
    con_noise=   (sqrt( sum(  p_fit(f>=105 & f<=300) ) )) ;
    
    minke_con_snr(i_time)= 20*log10(con_abs) - 20*log10(con_noise);
    if  minke_con_snr(i_time)>3
        minke_whale_con_p(i_time)= con_abs - con_noise ;
        minke_whale_con_db(i_time)=real(20*log10(minke_whale_con_p(i_time)));
    else
        minke_whale_con_p(i_time)=NaN;
        minke_whale_con_db(i_time)=NaN;
    end
    
    %%%% leopard seal contribution   
    con_abs=   (sqrt( sum(  pxx(f>=320 & f<=350) ) )) ;
    con_noise=   (sqrt( sum(  p_fit(f>=320 & f<=350) ) )) ;
    
    leo_con_snr(i_time)= 20*log10(con_abs) - 20*log10(con_noise);
    if  leo_con_snr(i_time)>3
        leo_whale_con_p(i_time)= con_abs - con_noise ;
        leo_whale_con_db(i_time)=real(20*log10(leo_whale_con_p(i_time)));
    else
        leo_whale_con_p(i_time)=NaN;
        leo_whale_con_db(i_time)=NaN;
    end
    
end

a69.time=a69.data.startdate;
a69.blue_whale_con_db=blue_whale_con_db;
a69.blue_whale_con_snr=blue_con_snr;
a69.fin_whale_con_db=fin_whale_con_db;
a69.fin_whale_con_snr=fin_con_snr;
a69.minke_whale_con_db=minke_whale_con_db;
a69.minke_whale_con_snr=minke_con_snr;
a69.leo_whale_con_db=leo_whale_con_db;
a69.leo_whale_con_snr=leo_con_snr;

clearvars -except a66 a69 

%% plot marine mammal contributions

%%% aural 66°S

[time,ia,ic]=unique(a66.time);
ia(1:(6*3))=[];
ia(end-(6*3):end)=[];

signal=a66.blue_whale_con_db(ia);
signalinterp=interp1(time(~isnan(signal)),signal(~isnan(signal)),a66.time,'nearest','extrap');

samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );
a66.lp_blue=filtfilt(b,a,signalinterp);


fincon_start=[7.334756e+05 7.338463333333334e+05 7.342003333333334e+05]
fincon_end= [datenum(2008,7,1) datenum(2009,7,1) datenum(2010,7,1)];

a66.lp_fin=nan([1 numel(a66.time)]);
    [time,ia,ic]=unique(a66.time);
    signal=a66.fin_whale_con_db(ia);
    
samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );

for i=1:3
    clear ix_time
    ix_time=time>fincon_start(i) & time<fincon_end(i);

    signalinterp=interp1(time(~isnan(signal) & ix_time'),signal(~isnan(signal)& ix_time'),a66.time,'nearest','extrap');
signal_filt=filtfilt(b,a,signalinterp);
a66.lp_fin(ix_time)=signal_filt(ix_time);

end

minkecon_start=[7.335325000000000e+05 7.338991666666666e+05 7.342801666666666e+05];
minkecon_end= [7.337543333333334e+05 7.341125000000000e+05 7.344701666666666e+05];


a66.lp_minke=nan([1 numel(a66.time)]);
    [time,ia,ic]=unique(a66.time);
    signal=a66.minke_whale_con_db(ia);
    
samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );

for i=1:3
    clear ix_time
    ix_time=time>minkecon_start(i) & time<minkecon_end(i);
signalinterp=interp1(time(~isnan(signal) & ix_time'),signal(~isnan(signal)& ix_time'),a66.time,'nearest','extrap');
signal_filt=filtfilt(b,a,signalinterp);
a66.lp_minke(ix_time)=signal_filt(ix_time);

end

leocon_start=[7.337426666666666e+05 734107 734467];
leocon_end= [7.337845000000000e+05 7.341458333333334e+05 7.344883333333334e+05];

a66.lp_leo=nan([1 numel(a66.time)]);
    [time,ia,ic]=unique(a66.time);
    signal=a66.leo_whale_con_db(ia);
    
samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );

for i=1:3
    clear ix_time
    ix_time=time>leocon_start(i) & time<leocon_end(i);
signalinterp=interp1(time(~isnan(signal) & ix_time'),signal(~isnan(signal)& ix_time'),a66.time,'nearest','extrap');
signal_filt=filtfilt(b,a,signalinterp);
a66.lp_leo(ix_time)=signal_filt(ix_time);

end

%%% aural 69°S

[time,ia,ic]=unique(a69.time);
ia(1:(6*3))=[];
ia(end-(6*3):end)=[];

signal=a69.blue_whale_con_db(ia);
signal(3869)=NaN;
signalinterp=interp1(time(~isnan(signal)),signal(~isnan(signal)),a69.time,'nearest','extrap');

samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );
a69.lp_blue=filtfilt(b,a,signalinterp);

fincon_start=[733500 7.338463333333334e+05 7.342003333333334e+05];
fincon_end= [733590 datenum(2009,7,1) datenum(2010,7,1)];

a69.lp_fin=nan([1 numel(a69.time)]);
    [time,ia,ic]=unique(a69.time);
    signal=a69.fin_whale_con_db(ia);
    
samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );

for i=1:3
    clear ix_time
    ix_time=time>fincon_start(i) & time<fincon_end(i);
signalinterp=interp1(time(~isnan(signal) & ix_time'),signal(~isnan(signal)& ix_time'),a69.time,'nearest','extrap');
signal_filt=filtfilt(b,a,signalinterp);
a69.lp_fin(ix_time)=signal_filt(ix_time);
end

minkecon_start=[7.335325000000000e+05 7.338991666666666e+05 7.342801666666666e+05];
minkecon_end= [7.337543333333334e+05 7.341125000000000e+05 7.344701666666666e+05];

a69.lp_minke=nan([1 numel(a69.time)]);
    [time,ia,ic]=unique(a69.time);
    signal=a69.minke_whale_con_db(ia);
    
samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );

for i=1:3
    clear ix_time
    ix_time=time>minkecon_start(i) & time<minkecon_end(i);
signalinterp=interp1(time(~isnan(signal) & ix_time'),signal(~isnan(signal)& ix_time'),a69.time,'nearest','extrap');
signal_filt=filtfilt(b,a,signalinterp);
a69.lp_minke(ix_time)=signal_filt(ix_time);

end

leocon_start=[7.337105000000000e+05 7.3409e+05 734467];
leocon_end= [7.337845000000000e+05 7.341458333333334e+05 7.344883333333334e+05];

a69.lp_leo=nan([1 numel(a69.time)]);
    [time,ia,ic]=unique(a69.time);
    signal=a69.leo_whale_con_db(ia);
    
samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );

for i=1:3
    clear ix_time
    ix_time=time>leocon_start(i) & time<leocon_end(i);
signalinterp=interp1(time(~isnan(signal) & ix_time'),signal(~isnan(signal)& ix_time'),a69.time,'nearest','extrap');
signal_filt=filtfilt(b,a,signalinterp);
a69.lp_leo(ix_time)=signal_filt(ix_time);
end


%%% prepare ticks

clear monthticks
i=1;
for iyear=[2008,2009,2010]
for imonth=1:12
monthticks(i)=datenum(iyear,imonth,1);
if imonth==1
    monthlabels{i}=datestr(monthticks(i),'yy-mm');
else
    monthlabels{i}=datestr(monthticks(i),'mm');
end

i=i+1;

end
end
monthlabels(monthticks<datenum(2008,3,1) | monthticks>datenum(2010,12,31))=[];
monthticks(monthticks<datenum(2008,3,1) | monthticks>datenum(2010,12,31))=[];

%%% compare species and recorders

figure(3)
set(gcf,'color',[1 1 1])
clf

subplot(411)
hold on
plot(a66.time,a66.lp_blue,'-r','linewidth',1.25)
plot(a69.time,a69.lp_blue,'-b','linewidth',1.25)

set(gca, 'Xtick',[datenum(2008,1,1):datenum(1,0,0):datenum(2011,1,1)],'XTickLabel', [])
% ylabel({'con intensity';'[db re 1 \muPa]'})
hold off
grid on
set(gca,'ylim',[60 110],'xlim',[datenum(2008,4,1) datenum(2010,12,31)])
text(datenum(2008,4,15),115,'a) Antarctic blue whales','fontweight','bold')
box on

%datetick('x','yyyy-mm','keeplimits')
set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'ylim',[60 110],'xlim',[datenum(2008,4,1) datenum(2010,12,31)])
set(gca,'xticklabelrotation',40) 

subplot(412)
hold on
plot(a66.time,a66.lp_fin,'-r','linewidth',1.25)
plot(a69.time,a69.lp_fin,'-b','linewidth',1.25)

set(gca, 'Xtick',[datenum(2008,1,1):datenum(1,0,0):datenum(2011,1,1)],'XTickLabel', [])
% ylabel({'con intensity';'[db re 1 \muPa]'})
hold off
grid on
set(gca,'ylim',[60 110],'xlim',[datenum(2008,4,1) datenum(2010,12,31)])
text(datenum(2008,4,15),115,'b) Fin whales','fontweight','bold')
box on
%datetick('x','yyyy-mm','keeplimits')
set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'ylim',[60 110],'xlim',[datenum(2008,4,1) datenum(2010,12,31)])
set(gca,'xticklabelrotation',40) 

subplot(413)
hold on
plot(a66.time,a66.lp_minke,'-r','linewidth',1.25)
plot(a69.time,a69.lp_minke,'-b','linewidth',1.25)

set(gca, 'Xtick',[datenum(2008,1,1):datenum(1,0,0):datenum(2011,1,1)],'XTickLabel', [])
% ylabel({'con intensity';'[db re 1 \muPa]'})
hold off
grid on
set(gca,'ylim',[60 110],'xlim',[datenum(2008,4,1) datenum(2010,12,31)])
text(datenum(2008,4,15),115,'c) Antarctic minke whales','fontweight','bold')
box on
%datetick('x','yyyy-mm','keeplimits')
set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'ylim',[60 110],'xlim',[datenum(2008,4,1) datenum(2010,12,31)])
set(gca,'xticklabelrotation',40) 

subplot(414)
hold on
plot(a66.time,a66.lp_leo,'-r','linewidth',1.25)
plot(a69.time,a69.lp_leo,'-b','linewidth',1.25)

set(gca, 'Xtick',[datenum(2008,1,1):datenum(1,0,0):datenum(2011,1,1)],'XTickLabel', [])
% ylabel({'con intensity';'[db re 1 \muPa]'})
hold off
grid on
text(datenum(2008,4,15),115,'d) Leopard seals','fontweight','bold')
box on
%datetick('x','yyyy-mm','keeplimits')
ylabel('SPL_M_M_C in dB re 1 \muPa')


set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'ylim',[60 110],'xlim',[datenum(2008,4,1) datenum(2010,12,31)])
set(gca,'xticklabelrotation',40) 

legend('66°S','69°S')

% set(gcf,'PaperPositionMode','auto')
% saveas(gcf,'marine_mammal_cones_spl_scale','fig')
%   print(gcf,'-dtiff ','marine_mammal_cones_spl_scale','-r400') 
%   print(gcf,'-dpdf ','marine_mammal_cones_spl_scale','-r400') 

%% compare marine mammal contribution SNR and SPL for minke whale contr.

[time,ia,ic]=unique(a66.time);
ia(1:(6*3))=[];
ia(end-(6*3):end)=[];

signal=a66.minke_whale_con_snr(ia);
signalinterp=interp1(time(~isnan(signal)),signal(~isnan(signal)),a66.time,'nearest','extrap');

samplerate=6;
[b,a] = butter(6, (1/(1*samplerate)) / samplerate );
a66.lp_minke_snr=filtfilt(b,a,signalinterp);


figure(4)
set(gcf,'color',[1 1 1])
clf

subplot(211)
hold on
box on
grid on
plot(a66.time,a66.minke_whale_con_snr,'.','color',[0.7 0.7 0.7])
plot(a66.time,a66.lp_minke_snr,'-k','linewidth',1.25)
ylabel({'Antarctic minke whale','SNR_M_M_C in dB re 1 \muPa'})
ylim([0 20])
% title('Minke whale - Signal to noise ratio')
 xlim([datenum(2008,4,1) datenum(2010,12,31)])
 text(datenum(2008,4,15),22,'a)','fontweight','bold')

set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'xticklabelrotation',40) 
% 

subplot(212)
hold on
box on
grid on
plot(a66.time,a66.minke_whale_con_db,'.','color',[0.7 0.7 0.7])
plot(a66.time,a66.lp_minke,'-k','linewidth',1.25)
ylabel({'Antarctic minke whale','SPL_M_M_C in dB re 1 \muPa'})
ylim([75 115])
%title('Minke whale con - Sound pressure level')
legend('Recordings in 4h intervals','1 day low pass filtered')
 xlim([datenum(2008,4,1) datenum(2010,12,31)])
 text(datenum(2008,4,15),119,'b)','fontweight','bold')

set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'xticklabelrotation',40) 

% set(gcf,'PaperPositionMode','auto')
% saveas(gcf,'minke_con_snr_spl_comparison','fig')
%  print(gcf,'-dtiff','minke_con_snr_spl_comparison','-r400') 
%  print(gcf,'-dpdf ','minke_con_snr_spl_comparison','-r400') 

clearvars -except a66 a69 

%% Minke whale contribution diel cycle plot

% requires function: shadedErrorBar from Rob Campbell - November 2009
% available on Matlab file exchange 
% Link: http://www.mathworks.com/matlabcentral/fileexchange/26311-shadederrorbar

% find Minke whale SPL for each day and each hour of day

ix_minke_days=a66.time>datenum(2008,4,20) & a66.time<datenum(2008,8,10) | a66.time>datenum(2009,4,20) & a66.time<datenum(2009,8,10) | a66.time>datenum(2010,4,20) & a66.time<datenum(2010,8,10) ;
minke_days=a66.minke_whale_con_db;

fs=6; % rec per day

starttime_loop=datenum(2008,4,1,12,00,00);
endtime_loop=datenum(2010,11,1,12,00,00);

clear diel_matrix_spl diel_matrix_time
i=1;
for timeinterval=starttime_loop:datenum(0,0,1):endtime_loop
    
    timetofind=timeinterval:datenum(0,0,0,4,0,0):timeinterval+datenum(0,0,0,24,0,0);
    %     datestr(timetofind)
    ix=find(ismember(a66.time,timetofind));
    
    diel_matrix_spl(i,:)=nan(1,7);
    dv= datevec(timetofind);
    diel_matrix_time(i,:)=datenum(0,0,0,dv(:,4),dv(:,5),dv(:,6));
    diel_matrix_date(i)=datenum(dv(4,1),dv(4,2),dv(4,3));
    
    diel_matrix_spl(i, ismember(timetofind,a66.time(ix)) )=minke_days(ix);
    
    i=i+1;
end

for i=1:size(diel_matrix_spl,1)

    a=diel_matrix_spl(i,:);
    b= (a - min(a)) ./ (max(a) -min(a)) ;
    
diel_matrix_spl_norm(i,:)=b;

end

% diel_matrix_spl_norm=normr(diel_matrix_spl);
diel_matrix_spl_norm(diel_matrix_spl_norm==1)=NaN;

% generate idealized dvm pattern
hour=0:4:20;
dvm=[1 0.9 0.1 0 0.1 0.9 ];
dv=datevec(a66.time);
a66.dvm=nan([1 numel(a66.time)]);

for i=1:6
    ix= dv(:,4)==hour(i);
    a66.dvm(ix)=dvm(i);
end

% calculate the contributions statistics over a 24h cycle

dv=datevec(diel_matrix_date);
ix_time=dv(:,2)>=5 & dv(:,2)<=7 ;
mean_dailynormspl_2008_dvm=nanmean(diel_matrix_spl_norm(ix_time,:),1);
std_dailynormspl_2008_dvm=nanstd(diel_matrix_spl_norm(ix_time,:),1);

mean_dailynormspl_2008_rand=nanmean(diel_matrix_spl_norm(~ix_time,:),1);
std_dailynormspl_2008_rand=nanstd(diel_matrix_spl_norm(~ix_time,:),1);

box_dvm_all=nan([4 numel(a66.time)]);
box_dvm_may_july=nan([4 numel(a66.time)]);
box_combi=nan([8 numel(a66.time)]);

dv=datevec(a66.time);
ix_time_may=dv(:,2)>=5 & dv(:,2)<=7 ;

levels=[0,0.1,0.9,1];
for i=1:4
    ix= a66.dvm==levels(i);
    ix_other=ix & ~ix_time_may';
    box_dvm_all(i,ix_other)=a66.minke_whale_con_db(ix_other);
    
    ix_mayjuly=ix & ix_time_may';
    box_dvm_may_july(i,ix_mayjuly)=a66.minke_whale_con_db(ix_mayjuly);   
end

boxcombi=[box_dvm_all;box_dvm_may_july];
levels={'0','0.1','0.9','1','0','0.1','0.9','1'}
group={'Other','Other','Other','Other','May-July','May-July','May-July','May-July'}

% plot the results
figure(5)
clf
set(gcf,'color','w')

subplot(211)
grid on
box on
hold on

shadedErrorBar(1:7,mean_dailynormspl_2008_rand,std_dailynormspl_2008_rand,'-b',1)

shadedErrorBar(1:7,mean_dailynormspl_2008_dvm,std_dailynormspl_2008_dvm,'-r',1)

p1=plot(1:7,mean_dailynormspl_2008_dvm,'-r')

p2= plot(1:7,mean_dailynormspl_2008_rand,'-b')


set(gca,'xtick',[1:7],'xticklabel',{'12:00','16:00','20:00','00:00','04:00','08:00','12:00'})
ylabel({'Normalized SPL_M_M_C'})
legend([p1,p2],'May-July','August-November')
ylim([0 1])

text(1,1.1,'a)','fontweight','bold')

subplot(212)
grid on
box on
hold on
boxplot(boxcombi',{levels,group},'colors',repmat('br',1,4),'factorgap',[50 0],'plotstyle','compact','symbol','');
set(gca,'XTick',[1.5,18],'XTickLabel',{'      12:00    \newlineZoopl. at depth','        00:00   \newlineZoopl. at surface'},'xgrid','Off')
ylim([75 115])
xlim([0 20])
ylabel({'SPL_M_M_C in dB re 1 \muPa'})
text(.1,119,'b)','fontweight','bold')

%  set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dtiff','minke_dvm_relation','-r400')
%   print(gcf,'-dpdf','minke_dvm_relation','-r400')
% saveas(gcf,'minke_dvm_relation_2','fig')

% plot supplementary figure 1
% minke whale chorus diel rythm over time

ix_useful=sum(isnan(diel_matrix_spl_norm),2)<3;
ixx=diel_matrix_date>datenum(2008,12,1) & diel_matrix_date<datenum(2009,5,1);
ix_useful(ixx)=1;
ixx=diel_matrix_date>datenum(2009,12,1) & diel_matrix_date<datenum(2010,5,1);
ix_useful(ixx)=1;

clear monthticks
i=1;
for iyear=[2008,2009,2010]
    for imonth=1:12
        monthticks(i)=datenum(iyear,imonth,1);
        if imonth==1
            monthlabels{i}=datestr(monthticks(i),'yy-mm');
        else
            monthlabels{i}=datestr(monthticks(i),'mm');
        end
        
        i=i+1;
        
    end
end
monthlabels(monthticks<datenum(2008,3,1) | monthticks>datenum(2010,12,31))=[];
monthticks(monthticks<datenum(2008,3,1) | monthticks>datenum(2010,12,31))=[];


figure(11)
clf
contourf(diel_matrix_date(ix_useful),1:7,diel_matrix_spl_norm(ix_useful,:)',10,'edgecolor','none')
colormap('gray')
colormap(flipud(colormap))
cb=colorbar
set(gca,'clim',[0.37 0.4])
set(gca,'ytick',[1,4,7],'yticklabel',{'12:00','00:00','12:00'},'tickdir','out')
cb=colorbar
ylabel(cb,'Normalized SPL')
grid on
box on
xlim([datenum(2008,4,1) datenum(2010,12,31)])

set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'xticklabelrotation',40) 


%% load sea ice concentration grid and calculate circle indicies
clearvars -except a66 a69 

url = 'http://www.iup.uni-bremen.de/seaice/amsredata/asi_daygrid_swath/l1a/s6250/grid_coordinates/LongitudeLatitudeGrid-s6250-Antarctic.hdf';
filename_grid = 'LongitudeLatitudeGrid-s6250-Antarctic.hdf';
outfilename = websave(filename_grid,url);

lon = hdfread(filename_grid, 'Longitudes');
lat = hdfread(filename_grid, 'Latitudes');

lon=lon-180;

latlim=double([min(min(lat)),max(max(lat))]);
lonlim=double([min(min(lon)),max(max(lon))]);

radius_steps_km=[50:50:2000];
center_point.lat=[-66,-69];
center_point.lon=[0,0];
r_earth=6371 ; %km

% find pixels in concentric circles relativ to recorder location

for i_point=1:numel(center_point.lat)
    
    for i_lat=1:size(lat,1)
        for i_lon=1:size(lat,2)            
            lat2=lat(i_lat,i_lon);
            lon2=lon(i_lat,i_lon);            
            circle_distance(i_lat,i_lon)= deg2km(distance(center_point.lat(i_point),center_point.lon(i_point),lat2,lon2));
        end
    end
    
    for i_r=1:numel(radius_steps_km)
        in_circle_grid(i_point,i_r,:,:)=circle_distance<radius_steps_km(i_r);
    end
end

%% load daily sea ice data from Uni Bremen
% 'http.//www.iup.uni-bremen.de/seaice/amsredata/asi_daygrid_swath/l1a/s6250/2008/jan/asi-s6250-20080101-v5.hdf'

timecounter=1;

for i_day=datenum(2008,04,1):datenum(2011,01,1)
    
month=lower(datestr(i_day,'mmm'));
year=datestr(i_day,'yyyy');
date=datestr(i_day,'yyyymmdd');

filename_daily=['http://www.iup.uni-bremen.de/seaice/amsredata/asi_daygrid_swath/l1a/s6250/',year,'/',month,'/asi-s6250-',date,'-v5.hdf'];
sea_ice_file='daily_sea_ice_file.hdf';

outfilename = websave(sea_ice_file,filename_daily);

ice = hdfread(sea_ice_file, 'ASI Ice Concentration');

for i_point=1:numel(center_point.lat)
    for i_r=1:numel(radius_steps_km)
        circle_sum_ice(i_point,timecounter,i_r)=nansum(ice( squeeze(in_circle_grid(i_point,i_r,:,:)) & ~isnan(ice)));
        circle_area_ice(i_point,timecounter,i_r)=numel(ice( squeeze(in_circle_grid(i_point,i_r,:,:)) & ~isnan(ice)));
    end
end

icecon.time(timecounter)=datenum(i_day);

disp(['sea ice conc at ' datestr(icecon.time(timecounter)) ' !'])
clear ice    
timecounter=timecounter+1;
end    

circle_mean_time=icecon.time;

mean_ice=circle_sum_ice./circle_area_ice;

radius_steps_km_ice=radius_steps_km;

%% plot SPL and sea ice concentration

a66.ice=interp1(circle_mean_time,mean_ice(1,:,4),a66.data.startdate);
a69.ice=interp1(circle_mean_time,mean_ice(2,:,4),a69.data.startdate);

Fs = 6;  % Sampling Frequency
windowsize=Fs*10;
b=ones(windowsize,1)./windowsize;
a=1;
splscale=90:0.5:125;

figure(6)
clf
set(gcf,'color','w')

subplot(221)
hold on
plot(a66.data.startdate,a66.data.spl_broadband_5min,'.','color',[0.7 0.7 0.7])
plot(a66.data.startdate,filtfilt(b,a,a66.data.spl_broadband_5min),'-k','linewidth',2)
set(gca,'xlim',[a66.data.startdate(50) a66.data.startdate(end-50)],'ylim',[100 125])
datetick('x','yyyy-mm','keeplimits')
ylabel('SPL_{RMS} [dB re 1 \muPa]')
box on
grid on

subplot(222)
hold on
histval=histc(a66.data.spl_broadband_5min,splscale);
histval=histval./numel(a66.data.spl_broadband_5min)*100;
area(splscale,histval,'facecolor',[0.7 0.7 0.7],'edgecolor','none');

histval=histc(a66.data.spl_broadband_5min(a66.ice>50),splscale);
histval=histval./numel(a66.data.spl_broadband_5min)*100;

plot(splscale,histval,'-k','linewidth',2)

histval=histc(a66.data.spl_broadband_5min(a66.ice<50),splscale);
histval=histval./numel(a66.data.spl_broadband_5min)*100;

plot(splscale,histval,'--k','linewidth',2)
view(90,90)

set(gca,'xlim',[100 125],'xdir','reverse')
box on
grid on

%-------------

subplot(223)
hold on
plot(a69.data.startdate,a69.data.spl_broadband_5min,'.','color',[0.7 0.7 0.7])
plot(a69.data.startdate,filtfilt(b,a,a69.data.spl_broadband_5min),'-k','linewidth',2)
set(gca,'xlim',[a69.data.startdate(50) a69.data.startdate(end-50)],'ylim',[100 125],'xtick',[datenum(2009,1,1),datenum(2010,1,1)])
datetick('x','yyyy-mm','keeplimits')
box on
grid on

subplot(224)
hold on
histval=histc(a69.data.spl_broadband_5min,splscale);
histval=histval./numel(a69.data.spl_broadband_5min)*100;
area(splscale,histval,'facecolor',[0.7 0.7 0.7],'edgecolor','none');

histval=histc(a69.data.spl_broadband_5min(a69.ice>50),splscale);
histval=histval./numel(a69.data.spl_broadband_5min)*100;

plot(splscale,histval,'-k','linewidth',2)

histval=histc(a69.data.spl_broadband_5min(a69.ice<50),splscale);
histval=histval./numel(a69.data.spl_broadband_5min)*100;

plot(splscale,histval,'--k','linewidth',2)

view(90,90)

set(gca,'xlim',[100 125],'xdir','reverse')
box on
grid on
ylabel('Normalized distribution [%]')

legend('All recordings','Sea ice conc. > 50 %','Sea ice conc. < 50 %')

%  set(gcf,'PaperPositionMode','auto')
%  print(gcf,'-dpdf','spl_time_series_and_histogram','-r400') 
%  print(gcf,'-dtiff','spl_time_series_and_histogram','-r400') 

%% time series plot comparing psd, sea ice cocentration, draft and extent

% load sea ice draft data

url='https://doi.pangaea.de/10.1594/PANGAEA.821257?format=textfile';
filename = 'AWI232-9_sea_ice_draft.tab';
outfilename = websave(filename,url);

filename = outfilename;
delimiter = '\t';
startRow = 25;
formatSpec = '%q%*q%*q%*q%*q%*q%f%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

time = dataArray{:, 1};
draft_model = dataArray{:, 2};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

sea_ice_draft_awi232_9.time=datenum(time,'yyyy-mm-ddTHH:MM');
sea_ice_draft_awi232_9.draft_model=draft_model;

%% load sea ice extent data
filename='DATASETS/NOAA/G02135/south/daily/data/SH_seaice_extent_final.csv';

mw = ftp('sidads.colorado.edu');
mget(mw, filename);
close(mw);

delimiter = ',';
startRow = 3;
formatSpec = '%f%f%f%f%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

Year = dataArray{:, 1};
Month = dataArray{:, 2};
Day = dataArray{:, 3};
Extent = dataArray{:, 4};

sea_ice_index.time=datenum(Year,Month,Day);
sea_ice_index.extent=Extent;

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% plot PSD and sea ice variables

ix_frec=a66.data.frec>30 & a66.data.frec<80;
a66.lf=mean(a66.data.spec_dB_5min(ix_frec,:),1);
% ix_frec=a66.data.frec>500 & a66.data.frec<10000;
ix_frec=a66.data.frec>500 & a66.data.frec<1000;
a66.hf=mean(a66.data.spec_dB_5min(ix_frec,:),1);

a66.lf_interp =interp1(a66.data.startdate,a66.lf,circle_mean_time) ;
a66.hf_interp =interp1(a66.data.startdate,a66.hf,circle_mean_time) ;

ix_frec=a69.data.frec>30 & a69.data.frec<80;
a69.lf=mean(a69.data.spec_dB_5min(ix_frec,:),1);
% ix_frec=a69.data.frec>500 & a69.data.frec<10000;
ix_frec=a69.data.frec>500 & a69.data.frec<1000;
a69.hf=mean(a69.data.spec_dB_5min(ix_frec,:),1);

a69.lf_interp =interp1(a69.data.startdate,a69.lf,circle_mean_time) ;
a69.hf_interp =interp1(a69.data.startdate,a69.hf,circle_mean_time) ;

samplerate=1;
[b,a] = butter(6, (1/(20*samplerate)) / samplerate );

%%%

clear monthticks
i=1
for iyear=[2008,2009,2010]
for imonth=1:12
monthticks(i)=datenum(iyear,imonth,1);
if imonth==1
    monthlabels{i}=datestr(monthticks(i),'yy-mm');
else
    monthlabels{i}=datestr(monthticks(i),'mm');
end

i=i+1;

end
end
monthlabels(monthticks<datenum(2008,3,1) | monthticks>datenum(2010,12,31))=[];
monthticks(monthticks<datenum(2008,3,1) | monthticks>datenum(2010,12,31))=[];

figure(7)
clf
set(gcf,'color','w')

subplot(511)
grid on
box on
hold on


plot(circle_mean_time,a66.lf_interp,'.','color',[1 0.7 0.7])
plot(circle_mean_time,a69.lf_interp,'.','color',[0.7 0.7 1])

ix_notnan=~isnan(a66.lf_interp);
p1=plot(circle_mean_time(ix_notnan),filtfilt(b,a,a66.lf_interp(ix_notnan)),'-r')
p2=plot(circle_mean_time(ix_notnan),filtfilt(b,a,a69.lf_interp(ix_notnan)),'-b')
ylabel('PSD in dB re 1 \muPa')
legend([p1,p2],'66°S','69°S')
xlim([datenum(2008,4,1) datenum(2010,12,31)])
ylim([65 95])
set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'xticklabelrotation',40)

text(datenum(2008,4,15),98,'a)','fontweight','bold')

subplot(512)
grid on
box on
hold on


plot(circle_mean_time,a66.hf_interp,'.','color',[1 0.7 0.7])
plot(circle_mean_time,a69.hf_interp,'.','color',[0.7 0.7 1])

ix_notnan=~isnan(a66.lf_interp);
p1=plot(circle_mean_time(ix_notnan),filtfilt(b,a,a66.hf_interp(ix_notnan)),'-r')
p2=plot(circle_mean_time(ix_notnan),filtfilt(b,a,a69.hf_interp(ix_notnan)),'-b')
ylabel('PSD in dB re 1 \muPa')
legend([p1,p2],'66°S','69°S')
xlim([datenum(2008,4,1) datenum(2010,12,31)])
ylim([45 80])
set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'xticklabelrotation',40)

text(datenum(2008,4,15),72,'b)','fontweight','bold')

subplot(513)
grid on
box on
hold on

plot(circle_mean_time,mean_ice(2,:,10),'.','color',[0.7 0.7 1])
plot(circle_mean_time,mean_ice(1,:,10),'.','color',[1 0.7 0.7])

p1=plot(circle_mean_time,filtfilt(b,a,double(mean_ice(1,:,10))),'-r')
p2=plot(circle_mean_time,filtfilt(b,a,double(mean_ice(2,:,10))),'-b')
ylabel('Sea ice conc. in %')
legend([p1,p2],'66°S','69°S')
xlim([datenum(2008,4,1) datenum(2010,12,31)])
set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'xticklabelrotation',40)
ylim([0 100])
text(datenum(2008,4,15),1.1,'c)','fontweight','bold')

subplot(515)

grid on
box on
hold on

plot(sea_ice_index.time,sea_ice_index.extent,'-k')

ylabel('Sea ice extent in km^2')
xlim([datenum(2008,4,1) datenum(2010,12,31)])
set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'xticklabelrotation',40)
% xticklabel_rotate([],45)
text(datenum(2008,4,15),21,'e)','fontweight','bold')

subplot(514)
grid on
box on
hold on

draft_interp_awi232_9=interp1(sea_ice_draft_awi232_9.time,sea_ice_draft_awi232_9.draft_model,circle_mean_time)

plot(circle_mean_time,draft_interp_awi232_9,'.','color',[0.7 0.7 1])

ix_notnan=~isnan(draft_interp_awi232_9);
p=plot(circle_mean_time(ix_notnan),filtfilt(b,a,draft_interp_awi232_9(ix_notnan)),'-b')

ylabel('Sea ice draft in m')
xlim([datenum(2008,4,1) datenum(2010,12,31)])
ylim([0 5])
set(gca,'xtick',monthticks)
set(gca,'xticklabel',monthlabels)
set(gca,'xticklabelrotation',40)
legend(p,'69°S')
text(datenum(2008,4,15),5.5,'d)','fontweight','bold')

%  set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dtiff','psd_time_series_30-80hz_and_sea_ice','-r400') 
%   print(gcf,'-dpdf','psd_time_series_30-80hz_and_sea_ice','-r400') 

%% load wind speed data

% get ERA-INTERIM reanalysis wind speed data (u10, v10) from ECMWF website
% http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/
% personal account (free) necessary for downloading data

% insert name of downloaded nc file here and unncomment:
% nc_file='netcdf-atls05-20141001165205-429-0777.nc';

info=ncinfo(nc_file)

for i=[1:5]    
    str=[ 'era_interim.' info.Variables(i).Name '=ncread(nc_file,''' info.Variables(i).Name  ''');'];
    eval(str);  
end

[era_interim.lat,era_interim.lon] = meshgrid(era_interim.latitude,era_interim.longitude);
% [Z, refvec] = geoloc2grid(era_interim.lat,era_interim.lon,era_interim.u10(:,:,1), 0.25);

latlim=double( [min(era_interim.latitude),max(era_interim.latitude)] );
lonlim=double( [min(era_interim.longitude),max(era_interim.longitude)] );
load coast

radius_steps_km=[50:50:1000];
center_point.lat=[-66,-69];
center_point.lon=[0,0];

r_earth=6371 ; %km

% calculate absolute wind speed

era_interim.uv= sqrt( era_interim.u10.^2 + era_interim.v10.^2 );
rmfield(era_interim,{'u10','v10'})

%% load topography and check for ocean and land

url='https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/binary/etopo1_ice_c_f4.zip';
filename='etopo1_ice_c_f4.zip';
websave(filename,url)
unzip(filename)

[etopo_z, etopo_r] = etopo('etopo1_ice_c_f4.flt', 1, latlim, lonlim);
etopo_ocean=etopo_z<=0;
era_interim.isocean = logical( ltln2val(etopo_ocean, etopo_r, era_interim.lat, era_interim.lon) );
% [isocean_z,isocean_r] = geoloc2grid(double(era_interim.lat),double(era_interim.lon),double( era_interim.isocean(:,:) ), 0.25);

clear etopo_ocean etopo_z etopo_r

%% calculate wind speed averages in circles around recorder location


timecounter=1;
for i_point=1:numel(center_point.lat)
    
    for i_lat=1:size(era_interim.latitude,1)
        for i_lon=1:size(era_interim.longitude,1)
            
            lat2=era_interim.latitude(i_lat);
            lon2=era_interim.longitude(i_lon);
            era_interim.circle_distance(i_lat,i_lon)= deg2km(distance(center_point.lat(i_point),center_point.lon(i_point),lat2,lon2));
            
        end
    end
    
    for i_r=1:numel(radius_steps_km)
        era_interim.in_circle_grid(i_point,i_r,:,:)=era_interim.circle_distance<radius_steps_km(i_r);
    end
    
end


for i_day=1:numel(era_interim.time)

    wind=era_interim.uv(:,:,i_day);
    wind(~era_interim.isocean)=NaN;
    
for i_point=1:numel(center_point.lat)  
     for i_r=1:numel(radius_steps_km)
era_interim.circle_sum_uv(i_point,timecounter,i_r)=nansum(wind( squeeze(era_interim.in_circle_grid(i_point,i_r,:,:))' & ~isnan(wind)));
era_interim.circle_area_uv(i_point,timecounter,i_r)=numel(wind( squeeze(era_interim.in_circle_grid(i_point,i_r,:,:))' & ~isnan(wind)));
     end   
end  

clear wind    
timecounter=timecounter+1;
end    

era_interim.circle_mean_time=linspace(datenum(2008,4,1),datenum(2011,1,1),numel(era_interim.time));

era_interim.mean_circle_uv=era_interim.circle_sum_uv./era_interim.circle_area_uv;
 
%% downsample wind data to match ice

mean_circle_time=circle_mean_time; 
mean_circle_uv_int(1,:) =interp1(era_interim.circle_mean_time,era_interim.mean_circle_uv(1,:,1),circle_mean_time) ;
mean_circle_uv_int(2,:) =interp1(era_interim.circle_mean_time,era_interim.mean_circle_uv(2,:,1),circle_mean_time) ;

%% PSD, wind speed and sea ice cocentration scatterplot

figure(8)
clf
set(gcf,'color','w')

subplot(221)
hold on
grid on
box on

dv=datevec(circle_mean_time);
month=dv(:,2);
cmp=parula(12);
for i=1:12
    ix=month==i;
    plot(mean_ice(1,ix,10),a66.lf_interp(ix),'.','color',cmp(i,:))
    plot(mean_ice(1,ix,10),a66.hf_interp(ix),'o','color',cmp(i,:))
end

xlabel('Sea ice concentration in %')
ylabel('PSD in dB re 1 \muPa^2 Hz^-^1')
legend('30-80 Hz','500-1000 Hz')
ylim([40 100])
text(0.05,103,'a) 66°S','fontweight','bold')
colormap(cmp)
set(gca,'clim',[0 1])
colorbar('ticks',linspace(0,1,12),'Ticklabels',{'Jan','Feb','Mar','Apr','Mai','Jun','Jul','Aug','Sep','Oct','Nov','Dez'})


subplot(223)
hold on
grid on
box on
dv=datevec(circle_mean_time);   
 month=dv(:,2);
cmp=parula(12);
for i=1:12
    ix=month==i;
plot(mean_ice(2,ix,10),a69.lf_interp(ix),'.','color',cmp(i,:))
plot(mean_ice(2,ix,10),a69.hf_interp(ix),'o','color',cmp(i,:))
end
xlabel('Sea ice concentration in %')
ylabel('PSD in dB re 1 \muPa^2 Hz^-^1')
legend off
ylim([40 100])
text(0.05,103,'b) 69°S','fontweight','bold')

colormap(cmp)
set(gca,'clim',[0 1])
colorbar('ticks',linspace(0,1,12),'Ticklabels',{'Jan','Feb','Mar','Apr','Mai','Jun','Jul','Aug','Sep','Oct','Nov','Dez'})

subplot(222)
hold on
grid on
box on

ix_ice_low=mean_ice(1,:,10)<50;
ix_ice_high=mean_ice(1,:,10)>50;
plot(mean_circle_uv_int(1,ix_ice_high),a66.lf_interp(ix_ice_high),'.b')
plot(mean_circle_uv_int(1,ix_ice_low),a66.lf_interp(ix_ice_low),'.r')

plot(mean_circle_uv_int(1,ix_ice_high),a66.hf_interp(ix_ice_high),'ob')
plot(mean_circle_uv_int(1,ix_ice_low),a66.hf_interp(ix_ice_low),'or')

xlabel('Wind speed in m s^-^1')
ylabel('PSD in dB re 1 \muPa^2 Hz^-^1')
legend('30-80 Hz, Ice con.>50%','30-80 Hz, Ice con.<50%','500-1000 Hz, Ice con.>50%','500-1000 Hz, Ice con.<50%')
ylim([40 100])
xlim([0 20])
text(1,103,'c) 66°S','fontweight','bold')

subplot(224)
hold on
grid on
box on

ix_ice_low=mean_ice(2,:,10)<50;
ix_ice_high=mean_ice(2,:,10)>50;
plot(mean_circle_uv_int(2,ix_ice_high),a69.lf_interp(ix_ice_high),'.b')
plot(mean_circle_uv_int(2,ix_ice_low),a69.lf_interp(ix_ice_low),'.r')

plot(mean_circle_uv_int(2,ix_ice_high),a69.hf_interp(ix_ice_high),'ob')
plot(mean_circle_uv_int(2,ix_ice_low),a69.hf_interp(ix_ice_low),'or')

xlabel('Wind speed in m s^-^1')
ylabel('PSD in dB re 1 \muPa^2 Hz^-^1')
legend off
ylim([40 100])
xlim([0 20])

text(1,103,'d) 69°S','fontweight','bold')


%  set(gcf,'PaperPositionMode','auto')
%  print(gcf,'-dpdf','scatterplot_wind_ice','-r400') 
%   print(gcf,'-dtiff','scatterplot_wind_ice','-r400') 

%% plot schematic spectra

local_era_interim=interp1(era_interim.circle_mean_time,era_interim.mean_circle_uv(1,:,1),a66.data.startdate);
local_ice=interp1(circle_mean_time,mean_ice(1,:,1),a66.data.startdate);

spldb=a66.data.spl_broadband_5min;
SPL_scale = min(spldb):0.25:max(spldb);
[counts_per_bin,bin] = histc(spldb,SPL_scale);
cumulative_cpb = cumsum(counts_per_bin);
cumulative_cpb=cumulative_cpb/numel(spldb);
[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;

figure(9)
set(gcf,'Color',[1 1 1]);
clf
hold on
box on
grid on
era_interimowsize=50;
view(0,90)
xlim([10 1000])
set(gca,'xscale','log')

ylabel('PSD [db re 1 \muPa^2 Hz^{-1}]')
xlabel('Frequency [Hz]')

ilow=find(a66.data.frec==10);
ihigh=find(a66.data.frec==1000);

y=nanmean(a66.data.spec_dB_5min(:,local_era_interim>0 & local_era_interim<10 & local_ice<50),2);
y=filtfilt(ones(1,era_interimowsize)/era_interimowsize,1,y(ilow:ihigh));
h.era_interim_noice=plot3(a66.data.frec(ilow:ihigh),y,ones(numel(a66.data.frec(ilow:ihigh)),1)-10,'--k','linewidth',1.5)

y=nanmean(a66.data.spec_dB_5min(:,local_era_interim>10 & local_era_interim<20 & local_ice<50),2);
y=filtfilt(ones(1,era_interimowsize)/era_interimowsize,1,y(ilow:ihigh));
plot3(a66.data.frec(ilow:ihigh),y,ones(numel(a66.data.frec(ilow:ihigh)),1)-10,'--k','linewidth',1.5)

y=nanmean(a66.data.spec_dB_5min(:,local_era_interim>20 & local_era_interim<30 & local_ice<50),2);
y=filtfilt(ones(1,era_interimowsize)/era_interimowsize,1,y(ilow:ihigh));
plot3(a66.data.frec(ilow:ihigh),y,ones(numel(a66.data.frec(ilow:ihigh)),1)-10,'--k','linewidth',1.5)

y=nanmean(a66.data.spec_dB_5min(:,local_era_interim>0 & local_era_interim<10 & local_ice>50),2);
y=filtfilt(ones(1,era_interimowsize)/era_interimowsize,1,y(ilow:ihigh));
plot3(a66.data.frec(ilow:ihigh),y,ones(numel(a66.data.frec(ilow:ihigh)),1)-10,'-k','linewidth',1.5)

y=nanmean(a66.data.spec_dB_5min(:,local_era_interim>10 & local_era_interim<20 & local_ice>50),2);
y=filtfilt(ones(1,era_interimowsize)/era_interimowsize,1,y(ilow:ihigh));
h.era_interim_ice=plot3(a66.data.frec(ilow:ihigh),y,ones(numel(a66.data.frec(ilow:ihigh)),1)-10,'-k','linewidth',1.5)

per=prctile(a66.data.spec_dB_5min,[1 99],2);

ylow=per(:,1);
ylow=filtfilt(ones(1,era_interimowsize)/era_interimowsize,1,ylow(ilow:ihigh));

yhigh=per(:,2);
yhigh=filtfilt(ones(1,era_interimowsize)/era_interimowsize,1,yhigh(ilow:ihigh));

xxx=[a66.data.frec(ilow:ihigh) ;  flipud(a66.data.frec(ilow:ihigh)) ];
yyy=[ylow ;  flipud(yhigh)];
ccc=[0.9 0.9 0.9];
fillhandle=patch(xxx,yyy,ones(1,numel(yyy))-20,ccc);
set(fillhandle,'Facecolor',ccc,'EdgeColor','none');

% blue

[~,ix_min]=min(abs(cumulative_cpb-0.10));
ix_percentile=bin==ix_min;
ylow=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);

[~,ix_min]=min(abs(cumulative_cpb-0.9));
ix_percentile=bin==ix_min;
yhigh=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);

ilow=find(a66.data.frec==15);
ihigh=find(a66.data.frec==30);

xxx=[a66.data.frec(ilow:ihigh) ;  flipud(a66.data.frec(ilow:ihigh)) ];
yyy=[ylow(ilow:ihigh) ;  flipud(yhigh(ilow:ihigh))];
ccc=[0 0.5 1];
fillhandle=patch(xxx,yyy,ccc);
set(fillhandle,'Facecolor',ccc,'EdgeColor','none');
h.blue=fillhandle

[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;
y=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);
plot3(a66.data.frec(ilow:ihigh),y(ilow:ihigh),ones(numel(a66.data.frec(ilow:ihigh)),1),'-','linewidth',1.5,'color',[0 0 1])

% fin

[~,ix_min]=min(abs(cumulative_cpb-0.10));
ix_percentile=bin==ix_min;
ylow=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);

[~,ix_min]=min(abs(cumulative_cpb-0.9));
ix_percentile=bin==ix_min;
yhigh=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);

ilow=find(a66.data.frec==93);
ihigh=find(a66.data.frec==103);

xxx=[a66.data.frec(ilow:ihigh) ;  flipud(a66.data.frec(ilow:ihigh)) ];
yyy=[ylow(ilow:ihigh) ;  flipud(yhigh(ilow:ihigh))];
ccc=[0 1 0.5];

fillhandle=patch(xxx,yyy,ccc);
set(fillhandle,'Facecolor',ccc,'EdgeColor','none');
h.fin=fillhandle

[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;
y=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);
plot3(a66.data.frec(ilow:ihigh),y(ilow:ihigh),ones(numel(a66.data.frec(ilow:ihigh)),1),'-','linewidth',1.5,'color',[0 0.5 0])

% minke
ilow=find(a66.data.frec==103);
ihigh=find(a66.data.frec==350);

spldb=nanmean(a66.data.spec_dB_5min(ilow:ihigh,:));
SPL_scale = min(spldb):0.25:max(spldb);
[counts_per_bin,bin] = histc(spldb,SPL_scale);
cumulative_cpb = cumsum(counts_per_bin);
cumulative_cpb=cumulative_cpb/numel(spldb);
[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;

[~,ix_min]=min(abs(cumulative_cpb-0.1));
ix_percentile=bin==ix_min;
ylow=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);
ylow=filtfilt(ones(1,5)/5,1,ylow);

[~,ix_min]=min(abs(cumulative_cpb-0.9));
ix_percentile=bin==ix_min;
yhigh=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);
yhigh=filtfilt(ones(1,5)/5,1,yhigh);

xxx=[a66.data.frec(ilow:ihigh) ;  flipud(a66.data.frec(ilow:ihigh)) ];
yyy=[ylow(ilow:ihigh) ;  flipud(yhigh(ilow:ihigh))];
ccc=[1 1 0];

fillhandle=patch(xxx,yyy,ccc);
set(fillhandle,'Facecolor',ccc,'EdgeColor','none');
h.minke=fillhandle

[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;
y=nanmean(a66.data.spec_dB_5min(:,ix_percentile),2);
y=filtfilt(ones(1,5)/5,1,y);
plot3(a66.data.frec(ilow:ihigh),y(ilow:ihigh),ones(numel(a66.data.frec(ilow:ihigh)),1),'-','linewidth',1.5,'color',[1 0.7 0])

% leo

ilow=find(a66.data.frec==280);
ihigh=find(a66.data.frec==370);

timesta=find(a66.data.startdate==datenum(2008,12,10));
timeend=find(a66.data.startdate==datenum(2009,1,1));

spldb=nanmean(a66.data.spec_dB_5min(ilow:ihigh,timesta:timeend));
SPL_scale = min(spldb):0.25:max(spldb);
[counts_per_bin,bin] = histc(spldb,SPL_scale);
cumulative_cpb = cumsum(counts_per_bin);
cumulative_cpb=cumulative_cpb/numel(spldb);
[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;

a=a66.data.spec_dB_5min(:,timesta:timeend);

[~,ix_min]=min(abs(cumulative_cpb-0.1));
ix_percentile=bin==ix_min;
ylow=nanmean(a(:,ix_percentile ),2);
ylow=filtfilt(ones(1,10)/10,1,ylow);

[~,ix_min]=min(abs(cumulative_cpb-0.9));
ix_percentile=bin==ix_min;
yhigh=nanmean(a(:,ix_percentile ),2);
yhigh=filtfilt(ones(1,10)/10,1,yhigh);

xxx=[a66.data.frec(ilow:ihigh) ;  flipud(a66.data.frec(ilow:ihigh)) ];
yyy=[ylow(ilow:ihigh) ;  flipud(yhigh(ilow:ihigh))];
ccc=[0.8 0 0.9];

fillhandle=patch(xxx,yyy,ccc);
set(fillhandle,'Facecolor',ccc,'EdgeColor','none');
h.leo=fillhandle

[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;
y=nanmean(a66.data.spec_dB_5min(:,a66.data.startdate>datenum(2008,12,10) & a66.data.startdate<datenum(2009,1,1)),2);
y=filtfilt(ones(1,5)/5,1,y);
plot3(a66.data.frec(ilow:ihigh),y(ilow:ihigh),ones(numel(a66.data.frec(ilow:ihigh)),1),'-','linewidth',1.5,'color',[0.5 0 0.5])

% crabeater

ilow=find(a66.data.frec==401);
ihigh=find(a66.data.frec==1000);

timesta=find(a66.data.startdate==datenum(2008,11,1));
timeend=find(a66.data.startdate==datenum(2008,11,31));

spldb=nanmean(a66.data.spec_dB_5min(ilow:ihigh,timesta:timeend));
SPL_scale = min(spldb):0.25:max(spldb);
[counts_per_bin,bin] = histc(spldb,SPL_scale);
cumulative_cpb = cumsum(counts_per_bin);
cumulative_cpb=cumulative_cpb/numel(spldb);
[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;

a=a66.data.spec_dB_5min(:,timesta:timeend);

[~,ix_min]=min(abs(cumulative_cpb-0.1));
ix_percentile=bin==ix_min;
ylow=nanmean(a(:,ix_percentile ),2);
ylow=filtfilt(ones(1,50)/50,1,ylow);

[~,ix_min]=min(abs(cumulative_cpb-0.9));
ix_percentile=bin==ix_min;
yhigh=nanmean(a(:,ix_percentile ),2);
yhigh=filtfilt(ones(1,50)/50,1,yhigh);

xxx=[a66.data.frec(ilow:ihigh) ;  flipud(a66.data.frec(ilow:ihigh)) ];
yyy=[ylow(ilow:ihigh) ;  flipud(yhigh(ilow:ihigh))];
ccc=[0 0.9 0.9];

fillhandle=patch(xxx,yyy,ccc);
set(fillhandle,'Facecolor',ccc,'EdgeColor','none');
h.crab=fillhandle

[~,ix_min]=min(abs(cumulative_cpb-0.5));
ix_percentile=bin==ix_min;
y=nanmean(a66.data.spec_dB_5min(:,a66.data.startdate>datenum(2008,11,1) & a66.data.startdate<datenum(2008,11,31)),2);
y=filtfilt(ones(1,100)/100,1,y);
plot3(a66.data.frec(ilow:ihigh),y(ilow:ihigh),ones(numel(a66.data.frec(ilow:ihigh)),1),'-','linewidth',1.5,'color',[0 0.5 0.8])

ax1=gca;
ax2 = axes('position',get(gca,'position'),'visible','off');

legend(ax1,[h.era_interim_noice,h.era_interim_ice],...
'Sea ice < 50% era_interim speed 0-10, 10-20, 20-30 m s^{-1}','Sea ice > 50% era_interim speed 0-10, 10-20 m s^{-1}')

legend(ax2,[h.blue,h.fin,h.minke,h.leo,h.crab],...
'Blue whales','Fin whales','Minke whales','Leopard seals','Crabeater seals')

%% plot correlation between sea ice and ambient sound in relation to frequency and spatial averaging radius
% Supplementary figure 2

frec_bins=logspace(log10(10),log10(16000),500);
frec_center= frec_bins(1:end-1) + diff(frec_bins)/2 ;

clear r

dv=datevec(icecon.time);
ix_time=dv(:,2)>1 & dv(:,2)<7;

for i_bin=2:numel(frec_bins)
    
    ix_frec=a66.data.frec>frec_bins(i_bin)-1 & a66.data.frec<frec_bins(i_bin);
    noise_time_series=nanmean(a66.data.spec_dB_5min(ix_frec,:),1);
    noise_time_series_int=interp1(a66.data.startdate,noise_time_series,icecon.time) ;
    ix_notnan=~isnan(noise_time_series_int) & ix_time';
    
    for i_radius=1:numel(radius_steps_km_ice)
        c=corrcoef(noise_time_series_int(ix_notnan),mean_ice(1,ix_notnan,i_radius));
        r_ice_a66(i_bin-1,i_radius)=c(1,2);
        clear c    
    end
    
    ix_frec=a69.data.frec>frec_bins(i_bin)-1 & a69.data.frec<frec_bins(i_bin);
    noise_time_series=nanmean(a69.data.spec_dB_5min(ix_frec,:),1);
    noise_time_series_int=interp1(a69.data.startdate,noise_time_series,icecon.time) ;
    ix_notnan=~isnan(noise_time_series_int)& ix_time';
    
    for i_radius=1:numel(radius_steps_km_ice)
        c=corrcoef(noise_time_series_int(ix_notnan),mean_ice(2,ix_notnan,i_radius));
        r_ice_a69(i_bin-1,i_radius)=c(1,2);
        clear c    
    end
    
    clear noise_time_series noise_time_series_int ix_frec ix_notnan 
     
end


frec=frec_center
nobio= frec>10 & frec<15 | frec>30 & frec<70 | frec>500 &  frec<10000 ;

figure(12)
set(gcf,'color','w')
clf

subplot(211)
hold on
contourf(frec_center(nobio),radius_steps_km_ice,r_ice_a66(nobio,:)',20,'edgecolor','none')

plot3([15,15],[60,1800],[10,10],'--w','linewidth',1.5)
plot3([15,30],[60,60],[10,10],'--w','linewidth',1.5)
plot3([15,30],[1800,1800],[10,10],'--w','linewidth',1.5)
plot3([30,30],[60,1800],[10,10],'--w','linewidth',1.5)

plot3([70,70],[60,1800],[10,10],'--w','linewidth',1.5)
plot3([70,500],[60,60],[10,10],'--w','linewidth',1.5)
plot3([70,500],[1800,1800],[10,10],'--w','linewidth',1.5)
plot3([500,500],[60,1800],[10,10],'--w','linewidth',1.5)

box on
cb=colorbar
set(gca,'xscale','log','yscale','log','xlim',[10 10000],'tickdir','out','clim',[-0.9 -0.4])
colormap('gray')
xlabel(cb,'r')
title('66°S')

subplot(212)
title('69°S')
hold on
contourf(frec_center(nobio),radius_steps_km_ice,r_ice_a69(nobio,:)',20,'edgecolor','none')

plot3([15,15],[60,1800],[10,10],'--w','linewidth',1.5)
plot3([15,30],[60,60],[10,10],'--w','linewidth',1.5)
plot3([15,30],[1800,1800],[10,10],'--w','linewidth',1.5)
plot3([30,30],[60,1800],[10,10],'--w','linewidth',1.5)

plot3([70,70],[60,1800],[10,10],'--w','linewidth',1.5)
plot3([70,500],[60,60],[10,10],'--w','linewidth',1.5)
plot3([70,500],[1800,1800],[10,10],'--w','linewidth',1.5)
plot3([500,500],[60,1800],[10,10],'--w','linewidth',1.5)

box on
cb=colorbar
set(gca,'xscale','log','yscale','log','xlim',[10 10000],'tickdir','out','clim',[-0.9 -0.4])
colormap('gray')

xlabel(cb,'r')

ylabel('Averaging radius in km')
xlabel('Frequeny in Hz')

%  set(gcf,'PaperPositionMode','auto')
%  print(gcf,'-dpdf','sea_ice_correlation','-r400') 
%   print(gcf,'-dtiff','sea_ice_correlation','-r400') 
