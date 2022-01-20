%clock_freq_stab_images

%vid_file='F:\WIN_20220117_17_46_12_Pro.mp4'
img_dir='F:\10s_gate\';
flist=dir(img_dir);
fnames={flist.name};
fnames(1)=[];
fnames(1)=[];

%% process names into posix time
iimax=numel(fnames);
times_vec=nan(iimax,1);
for ii=1:iimax
    str_tmp=fnames{ii};
    str_tmp=strsplit(str_tmp,'_');
    str_tmp=str_tmp{end};
    str_tmp=erase(str_tmp,'.jpg');
    str_tmp=strrep(str_tmp,'T',' ');
    dt_tmp=datetime(str_tmp,'InputFormat','yyyyMMdd HHmmss.SSS');
    time_posix_tmp=posixtime(dt_tmp);
    times_vec(ii)=time_posix_tmp;
end


%%
test_img=imread(fullfile(img_dir,fnames{3}));
rotation_angle=2;
test_img=imrotate(test_img,rotation_angle, 'bicubic');
figure(2)
imshow(test_img)
%roi = drawrectangle('Color','g');
% roi_pos=round(roi.Position);
roi_pos=[567   528   655    87];
im_roi=test_img(roi_pos(2):(roi_pos(2)+roi_pos(4)),roi_pos(1):(roi_pos(1)+roi_pos(3)),:);
imshow(im_roi)


%%

iimax=numel(fnames);
poll_frames_img=cell(iimax,1);
fprintf('getting frames')
for ii=1:iimax
    fprintf('%03uof%03u\n',ii,iimax)
    im_tmp=imread(fullfile(img_dir,fnames{ii}));;
    im_tmp=imrotate(im_tmp,rotation_angle, 'bicubic');
    im_tmp=im_tmp(roi_pos(2):(roi_pos(2)+roi_pos(4)),roi_pos(1):(roi_pos(1)+roi_pos(3)),:);
    poll_frames_img{ii}=im_tmp;
    %imwrite(im_tmp,sprintf('F:\\subregions\\%u.jpg',ii))
end

%% run ocr
trainedLanguage = 'E:\Dropbox\UNI\programs\Trap_freq_methods\lib\freq_count_ocr\sevensegment\tessdata\sevensegment.traineddata';

iimax=numel(poll_frames_img);
frames_text=cell(iimax,1);
frames_ocr_obj=cell(iimax,1);
fprintf('running ocr')
for ii=1:iimax
    fprintf('%03uof%03u\n',ii,iimax)
    im_ocr_tmp=poll_frames_img{ii};
    %im_roi_binarized=im_roi(:,:,2)>190;
    ocrResults = ocr(im_ocr_tmp,'CharacterSet','0123456789.','Language', trainedLanguage, ...
            'TextLayout', 'Block');
    frames_ocr_obj{ii}=ocrResults;
    frames_text{ii}=strtrim(ocrResults.Text);
end

%% process into a vector of values
freq_meas=[];
freq_meas.times=times_vec;
freq_meas.freqs=nan(iimax,1);
iimax=numel(freq_meas.times);
fprintf('converting to vector')
for ii=1:iimax
    str_tmp=frames_text{ii};
    str_tmp=strrep(str_tmp,' ','')
    if str_tmp(end)~='3'
        warning('not ending with 3')
        pause
    end
    if str_tmp(3)~='.'
        warning('does not have a . inserting')
        str_tmp=insertAfter(str_tmp,2,'.');
    end
    %imshow(poll_frames_img{ii})
    if numel(str_tmp)~=11
        freq_meas.freqs(ii)=nan;
    else
        value=str2num(strtrim(str_tmp(1:end-1)))*1e3;
        freq_meas.freqs(ii)=value;
    end
    
end

%%
window_size_sec=60;
window_size_idx=round(window_size_sec/mean(diff(freq_meas.times)));
freq_outlier=isoutlier(freq_meas.freqs,'movmedian',window_size_idx,'ThresholdFactor',4);
freq_meas.freqs(freq_outlier)=nan;

freq_meas.freqs=fillmissing(freq_meas.freqs,'pchip');

%%
stfig('freq time series')
clf

font_size = 15;
font_type='cmr10' 
yscale_fac=1e6;
xscale=1e-3;

freq_mean=mean(freq_meas.freqs);
freq_std=std(freq_meas.freqs);
freq_frac_dev=(freq_meas.freqs-freq_mean)/freq_mean;
freq_frac_dev_smooth=gaussfilt(freq_meas.times-freq_meas.times(1),freq_frac_dev,100)

plot((freq_meas.times-freq_meas.times(1))*xscale,freq_frac_dev*yscale_fac,'Color',[1,1,1]*0.7)
hold on
plot((freq_meas.times-freq_meas.times(1))*xscale,freq_frac_dev_smooth*yscale_fac,'b')
hold off
legend('raw','filtered $\tau$=100 s','Location','north')
legend('boxoff')

xlabel('time (ks)','FontSize',font_size,'FontName',font_type)
ylabel('Fractional Deviation (ppm)','FontSize',font_size,'FontName',font_type)

size_mult=100
set(gcf,'Position',[100,100,6*size_mult,4*size_mult])
export_fig('./figs/bthesis/clock_drift_time_series.svg')
exportgraphics(gca,'./figs//bthesis/clock_drift_time_series.pdf')


%%


allan_data=[];
allan_data.freq=freq_meas.freqs;
allan_data.time=freq_meas.times-freq_meas.times(1);
%tau_in=2.^(-12:0.1:14);%allan_data.time;%linspace(0.1,10^4);%
tau_in=logspace(log10(10),log10(allan_data.time(end)/4),2e2)
%tau_in=tau_in(tau_in<4000)

[sigma_tau, s, sigma_tau_uncert, tau]=allan(allan_data,tau_in);




%%
yscale_fac=1e9;
font_size = 15;
font_type='cmr10' 

sigma_tau_scaled=sigma_tau*yscale_fac/mean(freq_meas.freqs);
tau=tau;
sigma_tau_uncert_scaled=sigma_tau_uncert*yscale_fac/mean(freq_meas.freqs);

stfig('clock error')
%errorbar(tau,retval,errorb)


%plot(tau,sigma_tau_scaled-sigma_tau_uncert_scaled,'Color',[1,1,1]*0.7)
%plot(tau,sigma_tau_scaled+sigma_tau_uncert_scaled,'Color',[1,1,1]*0.7)
fill([tau, fliplr(tau)],...
    [sigma_tau_scaled-sigma_tau_uncert_scaled, fliplr(sigma_tau_scaled+sigma_tau_uncert_scaled)],...
    [1,1,1]*0.7,'EdgeAlpha',0)
hold on
plot(tau,sigma_tau_scaled,'b')
hold off
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on
xlabel('Averaging, $\tau$ (s)','FontSize',font_size,'FontName',font_type)
ylabel('Fractional Allan Deviation, $\sigma_{f}(\tau)/\bar{f}$ (ppb)','FontSize',font_size,'FontName',font_type)
grid on
size_mult=100
set(gcf,'Position',[100,100,6*size_mult,4*size_mult])
export_fig('./figs/bthesis/clock_allan_dev.svg')
export_fig(gca,'./figs/bthesis/clock_allan_dev.pdf')



%%

fit_options=[];
fit_options.windows=[[10 80];[200,2000];[3000,20000]];
fit_options.plot_options={{'--','Color','r'},{':','Color','b'},{'-.','Color','g'}}
fit_options.legend_location='best'


[ha,fit_dets]=plot_allan_with_fits(sigma_tau_scaled, sigma_tau_uncert_scaled, tau,fit_options,'Allan Deviation');
set(gca,'FontSize',16)

size_mult=100
set(gcf,'Position',[100,100,6*size_mult,4*size_mult])
export_fig('./figs/bthesis/clock_allan_dev_fits.svg')
exportgraphics(gca,'./figs/bthesis/clock_allan_dev_fits.pdf')


%%
function [h,fits]=plot_allan_with_fits(sigma_tau, errorb, tau,fit_options,dev_label)
plot_handles={};
windows=fit_options.windows;
font_size = 18;
font_type='cmr10' 
mask_outliers=~isinf(sigma_tau);
h=stfig('allan plot');
clf
set(gcf,'color','w')
hold on
plot_handles{1}=fill([tau, fliplr(tau)],...
    [sigma_tau-errorb, fliplr(sigma_tau+errorb)],...
    [1,1,1]*0.9,'EdgeAlpha',0)
plot_handles{2}=plot(tau(mask_outliers),sigma_tau(mask_outliers),'k-','LineWidth',1.6)
% plot(tau(mask_outliers),retval(mask_outliers)+errorb(mask_outliers),'-','Color',[1,1,1]*0.5)
% plot(tau(mask_outliers),retval(mask_outliers)-errorb(mask_outliers),'-','Color',[1,1,1]*0.5)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
pause(1e-6)
xl=xlim;
yl=ylim;
fits=[];
tausamp=logspace(log10(xl(1)),log10(xl(2)),1e3);
labels={dev_label,'uncertainty'};
for ii=1:size(windows,1)
    windows(ii,:)
    mask_fit=tau>windows(ii,1) & tau<windows(ii,2) & mask_outliers;
    fits{ii}=polyfit(log10(tau(mask_fit)),log10(sigma_tau(mask_fit)),1);
    fits{ii}
    mask_plot= tausamp>windows(ii,1) & tausamp<windows(ii,2);
    plot_handles{2+ii}=plot(tausamp(mask_plot),10.^polyval(fits{ii},log10(tausamp(mask_plot))),...
        fit_options.plot_options{ii}{:},'LineWidth',2);
    labels{end+1}=sprintf('Lin. Fit $\\sigma(\\tau)=%.2f\\tau^{%.2f}$',-fits{ii}(2),fits{ii}(1))
end
xlim(xl)
ylim(yl)
hold off
%title(dev_label,'FontSize',font_size*1.5,'FontName',font_type)
xlabel('Averaging, $\tau$ (s)','FontSize',font_size,'FontName',font_type)
ylabel('Allan Deviation, $\sigma(\tau)$ (Hz)','FontSize',font_size,'FontName',font_type)
ln=legend([plot_handles{2},plot_handles{1},plot_handles{3:end}],labels, ...
    'FontName',font_type,'FontSize',font_size*0.8, ...
    'location',fit_options.legend_location)
%ln.EdgeColor=[1,1,1];
legend('boxoff')
grid on


end

