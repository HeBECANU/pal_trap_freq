% monitor an experiment that uses measrument of the trap freq
% It would be usefull to get a decent idea of what the trap frequency is doing during a run without
% the need for a full processing of the data as in main
%  - processing each shot as it is made
%  - plot a history of the trap freq

% Known BUGS/ Possible Improvements
%
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-10-01
% BEGIN USER VAR-------------------------------------------------

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path
hebec_constants %call the constants function that makes some globals


%%

anal_opts=[];

anal_opts.tdc_import.dir = 'Z:\EXPERIMENT-DATA\2020_trap_freq_met\misc\aliasing_trap_freq\kicked_atom_laser_image';

%'Z:\EXPERIMENT-DATA\2020_trap_freq_met\stability_measures\20181017_two_probe_off_trap_freq_drift';

% 'Z:\EXPERIMENT-DATA\2020_trap_freq_met\20181017_probe_off_trap_freq_drift';
% 'Z:\EXPERIMENT-DATA\2020_trap_freq_met\aliasing_trap_freq\segmented_long\run1';
% 'Z:\EXPERIMENT-DATA\2020_trap_freq_met\20210616_trap_stability_measure';%20210622_good_trap_freq_mesaure';%20210622_low_fit_err_good_trap_freq_mesaure';%

%'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\'D:\trap_freq_data\normal_trap_freq\freq_space_method\20170410 80.00Hz';

%Z:\EXPERIMENT-DATA\2020_trap_freq_met\aliasing_trap_freq\aliasing\d 80 or
%160 or 90
%Z:\EXPERIMENT-DATA\2020_trap_freq_met\20210615_tune_trap_pal_good_atom_num_2\d
%Z:\EXPERIMENT-DATA\2020_trap_freq_met\aliasing_trap_freq\segmented_long\d

anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=false;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];


anal_opts.atom_laser.pulsedt=7.5e-3;%10.000e-3;%8.000e-3;%
anal_opts.atom_laser.t0=0.264;%0.41784; %0.8762;%center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=190;%130;%50;%300;%375;%260;%300;%
anal_opts.atom_laser.appr_osc_freq_guess=[11,11,11];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;


anal_opts.osc_fit.binsx=1000;
anal_opts.osc_fit.blur=0.3;
anal_opts.osc_fit.xlim=[-20,20]*1e-3;
anal_opts.osc_fit.tlim=[0.86,1.08];
anal_opts.osc_fit.dimesion=2; %Select coordinate to bin. 1=X, 2=Y.

anal_opts.history.shots=200;

% END USER VAR-----------------------------------------------------------

%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
% addpath(genpath(folder));

hebec_constants
anal_opts.tdc_import.mat_save=false;
anal_opts.global.velocity=const.g0*anal_opts.global.fall_time;
anal_opts.global.fall_velocity=const.g0*anal_opts.global.fall_time; %velocity when the atoms hit the detector
% fall_dist=1/2 a t^2 
%TODO get from engineering documents
anal_opts.global.fall_dist=(1/2)*const.g0*anal_opts.global.fall_time^2;
anal_opts.atom_laser.plot=false;
%anal_opts.atom_laser.t0=0.417770;
anal_opts.atom_laser.global=anal_opts.global; %coppy global into the options structure
const.fall_disntace = anal_opts.global.fall_dist;
anal_opts.atom_laser.c = const;

% if anal_opts.tdc_import.dir(end) ~= '\', dirpath = [dirpath '\']; end
% if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end
%  
% anal_out.dir=sprintf('%sout\\monitor\\',...
%     anal_opts.tdc_import.dir);
% if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
% anal_opts.global.out_dir=anal_out.dir;




%%
loop_num=0;
sfigure(1); 
set(gcf,'color','w')



anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
%anal_opts.tdc_import.shot_num=47:65;

%max_shot_num=max(anal_opts.tdc_import.shot_num);
%anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
%    anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));

data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);


%% CHECK ATOM NUMBER
sfigure(1)
%create a list of indicies (of the mcp_tdc) that have an ok number of counts
%exclude the very low and then set the thresh based on the sd of the remaining
not_zero_files=data.mcp_tdc.num_counts>1e3; 
num_thresh=mean(data.mcp_tdc.num_counts(not_zero_files))-4*std(data.mcp_tdc.num_counts(not_zero_files));
data.mcp_tdc.num_ok=data.mcp_tdc.num_counts>num_thresh;
 %   (data.mcp_tdc.time_create_write(:,1)'-data.mcp_tdc.time_create_write(1,1))<(anal_opts.max_runtime*60*60);
fprintf('shots number ok %u out of %u \n',sum(data.mcp_tdc.num_ok),numel(data.mcp_tdc.num_ok))

plot((data.mcp_tdc.time_create_write(:,2)-data.mcp_tdc.time_create_write(1,2))/(60*60),data.mcp_tdc.num_counts)
xlabel('time (h)')
ylabel('total counts')
title('num count run trend')


%% Bin pulses

data.mcp_tdc.all_ok=data.mcp_tdc.num_ok;
data.mcp_tdc.al_pulses=bin_al_pulses(data.mcp_tdc,anal_opts.atom_laser);

%%
anal_opts.osc_fit.adaptive_freq=true; %estimate the starting trap freq 
anal_opts.osc_fit.appr_osc_freq_guess=[52,46.7,40];%[50,12,37];%
anal_opts.osc_fit.plot_fits=false;
anal_opts.osc_fit.plot_err_history=true;
anal_opts.osc_fit.plot_fit_corr=false;
anal_opts.osc_fit.freq_fit_tolerance = 1;

anal_opts.osc_fit.global=anal_opts.global;
data.osc_fit=fit_trap_freq(anal_opts.osc_fit,data);
%%
data.osc_fit.trap_freq_recons=nan*data.osc_fit.ok.did_fits;
mask=data.osc_fit.ok.did_fits & (data.osc_fit.model_coefs(:,2,1)<49).' & (data.osc_fit.model_coefs(:,2,1)>20).';
data.osc_fit.trap_freq_recons(mask)=3*(1/anal_opts.atom_laser.pulsedt)+data.osc_fit.model_coefs(mask,2,1);
data.osc_fit.trap_freq_unc=nan*mask;
data.osc_fit.trap_freq_unc(mask)=data.osc_fit.model_coefs(mask,2,2)';

%%
n=1;
%[800 2375];% or [450 2378];%
shot_lims = [450 2378];%[800 2375];%[800 2375];%[450 2378];%[155 2375];%[0 inf];
mask2 = mask & (data.mcp_tdc.shot_num>shot_lims(1)) & (data.mcp_tdc.shot_num<shot_lims(2));
std_2 = nanstd(data.osc_fit.trap_freq_recons(mask2));
mean_2 = nanmean(data.osc_fit.trap_freq_recons(mask2));
%nanmean(data.osc_fit.trap_freq_unc(mask2));

sfigure(10099+n);
clf
errorbar(data.osc_fit.dld_shot_num(mask2),...
    data.osc_fit.trap_freq_recons(mask2),data.osc_fit.trap_freq_unc(mask2),...
    'kx-','MarkerSize',7,'CapSize',0,'LineWidth',1)
hold on
t_vec = linspace(data.osc_fit.dld_shot_num(1).*0.8,data.osc_fit.dld_shot_num(end).*1.2,1e4);
mean_vec = mean_2.*ones(size(t_vec));
plot(t_vec,mean_vec+std_2,'r-','LineWidth',1.5)
plot(t_vec,mean_vec-std_2,'r-','LineWidth',1.5)
plot(t_vec,mean_vec,'b--','LineWidth',1.5)
xlabel('Shot Number')
ylabel('Fit Trap Freq')
pause(1e-6)


data_time=data.mcp_tdc.time_create_write(data.osc_fit.dld_shot_idx,2);
data_time=(data_time-min(data_time))./3600;
h1=stfig('drift');
clf
time_vec = data_time(mask2);
% errorbar(time_vec-time_vec(1),...
%     data.osc_fit.trap_freq_recons(mask2)-mean_2,data.osc_fit.trap_freq_unc(mask2),...
%     'ko','MarkerSize',1.1,'CapSize',0,'LineWidth',1,'MarkerFaceColor','k')
plot(time_vec-time_vec(1),...
    1e3.*(data.osc_fit.trap_freq_recons(mask2)-mean_2),...
    '-','MarkerSize',1.1,'LineWidth',0.41,'Color',[1,1,1]*0.7)
hold on
delta_freq_smooth=gaussfilt(time_vec-time_vec(1),1e3.*(data.osc_fit.trap_freq_recons(mask2)'-mean_2),5/60)
plot(time_vec-time_vec(1),...
    delta_freq_smooth,...
    '-','MarkerSize',1.1,'LineWidth',0.41,'Color',[1,1,1]*0)
t_vec = linspace(data_time(1).*0.8-data_time(1),data_time(end).*1.2,1e4);
mean_vec = 0.*ones(size(t_vec));
plot(t_vec,1e3.*mean_vec,'b--','LineWidth',1.5)
plot(t_vec,1e3.*(mean_vec+std_2),'r-.','LineWidth',1.5)
plot(t_vec,1e3.*(mean_vec-std_2),'r-.','LineWidth',1.5)
xlabel('Time (hours)')
ylabel('Trap Drift, $f-\bar{f}$ (mHz)')
xlim([0 max(time_vec-time_vec(1))])
set(gca,'FontSize',16)
pause(1e-6)
legend('measurement','gauss. filt ($\tau=300$ s)','$\bar{f}$','$\bar{f}\pm\sigma_f$',...
    'Location','southwest')
legend('FontSize',12,'NumColumns',2)
legend('boxoff')

nanstd(data.osc_fit.trap_freq_recons(mask2));
nanmean(data.osc_fit.trap_freq_unc(mask2));
set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h1,'trap_drift_vs_time','-dpdf','-r0')
set(h1,'Position',[2,2,5*1.3,4*1.3])
export_fig('./figs/trap_drift_vs_time.pdf')
%export_fig('./figs/trap_drift_vs_time.svg')


%% histogram
h2=stfig('freq hist')
clf
shot_lims = [450 2378];%[18 1500];
mask3 = mask & (data.mcp_tdc.shot_num>shot_lims(1)) & (data.mcp_tdc.shot_num<shot_lims(2));
mean_3 = nanmean(data.osc_fit.model_coefs(mask3,2,1));

h=histogram(1e3.*(data.osc_fit.model_coefs(mask3,2,1)-mean_3),100);
pd = fitdist(-mean_3+data.osc_fit.model_coefs(mask3,2,1),'Normal');
x_values = linspace(46-mean_3,46.5-mean_3,1e5);
y = pdf(pd,x_values);
hold on
plot(x_values.*1e3,y.*2,'LineWidth',1.8)
xlim(1e3.*([46.23 46.37]-mean_3))
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [0.5 0.5 0.5];
set(gca,'FontSize',14)
xlabel('Trap Drift, \(f-\bar{f}\) (mHz)')
ylabel('Counts')
legend('measurements','gaussian fit')
legend('boxoff')

set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h2,'trap_drift_hist','-dpdf','-r0')
set(h2,'Position',[2,2,5*1.3,4*1.3])
export_fig('./figs/trap_drift_hist.svg')
exportgraphics(h2,'./figs/trap_drift_hist.pdf')


h3=stfig('single shot uncert hist')
clf
shot_lims = [450 2378];%[18 1500];
mask3 = mask & (data.mcp_tdc.shot_num>shot_lims(1)) & (data.mcp_tdc.shot_num<shot_lims(2));
h2=histogram(data.osc_fit.model_coefs(mask3,2,2).*1e3,100);

pd2 = fitdist(data.osc_fit.model_coefs(mask3,2,2).*1e3,'Normal');
x_values = linspace(4,9,1e5);
y = pdf(pd2,x_values);
hold on
plot(x_values,y.*90,'LineWidth',1.8)
xlim([4.3,9])
xlabel('Single Shot Uncertainty, \(\delta f\) (mHz)')
ylabel('Counts')
set(gca,'FontSize',14)
h2.FaceColor = [0.5 0.5 0.5];
h2.EdgeColor = [0.5 0.5 0.5];

set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h3,'unc_hist','-dpdf','-r0')

%% allan

shot_lims = [450 2378];%[800 2375];%[155 500];%[155 800];%[155 2375];%[18 802];
allan_data=[];
allan_data.freq=data.osc_fit.trap_freq_recons-mean(data.osc_fit.trap_freq_recons);
allan_data.time=data.mcp_tdc.time_create_write(data.osc_fit.dld_shot_idx,2);
mask=~isnan(allan_data.freq) & ~isnan(allan_data.time)' & (data.mcp_tdc.shot_num>shot_lims(1)) & (data.mcp_tdc.shot_num<shot_lims(2));
allan_data.freq=allan_data.freq(mask);
allan_data.time=allan_data.time(mask);
allan_data.time=allan_data.time-min(allan_data.time);
tau_in=logspace(log10((mean(diff(allan_data.time)))), ...
             log10(allan_data.time(end)/3),2e2);%allan_data.time;%linspace(0.1,10^4);%


[sigma_tau, s, errorb, tau]=allan(allan_data,tau_in);

%%
fit_options=[];
fit_options.windows=[[30 3000];[3000,20000]];
fit_options.plot_options={{'--','Color','r'},{':','Color','b'}}


[ha,fit_dets]=plot_allan_with_fits(sigma_tau, errorb, tau,fit_options,'Allan Deviation');
set(gca,'FontSize',16)


%%
fit_options=[];
fit_options.windows=[[30 3000];[3000,20000]];
fit_options.plot_options={{'--','Color','r'},{':','Color','b'}}


[ha,fit_dets]=plot_allan_with_fits(sigma_tau, errorb, tau,fit_options,'Normal Allan Deviation');
set(gca,'FontSize',16)

%%

export_fig('./figs/bthesis/trap_allan.svg')
exportgraphics(ha,'./figs/bthesis/trap_allan.pdf')


%%
set(ha,'Units','Inches');
pos = get(ha,'Position');
set(ha,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ha,'allan_deviation','-dpdf','-r0')

%%
% [retval, s, errorb, tau]=allan_overlap(allan_data,2.^(-10:0.005:20))

%%
fit_options=[];
fit_options.windows=[[30 3000]];
fit_options.plot_options={{'--','Color','r'}}


[ha,fit_dets]=plot_allan_with_fits(sigma_tau, errorb, tau,fit_options,'Allan Deviation');
set(gca,'FontSize',16)
set(ha,'Units','Inches');
pos = get(ha,'Position');
set(ha,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%%
set(gcf,'Position',[2,2,6*1.3,4*1.3])
export_fig('./figs/trap_allan.svg')
exportgraphics(gca,'./figs/trap_allan.pdf')


%%
function [h,fits]=plot_allan_with_fits(sigma_tau, errorb, tau,fit_options,dev_label)
plot_handles={};
windows=fit_options.windows;
font_size = 18;
font_type='cmr10' 
mask_outliers=sigma_tau<1;
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
labels={dev_label,'Uncertainty'};
for ii=1:size(windows,1)
    windows(ii,:)
    mask_fit=tau>windows(ii,1) & tau<windows(ii,2) & mask_outliers;
    fits{ii}=polyfit(log10(tau(mask_fit)),log10(sigma_tau(mask_fit)),1);
    fits{ii}
    mask_plot= tausamp>windows(ii,1) & tausamp<windows(ii,2);
    plot_handles{2+ii}=plot(tausamp(mask_plot),10.^polyval(fits{ii},log10(tausamp(mask_plot))),...
        fit_options.plot_options{ii}{:},'LineWidth',2);
    %labels{end+1}=sprintf('Lin. Fit $\\sigma(\\tau)=%.2f\\tau^{%.2f}$',-fits{ii}(2),fits{ii}(1))
    labels{end+1}=sprintf('Lin. Fit')
end
xlim(xl)
ylim(yl)
hold off
%title(dev_label,'FontSize',font_size*1.5,'FontName',font_type)
xlabel('Averaging, $\tau$ (s)','FontSize',font_size,'FontName',font_type)
ylabel('Allan Deviation, $\sigma(\tau)$ (Hz)','FontSize',font_size,'FontName',font_type)
ln=legend([plot_handles{2},plot_handles{1},plot_handles{3:end}],labels,'FontName',font_type,'FontSize',font_size*0.8,'location','southwest')
ln.EdgeColor=[1,1,1];
grid on
box on


end


