%pal_osc_bec_hist

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
anal_opts.atom_laser.t0=0.2641  ;%center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=130;%130;%50;%300;%375;%260;%300;%
anal_opts.atom_laser.appr_osc_freq_guess=[11,11,11];
anal_opts.atom_laser.pulse_twindow=[-0.5,0.5]*anal_opts.atom_laser.pulsedt*0.9;
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
anal_opts.osc_fit.dimesion=3; %Select coordinate to bin. 1=X, 2=Y.

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
%num_thresh=mean(data.mcp_tdc.num_counts(not_zero_files))-2*std(data.mcp_tdc.num_counts(not_zero_files));
num_thresh=1e4;
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
%anal_opts.osc_fit.dimesion=3; %Select coordinate to bin. 1=X, 2=Y.

anal_opts.osc_fit.global=anal_opts.global;
data.osc_fit=fit_trap_freq(anal_opts.osc_fit,data);
%%
data.osc_fit.trap_freq_recons=nan*data.osc_fit.ok.did_fits;
mask=data.osc_fit.ok.did_fits & (data.osc_fit.model_coefs(:,2,1)<52).' & (data.osc_fit.model_coefs(:,2,1)>20).' & ...
                                (abs(data.osc_fit.model_coefs(:,1,1))<50e-3).' ;
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


%%
mean_fit_params=squeeze(nanmean(data.osc_fit.model_coefs(mask,:,:),1))


% mean_fit_params =
% 
%     0.0111    0.0004
%    46.2844    0.0099
%     1.5425    0.0074
%     0.0001    0.0004
%    -0.1924    0.1812
%    -0.0279    0.0688
%     0.0414    0.0408
%    -0.0001    0.0005

hist(wrapTo2Pi(data.osc_fit.model_coefs(mask,3,1)*2*pi))
mean(wrapTo2Pi(data.osc_fit.model_coefs(mask,3,1)*2*pi))/pi

mean_pos=nanmean(squeeze(nanmean(data.mcp_tdc.al_pulses.vel_zxy.mean,1)),1);

%%

anal_opts.atom_laser.pulse_twindow=[-0.5,0.5]*anal_opts.atom_laser.pulsedt*2;
data.mcp_tdc.all_ok=data.mcp_tdc.num_ok;
data.mcp_tdc.al_pulses=bin_al_pulses(data.mcp_tdc,anal_opts.atom_laser);


%%


font_name='cmr10';
linewidth=1.5;
font_size=12;

fh=stfig('hists vel');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

iimax=4;
for ii=1:iimax
    fh_s=subplot(1,iimax,ii);

    zxy_counts=size(data.mcp_tdc.al_pulses.vel_zxy.counts);
    zxy_tmp=cat(1,data.mcp_tdc.al_pulses.vel_zxy.counts{data.mcp_tdc.num_ok,ii});
    zxy_tmp=zxy_tmp-repmat(mean_pos,[size(zxy_tmp,1),1]);
    xyz_tmp=zxy_tmp(:,[2,3,1]);


    % plots
    hist_opt=[];
    hist_opt.cords_to_plot=[2,3];
    hist_opt.cords_labels={'x','y','z'};
    hist_opt.cords_units=repmat({'m/s'},[1,3]);
    hist_opt.mean_marker=true;
    hist_opt.disp_scale=[1,1,1]*1e-3;
    hist_opt.filter_in_all_axis=true;

    hist_opt.kernel_factor=10;
    % opts.std_norm_bounds=true;
    % opts.hist_bounds=[5,1,1]*3;
    hist_opt.std_norm_bounds=false;
    %opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
    hist_opt.cmap=viridis(512);
    hist_opt.density_power=0.4;
    hist_opt.hist_bounds=[[-40,40];[-40,40];[-25,47]]*1.1*1e-3;
    hist_opt.blur_size=[1,1,1]*5e-4; %0.3e-3
    % optional
    hist_opt.filter_in_other_axis=false;
    hist_opt.norm_den_to_num_in=true;
    hist_opt.scale_den_disp_fac=true;
    hist_opt.norm_den_unity=true;
    %opts.cbar_lims=[0,2e-3];
    % opts.bin_factor;
    hist_opt.save_hist_as_image=false;
    hist_opt.open_fig=false;
    %hist_opt.save_hist_name=strcat(save_path,sprintf('raw_%0.4u.png',jj));
    
    
    hist_fig=hist_2d_data( xyz_tmp,hist_opt);

    if ii~=iimax
        ch=colorbar;
        delete(ch)
    end

    string_time=sprintf('t=%.1f ms',(ii-1)*anal_opts.atom_laser.pulsedt*1e3);
    sp_pos=get(gca,'Position');
    an=annotation('textbox', [sp_pos(1)+sp_pos(3)*0.2, sp_pos(2)+sp_pos(4)*0.95, 0.08, 0.001], 'string', string_time,...
        'FontSize',round(font_size*1.2),'FontName',font_name,'Color',[1,1,1]);
    an.FitBoxToText='on';
    an.LineStyle='none';
    set(gca,'linewidth', 1.1)
    set(gca,'TickLength',[0.04,0])
    set(gca, {'XColor', 'YColor'}, {[1,1,1]*0.2, [1,1,1]*0.2});

end



set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=2000;
fig_aspect_ratio=0.15; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])

%grid on
fig_name='pal_hist_vel_out_times_exp';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))
