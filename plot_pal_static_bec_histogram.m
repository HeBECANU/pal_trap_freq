%pal_osc_bec_hist
% plot the 2d histogram for a static pulsed atom laser;

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path
hebec_constants %call the constants function that makes some globals


%%

anal_opts=[];

%anal_opts.tdc_import.dir = 'Z:\EXPERIMENT-DATA\2020_trap_freq_met\misc\aliasing_trap_freq\kicked_atom_laser_image';
anal_opts.tdc_import.dir = 'F:\scratch\Run 2 background no raman'

anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=false;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];


anal_opts.atom_laser.pulsedt=16e-3;%10.000e-3;%8.000e-3;%
anal_opts.atom_laser.t0=1.018  ;%center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=10;%130;%50;%300;%375;%260;%300;%
anal_opts.atom_laser.pulse_twindow=[-0.5,0.5]*anal_opts.atom_laser.pulsedt*1;
anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;




% END USER VAR-----------------------------------------------------------
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



%anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
anal_opts.tdc_import.shot_num=1:4807;

%max_shot_num=max(anal_opts.tdc_import.shot_num);
%anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
%    anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));

data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);


%% CHECK ATOM NUMBER
sfigure(1)
%create a list of indicies (of the mcp_tdc) that have an ok number of counts
%exclude the very low and then set the thresh based on the sd of the remaining
not_zero_files=data.mcp_tdc.num_counts>1e3; 
num_thresh=mean(data.mcp_tdc.num_counts(not_zero_files))-1*std(data.mcp_tdc.num_counts(not_zero_files));
%num_thresh=1e4;
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

pulse_vel_mean_windowed=squeeze(nanmean(data.mcp_tdc.al_pulses.vel_zxy.mean,1));

%%

anal_opts.atom_laser.pulse_twindow=[-0.5,0.5]*anal_opts.atom_laser.pulsedt*2;
data.mcp_tdc.all_ok=data.mcp_tdc.num_ok;
data.mcp_tdc.al_pulses=bin_al_pulses(data.mcp_tdc,anal_opts.atom_laser);


%%

tf_param=[];
tf_param.trap_freq=2*pi*[51,550,550]; 
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;
tf_param.num=7e5;
tf_details=bec_properties(tf_param.trap_freq/(2*pi),tf_param.num,tf_param.mass,tf_param.a_scat_len);
%%

font_name='cmr10';
linewidth=1.5;
font_size=12;

fh=stfig('hists vel');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

% plot pulse number
ii=1;

%mean_vel=nanmean(squeeze(nanmean(data.mcp_tdc.al_pulses.vel_zxy.mean,1)),1);

zxy_counts=size(data.mcp_tdc.al_pulses.vel_zxy.counts);
zxy_tmp=cat(1,data.mcp_tdc.al_pulses.vel_zxy.counts{data.mcp_tdc.num_ok,ii});
zxy_tmp=zxy_tmp -repmat(pulse_vel_mean_windowed(ii,:),[size(zxy_tmp,1),1]);
xyz_tmp=zxy_tmp(:,[2,3,1]);


% plots
hist_opt=[];
hist_opt.cords_to_plot=[2,3];
hist_opt.cords_labels={'x','y','z'};
hist_opt.cords_units=repmat({'m/s'},[1,3]);
hist_opt.mean_marker=true;
hist_opt.disp_scale=[1,1,1]*1e-3;
hist_opt.filter_in_all_axis=true;

hist_opt.kernel_factor=4;
% opts.std_norm_bounds=true;
% opts.hist_bounds=[5,1,1]*3;
hist_opt.std_norm_bounds=false;
%opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
hist_opt.cmap=viridis(512);
hist_opt.density_power=0.4;
hist_opt.hist_bounds=[[-45,45];[-45,45];[-40,100]]*1.1*1e-3;
hist_opt.blur_size=[1,1,1]*2.5e-4; %0.3e-3
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


[hist_fig,ch]=hist_2d_data( xyz_tmp,hist_opt);
ch.Ticks=linspace(0,ch.Ticks(end),6);

hold on
cir_r=tf_details.pal_vmax;
cir_theta = linspace(0.5*pi,1.5*pi,100);
cir_cen=[0,4]*1e-3;
plot((cir_cen(1) + cir_r*cos(cir_theta))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
    (cir_cen(2) + cir_r*sin(cir_theta))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),'--','Color',[209,75,175]./255,'LineWidth',2)
hold off

% if ii~=iimax
%     ch=colorbar;
%     delete(ch)
% end

%anotate
% string_time=sprintf('t=%.1f ms',(ii-1)*anal_opts.atom_laser.pulsedt*1e3);
% sp_pos=get(gca,'Position');
% an=annotation('textbox', [sp_pos(1)+sp_pos(3)*0.2, sp_pos(2)+sp_pos(4)*0.95, 0.08, 0.001], 'string', string_time,...
%     'FontSize',round(font_size*1.2),'FontName',font_name,'Color',[1,1,1]);
% an.FitBoxToText='on';
% an.LineStyle='none';

set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.04,0])
set(gca, {'XColor', 'YColor'}, {[1,1,1]*0.2, [1,1,1]*0.2});




%
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=500;
fig_aspect_ratio=0.8; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])

%%
fig_name='pal_hist_static_bcr_exp';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))
