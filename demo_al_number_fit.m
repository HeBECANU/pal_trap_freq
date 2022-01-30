%demo_al_number_fit

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
anal_opts.tdc_import.shot_num=1:10;

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
data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,data);


%% FITTING THE ATOM NUMBER
%use the inital few atom laser pulses in order to determine the atom number
%not of that much benifit TBH
anal_opts.atom_num_fit=[];
anal_opts.atom_num_fit.pulses=[3,40]; %min,max index of pulses
anal_opts.atom_num_fit.plot.each_shot=true;
anal_opts.atom_num_fit.plot.history=false;
anal_opts.atom_num_fit.qe=anal_opts.global.qe;

data.num_fit=fit_pal_atom_number(anal_opts.atom_num_fit,data);


%% do it manualy and make a nice figure
anal_opts.atom_num_fit.pulses=[3,40]; %min,max index of pulses
%anal_opts.atom_num_fit.pulses=[3,160]; %min,max index of pulses


colors_main=[[98,136,63];[95,109,187];[180,72,117]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2)-0.2;
hsv(:,3)=hsv(:,3)+0.2;
colors_shaded=colorspace('HSV->RGB',hsv);



x_lim=[-1,1]*30e-6;
y_lim=x_lim;
xy_factor=1e6;

font_name='cmr10';
linewidth=1.5;
font_size=12;

stfig('Atom Number Fit')
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

ii=3;
pulse_num=anal_opts.atom_num_fit.pulses(1):anal_opts.atom_num_fit.pulses(2);
pulse_counts_val=data.mcp_tdc.al_pulses.num_counts(ii,pulse_num);
pulse_counts_err=sqrt(pulse_counts_val);
% add noise to estimate sensitivity
%pulse_counts_val=pulse_counts_val+normrnd(0,1,size(pulse_counts_val)).*pulse_counts_err; 

% fits work best when the scale of the fit parmeters are the same
% in this fit the atom number is much larger than the outcoupling fraction. So lets add a factor to keep them the same 
% order of mag.
fit_frac_factor=5e5/1e-2;
atom_fit_mdl = @(b,x) b(1)*(b(2)/fit_frac_factor)*((1-(b(2)/fit_frac_factor)).^(x(:,1)-1));
exp_qe=anal_opts.atom_num_fit.qe;

coefs_guess=[data.mcp_tdc.num_counts(ii),0.1*fit_frac_factor];
opts = statset('TolFun',1e-3,'TolX',1e-6,...
    'MaxIter',1e4,... %1e4
    'UseParallel',1);
fitobject=fitnlm(pulse_num,pulse_counts_val,atom_fit_mdl,coefs_guess,'Options',opts,...
                'ErrorModel','combined') %constant,proportional,combined

fit_out_frac_ntot=[[fitobject.Coefficients.Estimate(2)/fit_frac_factor,fitobject.Coefficients.SE(2)/fit_frac_factor];...
                  [fitobject.Coefficients.Estimate(1)/exp_qe,...
                  fitobject.Coefficients.SE(1)/exp_qe]];


 x_samp_mdl= linspace(min(pulse_num),max(pulse_num),1e3)';

[~,ci]=predict(fitobject,x_samp_mdl,'Alpha',1-erf(1/sqrt(2)),'Prediction','curve');
[prediction,oi]=predict(fitobject,x_samp_mdl,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');

hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2)-0.25;
hsv(:,3)=hsv(:,3)+0.65;
colors_shaded=colorspace('HSV->RGB',hsv);

hold on
ci_up=oi(:,1);
ci_down=oi(:,2);
oih=patch([x_samp_mdl', fliplr(x_samp_mdl')], ...
    [ci_up' , fliplr(ci_down')], colors_shaded(3,:),'EdgeColor','none',...
    'FaceAlpha',1);  %[1,1,1]*0.80

hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2)-0.2;
hsv(:,3)=hsv(:,3)+0.2;
colors_shaded=colorspace('HSV->RGB',hsv);

ci_up=ci(:,1);
ci_down=ci(:,2);
cih=patch([x_samp_mdl', fliplr(x_samp_mdl')], ...
    [ci_up' , fliplr(ci_down')], colors_shaded(3,:),'EdgeColor','none',...
    'FaceAlpha',1) ; %[1,1,1]*0.80


ebh=errorbar(pulse_num,pulse_counts_val,sqrt(pulse_counts_val),'o','LineWidth',linewidth,...
    'color',colors_main(1,:),'MarkerEdgeColor',colors_main(1,:), ...
    'MarkerFaceColor',colors_shaded(1,:),...
    'CapSize',0,'MarkerSize',4); %6
predh=plot(x_samp_mdl,prediction,'--','LineWidth',linewidth,'Color',colors_main(3,:));
hold off
xlabel('Atom Laser Pulse Index')
ylabel('Detected Atoms')
box on
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.025,0])
set(gca, {'XColor', 'YColor'}, {[1,1,1]*0.2, [1,1,1]*0.2});

legend([ebh,predh,cih,oih],'Meas.','Fit','CI','OI')

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=500;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off

fprintf('frac predicted %s\n',...
    string_value_with_unc(fit_out_frac_ntot(1,1),fit_out_frac_ntot(1,2),'type','b'))
fprintf('Total no. fit %se+05 measured in file %.2e (qe %.1e corrected)\n',...
    string_value_with_unc(fit_out_frac_ntot(2,1)*1e-5,fit_out_frac_ntot(2,2)*1e-5,'type','b')...
    ,data.mcp_tdc.num_counts(ii)/exp_qe,exp_qe)
pause(1e-9)

%%
%grid on
%fig_name='pal_atom_num_fit_short';
fig_name='pal_atom_num_fit_long';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))



% 


