%trap_curr_stab
%initalize
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));
hebec_constants
close('all')

%%
data_file='data/curr_stab_meas/imes0003.csv';
fp=fopen(data_file,'r');
%fgetl(fp);
raw_line1=fgetl(fp);
line1=split(raw_line1,',');
time_start=line1{7};
date_start=line1{6};
datetime_start=[date_start,' ',time_start];
datetime_start=datetime(datetime_start,'InputFormat', 'dd-MMM-yyyy HH:mm:ss',...
    'TimeZone','local','Format', 'yyyy-MM-dd HH:mm:ss.SSSxxxxx');
posix_start=posixtime(datetime_start);

while 1
  line_tmp = fgetl(fp);
  if ~ischar(line_tmp), break, end
  last_line=line_tmp;
end

fclose(fp);
line_last=split(last_line,',');
time_end=line_last{7};
date_end=line_last{6};
datetime_end=[date_end,' ',time_end];
datetime_end=datetime(datetime_end,'InputFormat', 'dd-MMM-yyyy HH:mm:ss',...
    'TimeZone','local','Format', 'yyyy-MM-dd HH:mm:ss.SSSxxxxx');
posix_end=posixtime(datetime_end);



%
%
array=readtable(data_file);
meas_v = array{:, 2};
meas_v=meas_v(1:end-1);
samp_rate=(numel(meas_v))/(posix_end-posix_start);
samp_time=(0:(numel(meas_v)-1))*(1/samp_rate);


stfig('raw_data')
plot(samp_time,meas_v)

resistance=0.00962;
meas_curr=meas_v/resistance;

%%

stfig('current')
plot(samp_time,meas_curr)


%% lets find the trap off times

font_name='cmr10';
linewidth=1.5;
font_size=12;
colors_main=[[[214,72,154];
[102,181,69];
[151,90,214];
[208,153,44];
[212,73,58]]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2)-0.3;
hsv(:,3)=hsv(:,3)+0.25;
hsv(1,2)=hsv(1,2)+0.1;
colors_shaded=colorspace('HSV->RGB',hsv);

%

trap_on_logical=meas_curr>10;
stfig('trap on logic');
plot(trap_on_logical)
trap_off_idxs=strfind(trap_on_logical',[1 0]);
trap_off_times=samp_time(trap_off_idxs);

start_samp=-28.75;
end_samp=-0.1;
filt_sigma=0.01;

shot_data=[];
shot_data.start_time=trap_off_times*nan;
shot_data.end_time=trap_off_times*nan;
shot_data.mean_currents=trap_off_times*nan;
shot_data.currents=cell(1,numel(trap_off_times));
shot_data.samp_times=cell(1,numel(trap_off_times));
shot_data.mean_current_const=cell(1,numel(trap_off_times));
for ii=1:numel(trap_off_times)
    mask_tmp=samp_time>trap_off_times(ii)+start_samp & samp_time<trap_off_times(ii)+end_samp;
    shot_data.currents{ii}=meas_curr(mask_tmp);
    shot_data.start_time(ii)=samp_time(find(mask_tmp,1,'first'));
    shot_data.end_time(ii)=samp_time(find(mask_tmp,1,'last'));
    shot_data.samp_times{ii}=(0:(numel(shot_data.currents{ii})-1))*(1/samp_rate);
    shot_data.mean_currents(ii)=mean(meas_curr(mask_tmp));
    shot_data.mean_current_const{ii}=shot_data.currents{ii}*0+shot_data.mean_currents(ii);
    shot_data.currents_smooth{ii}=gaussfilt(shot_data.samp_times{ii}, shot_data.currents{ii},filt_sigma);
end

cat_time_breaks=cellfun(@numel,shot_data.samp_times);
cat_time_breaks=cumsum(cat_time_breaks);

study_cat_currents=cat(1,shot_data.currents{:});
study_cat_smooth_currents=cat(1,shot_data.currents_smooth{:});
study_cat_mean_current=cat(1,shot_data.mean_current_const{:});
study_cat_times=(0:(numel(study_cat_currents)-1))*(1/samp_rate);
cat_time_breaks=study_cat_times(cat_time_breaks);
fh=stfig('cat curr')
clf
colororder({'k','k'})

frac_change=(study_cat_currents-mean(study_cat_currents))/mean(study_cat_currents);
frac_change_means_shots=(study_cat_mean_current-mean(study_cat_currents))/mean(study_cat_currents);
frac_change_means=(shot_data.mean_currents-mean(study_cat_currents))/mean(study_cat_currents);
frac_change_smooth_currents=(study_cat_smooth_currents-mean(study_cat_currents))/mean(study_cat_currents);
plot(study_cat_times,1e6*frac_change,'LineWidth',linewidth,'Color',colors_shaded(2,:))

hold on
plot(study_cat_times,1e6*frac_change_smooth_currents,'LineWidth',linewidth,'Color',colors_shaded(1,:))
plot(study_cat_times,1e6*frac_change_means_shots,':','LineWidth',linewidth*2,'Color',colors_main(3,:))
%xline(cat_time_breaks,'--','LineWidth',linewidth,'Color',colors_main(3,:))
xlim([study_cat_times(1),study_cat_times(end)]+[-0.1,0.1])
ylim([-110,130])
ylim
hold off
xlabel('Time (s)')
ylabel('Frac. Change in Current, $\Delta I/\langle I \rangle$ (ppm) ')
legend('Measurement','Gauss. Smooth (0.01 s)','Shot Mean')
%at %.2f Amps mean(study_cat_currents)
legend('FontSize',font_size)
legend('Location','North')
set(gca,'linewidth', 1.2)
set(gca,'TickLength',[0.02,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
yl=ylim;
yyaxis right
plot(1,1,'.w')
ylabel('Change in Current, $\Delta I$ (mA)')
ylim(yl.*mean(study_cat_currents)*1e-3)

% for gaussfilt sigma*BW(-3db) =  0.1876 see bw_of_gaussfilt
fprintf('filt BW 3db %f\n',0.1876/filt_sigma)

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=500;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[500,355,fig_width_px,fig_width_px*fig_aspect_ratio])

%%
fig_name='trap_current_single_shot';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))


%%
fh=stfig('current spectrum')
clf
colororder({'k','k'})
ftout=fft_tx(study_cat_times,frac_change,'padding',2,'window','chebyshev','win_param',{200})
plot(ftout(1,:),abs(ftout(2,:))*1e6,'LineWidth',linewidth,'Color',colors_main(2,:))
xlim([0,125])
xlabel('Frequency (Hz)')
ylabel('Frac. Current Amplitude (ppm)')
yl=ylim;
yyaxis right
plot(1,1,'.w')
ylabel('Current Amplitude (mA)')
ylim(yl.*mean(study_cat_currents)*1e-3)
legend('','')
legend('Box','off')
legend('Location','best')

%at %.2f Amps mean(study_cat_currents)
legend('FontSize',font_size)
legend('Location','North')
set(gca,'linewidth', 1.2)
set(gca,'TickLength',[0.02,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
grid on

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=500;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[500,355,fig_width_px,fig_width_px*fig_aspect_ratio])




%%
fig_name='trap_current_spectrum_lin';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))

%%
yyaxis left
set(gca,'YScale','log')
ylim([1,100])
yl=ylim;
yyaxis right
set(gca,'YScale','log')
ylim(yl.*mean(study_cat_currents)*1e-3)
pause(1)
%%
fig_name='trap_current_spectrum_log';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))

%%

fh=stfig('current')
clf


plot(samp_time,meas_curr,'LineWidth',linewidth,'Color',colors_main(2,:))
hold on
%yline(0)
ylim([0,max(meas_curr)]+[-0.5,0.1])
yl=ylim;

for ii=1:numel(shot_data.start_time)
    patch([shot_data.start_time(ii),shot_data.end_time(ii), shot_data.end_time(ii),shot_data.start_time(ii),], ...
    [yl(1),yl(1), yl(2),yl(2)],colors_shaded(3,:),'EdgeColor','none',...
    'FaceAlpha',0.3)
    %')  %[1,1,1]*0.80 , 
end

legend('Measurment','Subset')
legend('FontSize',font_size)
hold off
set(gca,'TickLength',[0.02,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
xlabel('Time (s)')
ylabel('Trap Current (Amps) ')
xlim([min(samp_time),max(samp_time)])

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=700;
fig_aspect_ratio=0.5; %0.67;
set(gcf,'Position',[500,355,fig_width_px,fig_width_px*fig_aspect_ratio])
%%
fig_name='trap_current_single';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))


