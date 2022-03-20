%% trap freq plot

anal_opts=[];

anal_opts.tdc_import.dir = 'Z:\EXPERIMENT-DATA\2020_trap_freq_met\misc\aliasing_trap_freq\aliasing\';
shot=75;%33;%85;%6;%
%aliased signals
%Z:\EXPERIMENT-DATA\2020_trap_freq_met\aliasing_trap_freq\aliasing\d 75 (9.500e-3,0.262)
%80 (9.500e-3,0.262) or
%160 (13.000e-3,0.262) or 90 (10.000e-3,0.262)
%Z:\EXPERIMENT-DATA\2020_trap_freq_met\20210615_tune_trap_pal_good_atom_num_2\d
%(18.000e-3,0.41784)

%unaliased signal
%Z:\EXPERIMENT-DATA\2020_trap_freq_met\aliasing_trap_freq\segmented_long\run1\d
%(1.250e-3,0.87606)
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=false;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];


anal_opts.atom_laser.pulsedt=9.500e-3;%8.45e-3;%10e-3;%4e-3;%13.000e-3;%10.000e-3;%8.000e-3;%
anal_opts.atom_laser.t0=0.262;%0.87606;%0.41784; %0.262;%0.8762;%center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=156;%174;%375;%260;%300;%
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

anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
anal_opts.tdc_import.shot_num=shot;

%max_shot_num=max(anal_opts.tdc_import.shot_num);
%anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
%    anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));

data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);


%% CHECK ATOM NUMBER
% sfigure(1)
%create a list of indicies (of the mcp_tdc) that have an ok number of counts
%exclude the very low and then set the thresh based on the sd of the remaining
not_zero_files=data.mcp_tdc.num_counts>1e3;
num_thresh=mean(data.mcp_tdc.num_counts(not_zero_files))-4*std(data.mcp_tdc.num_counts(not_zero_files));
data.mcp_tdc.num_ok=not_zero_files;%data.mcp_tdc.num_counts>num_thresh;
%   (data.mcp_tdc.time_create_write(:,1)'-data.mcp_tdc.time_create_write(1,1))<(anal_opts.max_runtime*60*60);
fprintf('shots number ok %u out of %u \n',sum(data.mcp_tdc.num_ok),numel(data.mcp_tdc.num_ok))

% plot((data.mcp_tdc.time_create_write(:,2)-data.mcp_tdc.time_create_write(1,2))/(60*60),data.mcp_tdc.num_counts)
% xlabel('time (h)')
% ylabel('total counts')
% title('num count run trend')


%% Bin pulses

data.mcp_tdc.all_ok=data.mcp_tdc.num_ok;
data.mcp_tdc.al_pulses=bin_al_pulses(data.mcp_tdc,anal_opts.atom_laser);

%%
anal_opts.osc_fit.adaptive_freq=true; %estimate the starting trap freq
anal_opts.osc_fit.appr_osc_freq_guess=[52,46.7,40];%[50,12,37];%
anal_opts.osc_fit.plot_fits=true;
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
shaded_ci_lines=0;
color_shaded=[100,100,100]./255;
h=stfig('oscillation plot');
clf
subplot(1,2,1)
box on
hold on
tplotvalues=linspace(0.2,1.8,1e5)';
% tplotvalues=linspace(0.7,1.2,1e5)';
predictor = squeeze(data.mcp_tdc.al_pulses.pos.mean(1,:,:));
% predictor(isnan(predictor))= 0 ;
predictorplot=[tplotvalues,...
    interp1(predictor(:,1),predictor(:,2),tplotvalues),...
    interp1(predictor(:,1),predictor(:,3),tplotvalues)];
[prediction,ci]=predict(data.osc_fit.model{1},predictorplot,'Alpha',1-erf(1/sqrt(2)));
sqrtn=sqrt(col_vec(data.mcp_tdc.al_pulses.num_counts)); %find the statistical uncert in a single shot
xyzerr_tmp=squeeze(data.mcp_tdc.al_pulses.vel_zxy.std(:,:,:))./sqrtn;
xyzerr_tmp(:,sqrtn<2)=nan;

% if shaded_ci_lines
%     patch([predictorplot(:,1)', fliplr(predictorplot(:,1)')], [ci(:,1)', fliplr(ci(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
% else
%     plot(predictorplot(:,1)-data.mcp_tdc.al_pulses.time_cen(1,1),ci(:,1)*1e3,'-','LineWidth',1.5,'Color',color_shaded)
%     plot(predictorplot(:,1)-data.mcp_tdc.al_pulses.time_cen(1,1),ci(:,2)*1e3,'-','LineWidth',1.5,'Color',color_shaded)
% end
time_cen = data.mcp_tdc.al_pulses.time_cen(1,1);
plot(predictorplot(:,1)-data.mcp_tdc.al_pulses.time_cen(1,1),prediction*1e3,'-','LineWidth',1.8,'Color',color_shaded)
errorbar(data.mcp_tdc.al_pulses.time_cen-data.mcp_tdc.al_pulses.time_cen(1,1),1e3.*(data.mcp_tdc.al_pulses.vel_zxy.mean(1,:,anal_opts.osc_fit.dimesion+1)-mean(data.mcp_tdc.al_pulses.vel_zxy.mean(1,:,anal_opts.osc_fit.dimesion+1))),1e3.*xyzerr_tmp(:,anal_opts.osc_fit.dimesion+1),'ko','CapSize',0,'MarkerSize',3,'LineWidth',1.5);
ylabel('$y$ velocity (mm/s)')
% ylabel('$$ velocity (mm/s)')
xlabel('Time (s)')
%title('Velocity Space')
set(gca,'FontSize',17)
xlim([0 1.46])
ylim([-15 20])
% ylim([-45 70])
% xlim([0 0.19])
%%
%stfig('fourier transform oscillations');
stfig('oscillation plot');
tdat = data.mcp_tdc.al_pulses.time_cen-data.mcp_tdc.al_pulses.time_cen(1,1);
tint = linspace(0,1.46911,1e4);
n=5e3;
S=[];
%
for ii = 1:1
S=[S,(data.mcp_tdc.al_pulses.vel_zxy.mean(ii,:,anal_opts.osc_fit.dimesion+1)-mean(data.mcp_tdc.al_pulses.vel_zxy.mean(ii,:,anal_opts.osc_fit.dimesion+1)))];
end
S= exp(data.mcp_tdc.al_pulses.time_cen.'.*0.2964).*S;
% S=prewhiten(S.').';
% n=size(S,2);
% size(S)
% Sinterp=smooth(tdat,S,0.1);%interp1(tdat,S,tint,'spline');
% S = Sinterp;
% S = S(1:101);
Fs = 1./anal_opts.atom_laser.pulsedt;   %1./(tint(2)-tint(1));%         % Sampling frequency
T = 1/Fs;             % Sampling period
L = length(S);             % Length of signal
t = (0:L-1)*T;%(1:L);%        % Time vector
Y = fft(S,n);%fft(S);%
P2 = abs(Y/L);
P1 = P2(1:n/2+1);%P2(1:L/2+1);%
P1(2:end-1) = 2*P1(2:end-1);
f = linspace(0,Fs*0.5,n/2+1);%%Fs*(0:(L/2))/L;%

% Up-Sampling factor:
k = 4;
ias = S.';
[M,N] = size(ias);       
% Zero pad the signal:
ias = [ias ; zeros((k-1)*M,N) ];
M = k*M;
% Time domain:
dt = anal_opts.atom_laser.pulsedt;
t = dt*(0:M-1)';
% Frequency domain:
Fs = 1/dt;
dF = Fs/M;
f2 = (-Fs/2:dF:Fs/2-dF)';
% Discrete Fourier Transform:
Y = fftshift(fft(ias))/M;
% Plot the spectrum:
figure(61);
clf
plot(f2,abs(Y));
hold on
plot(f,P1)

fwhm(f2,abs(Y),0.5)
fwhm(f,P1,0.5)./2.35
fwhm(f,P1,0.667)./2

%%
% subplot(1,panes,pane_counter);
%         pane_counter=pane_counter+1;
% stfig('position hist');
stfig('oscillation plot');
% subplot(1,2,1)
axes('Position',[.25 .73 .2 .16])
box on
% anal_opts.tdc_import.dir = 'Z:\EXPERIMENT-DATA\2020_trap_freq_met\aliasing_trap_freq\aliasing';
shot=72:80;%1:10;
anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
anal_opts.tdc_import.shot_num=shot;
data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);
% clf
bins=502;
spatial_blur=2;
tmin = 0.26;%0.9;%
tmax = 0.56;%1.0;%
ymin = -10e-3;
ymax = 20e-3;
xmin = -20;
xmax =  20;
XEdges=linspace(xmin,xmax,bins);
YEdges=linspace(ymin,ymax,bins);
TEdges=linspace(tmin,tmax,bins);
counts_txy=cell2mat(data.mcp_tdc.counts_txy(:));
bin_area=((ymax-ymin)/bins)*((tmax-tmin)/bins);
[counts,centers]=hist3(counts_txy(:,[3,1]),'edges',{YEdges,TEdges});
counts=counts/bin_area;

if  ~spatial_blur==0
    counts=imgaussfilt(counts,spatial_blur);
end
imagesc(centers{2}-time_cen,centers{1}.*1e3,counts.^1)
% caxis([0 1e7])
% colormap(gca,cmap)
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
title('Position Space Density','FontSize',13)
ylabel('$y$ (mm)','FontSize',13)
xlabel('T (s)','FontSize',13)
%h=colorbar;
%xlabel(h,'Count Density (m^{-1}s^{-1})')
%%
subplot(1,2,2)
plot(f,P1.*1e3,'k-','LineWidth',1.8)
%title('Fourier Space')
xlabel('Frequency (Hz)')
ylabel('Fourier Amplitude (arb.)')
set(gca,'FontSize',17)
set(gca,'XScale','log')
xlim([0 max(f)])
fwhm(f,P1,0.5)./2.35
fwhm(f,P1,0.667)./2
xticks([10 50 100 500 1000])
%%

%%
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'oscillations_unaliased_pdf_2','-dpdf','-r0')

%%

function width = fwhm(x,y,lvl)
% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)
y = y / max(y);
N = length(y);
lev50 = lvl;%0.5;
if y(1) < lev50                  % find index of center (max or min) of pulse
    [garbage,centerindex]=max(y);
    Pol = +1;
    disp('Pulse Polarity = Positive')
else
    [garbage,centerindex]=min(y);
    Pol = -1;
    disp('Pulse Polarity = Negative')
end
i = 2;
while sign(y(i)-lev50) == sign(y(i-1)-lev50)
    i = i+1;
end                                   %first crossing is between v(i-1) & v(i)
interp = (lev50-y(i-1)) / (y(i)-y(i-1));
tlead = x(i-1) + interp*(x(i)-x(i-1));
i = centerindex+1;                    %start search for next crossing at center
while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) & (i <= N-1))
    i = i+1;
end
if i ~= N
    Ptype = 1;
    disp('Pulse is Impulse or Rectangular with 2 edges')
    interp = (lev50-y(i-1)) / (y(i)-y(i-1));
    ttrail = x(i-1) + interp*(x(i)-x(i-1));
    width = ttrail - tlead;
else
    Ptype = 2;
    disp('Step-Like Pulse, no second edge')
    ttrail = NaN;
    width = NaN;
end
end