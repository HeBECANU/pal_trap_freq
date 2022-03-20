clear all
% if contains(getenv('COMPUTERNAME'),'RSPE')
%     data_root = '\\AMPLPC29\He_BEC_Archive\EXPERIMENT-DATA\';
% else
%     data_root = 'E:\Data';
% end
% data_root = '\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
data_root = 'Z:\EXPERIMENT-DATA\2020_trap_freq_met';
% data_dir = fullfile(data_root,'2019_QD\tests\20190917_hf_test_freq_cal_trials_1');
data_folder = 'sample_freq_var\20210611_tune_trap_pal_good_atom_num';
%sample_freq_var\20210611_tune_trap_pal_good_atom_num
%'sample_freq_var\20210525_halo_trap_freq_measure';
%'sample_freq_var\20210615_tune_trap_pal_good_atom_num_2';
%'sample_freq_var\20210611_tune_trap_pal_good_atom_num'
%'sample_freq_var\20210309_tune_trap_freq_met';
%'sample_freq_var\20210317_tighter_trap_scan_improved';%';
%'sample_freq_var\20210610_pal_low_evap';
%'20210609_testing_pal';
%'20210309_tune_trap_freq_met';
%'20210525_halo_trap_freq_measure';%
% data_folder = '';%
data_dir = fullfile(data_root,data_folder);

%%
log_tabl = readtable([data_dir,'\log_freq.txt']);
log_freq=table2array(log_tabl(:,[1,4]));
log_samp=table2array(log_tabl(:,[3]));

unique_freq=unique(log_freq(:,2));
numfiles=length(unique_freq);

tdc_data = cell(numfiles,1);
sample_period = zeros(numfiles,1);

indx_vec = log_freq(:,1);%(1:size(log_freq,1))';%log_freq

for freq_idx = 1:numfiles
    import_opts.dir = data_dir;
    import_opts.force_reimport = true;
    import_opts.shot_num = unique(indx_vec(log_freq(:,2)==unique_freq(freq_idx),1))';
    tdc_data{freq_idx} = import_mcp_tdc_data(import_opts);
    sample_period(freq_idx) = unique(log_samp(log_freq(:,2)==unique_freq(freq_idx),1))+1e-3;%str2double(this_dir(27:end-8))*1e-3;
end
%%
% allfiles = dir(data_dir);
% fnames = {allfiles.name};
% fnames = fnames(~ismember(fnames,{'.','..'}));
% numfiles = numel(fnames);
% tdc_data = cell(numfiles,1);
% sample_period = zeros(numfiles,1);
% for dir_idx = 1:numfiles
%     this_dir = fnames{dir_idx};
%     import_opts.dir = fullfile(data_dir,this_dir);
%     % import_opts.shot_num = 6:10;
%     tdc_data{dir_idx} = import_mcp_tdc_data(import_opts);
%     sample_period(dir_idx) = 1/str2double(this_dir(6:end));%str2double(this_dir(27:end-8))*1e-3;
% end
% 
%%
sample_freq_vec = 1./sample_period;
for ii=1:length(tdc_data)
    num_mask=tdc_data{ii}.num_counts>0;%2e4;
    tdc_data{ii} = struct_mask(tdc_data{ii},num_mask);
end

bin_data = cell(numfiles,1);

for dir_idx = 1:numfiles
    % hist_data = show_txy_raw(tdc_data,'lims',[0.4,.55;-.002,.0045;-.02*1.2,.04*1.2],'num_bins',[300,300,300]);
    % % Stuff for illustration fig
    
    
    % opts - structure with fields
    bin_opts.pulses = 150;%70;%130;%150;%300;%105;%
    bin_opts.pulsedt = sample_period(dir_idx);
    bin_opts.t0 = 0.41778;%1.777;%1.7184;%
    bin_opts.xylim = [-.02,.02;-.025,.025];
    bin_opts.plot = 0;
    bin_opts.start_pulse = 1;
    bin_data{dir_idx} = bin_al_pulses(tdc_data{dir_idx},bin_opts);
end

%%

trap_f_guess = [415,52.8,411.6];%[536,51,534];%[200,50,200];%[407,52,405];%[402,52.2,411.5];%[402,52,411.5]%[75,830,830];
t_max = 3;%1.25;%1;%
num_lim = 30;%30;

% tiledlayout(4,6)

do_ax = [1:3];

sample_freqs=[];
fitmdl = cell(numfiles,3);
fits_done = cell(numfiles,3);

p_stash = zeros(numfiles,3,5);
do_files = 1:numfiles;
ax=1;
for dir_idx = do_files % doing fits
    %
    sample_freq = 1./sample_period(dir_idx);
    nyquist_guess = ceil(trap_f_guess./(.5*sample_freq)); %technically offset by 1
    if abs(sample_freq-52)<0.1
        dum=0;
    end
    for jj = 1:3
        if mod(nyquist_guess(jj),2)==0
            f_apparent_guess(jj) =  (nyquist_guess(jj)).*sample_freq/2-trap_f_guess(jj);
            f_recover_guess(jj) =  (nyquist_guess(jj)).*sample_freq/2-f_apparent_guess(jj);
        else
            f_apparent_guess(jj) = trap_f_guess(jj) -  (nyquist_guess(jj)-1).*sample_freq/2;
            f_recover_guess(jj) = f_apparent_guess(jj) +  (nyquist_guess(jj)-1).*sample_freq/2;
        end
    end
    if isempty(bin_data{dir_idx}.num_counts)
       fitmdl{dir_idx,ax}{ii} = nan;%fit_obj;
       fits_done{dir_idx,ax}{1} = 0;
       continue
    end
    bin_data_curr = bin_data{dir_idx}.vel_zxy.mean-nanmean(bin_data{dir_idx}.vel_zxy.mean,[1,2]);
    bin_data_std = bin_data{dir_idx}.vel_zxy.std;
    
    num_curr = mean(bin_data{dir_idx}.num_counts);
    
    %     com_disp = 1e3*squeeze(nanmean(bin_data{dir_idx}.vel_zxy.mean-nanmean(bin_data{dir_idx}.vel_zxy.mean,[1,2]),1));
    for ii = 1:size(bin_data{dir_idx}.vel_zxy.mean,1)
        com_disp = 1e3*squeeze(bin_data_curr(ii,:,:));
        com_std = 1e3*squeeze(bin_data_std(ii,:,:));
        
        t_vals = bin_data{dir_idx}.time_cen-min(bin_data{dir_idx}.time_cen);
        t_fine = linspace(0,t_max,300)';
        data_mask = t_vals<t_max & num_curr'>num_lim;
        
        t_vals = t_vals(data_mask);
        com_disp = (com_disp(data_mask,:));
        com_std = (com_std(data_mask,:));
        
        weights = 1./com_std(:,3).^2;
        w_mask = 1.5.*mean(weights)>weights;
        
        A_guess = max(com_disp,[],1);
        time_constant = [.2,.2,.2];
        if isempty(t_vals) || size(t_vals,1)<8
            for ax = do_ax
            fitmdl{dir_idx,ax}{ii} = fit_obj;%nan;%
            fits_done{dir_idx,ax}{ii} = 0;
            peakfourier{dir_idx,ax}{ii} = nan;
            end
        else
            
        mdl_plot = A_guess.*exp(-t_fine./time_constant).*sin(2*pi*f_apparent_guess.*t_fine);
        
        % amp tau freq phase trend
        
        for ax = do_ax
            %                p_guess = [A_guess(ax), time_constant(ax),100,0,0];
            p_guess = [A_guess(ax), time_constant(ax),f_apparent_guess(ax),0,0];
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0,0,0,-pi,-inf],...
                'Upper',[30,15,sample_freq/2,pi,inf],...
                'StartPoint',p_guess,...
                'MaxFunEvals',3e3);
            p_stash(dir_idx,ax,:) = p_guess;
            ft = fittype('a.*exp(-x./b).*sin(2*pi*c*x+d)+e','options',fo);
            com_disp(isnan(com_disp(:,ax)),ax)=0;
            
            sin_decay = @(b,x) b(1).*exp(-x./b(2)).*sin(2*pi*b(3)*x+b(4))+b(5)+b(6).*x;
            
            %         fit_obj_nlm = fitnlm(t_vals,com_disp(:,ax),sin_decay,p_guess,'Weights',weights);
            [fit_obj,res,output] = fit(t_vals(2:end),com_disp(2:end,ax),ft,fo);
            
%             p_guess = [A_guess(ax), time_constant(ax),f_apparent_guess(ax),0,0,0];
%             fo = fitoptions('Method','NonlinearLeastSquares',...
%                 'Lower',[0,0,0,-pi,-inf,-inf],...
%                 'Upper',[30,15,sample_freq/2,pi,inf,inf],...
%                 'StartPoint',p_guess,...
%                 'MaxFunEvals',3e3);
%             ft = fittype('a.*exp(-x./b).*sin(2*pi*c*x+d)+e+f*x','options',fo);
%             [fit_obj,res,output] = fit(t_vals(2:end),com_disp(2:end,ax),ft,fo);
%             
            
            n=5e4;%50e3;%15e3;%5e3;
    S=com_disp(2:end,ax);
%     if ax == 2
%     S=com_disp(2:end,ax)+6.67.*t_vals(2:end);
%     end
Fs = sample_freq;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(S);             % Length of signal
t = (0:L-1)*T;        % Time vector
Y = fft(S,n);
P2 = abs(Y/L);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = linspace(0,Fs*0.5,n/2+1);%Fs*(0:(L/2))/L;
[val,indx]=max(P1);
peakfourier{dir_idx,ax}{ii} = f(indx);

P{dir_idx,ax}{ii} = P1;
freq_vec{dir_idx,ax} = f;
            
%             size(com_disp(2:end,ax))
            
            ft2 = fittype('a.*sin(2*pi*c*x+d)+b','options',fo);
            p_guess = [A_guess(ax), 0,f_apparent_guess(ax),0];
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0,-15,0,-pi],...
                'Upper',[30,15,sample_freq/2,pi],...
                'StartPoint',p_guess,...
                'MaxFunEvals',3e3);
            [fit_obj_2,res_2] = fit(t_vals(:),com_disp(:,ax),ft2,fo);
            
            ci1 = confint(fit_obj,0.95);
            ci2 =confint(fit_obj_2,0.95);
            
            if res_2.rsquare>res.rsquare || (range(ci1(:,3))>range(ci2(:,3)) && res_2.rsquare>0.9)
                res=res_2;
                fit_obj = fit_obj_2;
            end
            
            if res.rsquare<0.9
                p_guess = [A_guess(ax), time_constant(ax),f_apparent_guess(ax)./2,0,0];
                fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',[0,0,0,-pi],...
                    'Upper',[30,15,sample_freq/2,pi],...
                    'StartPoint',p_guess,...
                    'MaxFunEvals',3e3);
                [fit_obj_half,res_half] = fit(t_vals,com_disp(:,ax),ft,fo);
                if res_half.rsquare>res.rsquare
                    fit_obj = fit_obj_half;
                    res=res_half;
                else
                    p_guess = [A_guess(ax), time_constant(ax),f_apparent_guess(ax)-10,0,0];
                    fo = fitoptions('Method','NonlinearLeastSquares',...
                        'Lower',[0,0,0,-pi,-inf],...
                        'Upper',[30,15,sample_freq/2,pi,inf],...
                        'StartPoint',p_guess,...
                        'MaxFunEvals',3e3);
                    [fit_obj_half,res_half] = fit(t_vals,com_disp(:,ax),ft,fo);
                    if res_half.rsquare>res.rsquare
                        fit_obj = fit_obj_half;
                        res=res_half;
                    end
                end
            end
%             if sample_freq>130 && ax==2 && res.rsquare>0.8
%                 dum=0;
%             end
            fitmdl{dir_idx,ax}{ii} = fit_obj;
            fits_done{dir_idx,ax}{ii} = 1;
        end
        end
        sample_freqs = [sample_freqs,sample_freq_vec(dir_idx)];
    end
end
cli_header('Done.');
%%
C=fitmdl(:,do_ax);
D=peakfourier(:,do_ax);
clear AllC freq_fits fourier_val
for ii = 1:length(do_ax)
    AllC(ii) = {cat(2, C{:,ii})};
    AllD(ii,:) = cat(2, D{:,ii});
    freq_fits(ii,:) = cellfun(@(x) x.c, AllC{ii});
    fourier_val(ii,:) = cell2mat(AllD(ii,:));
    
    for jj = 1:size(P,1)
        freq_spec(ii,jj,:) = mean(cell2mat(P{jj,ii}),2);
        ft = freq_vec{jj,ii};
        if ii == 2 && (1./sample_period(jj)<104) && (1./sample_period(jj)>51) && ~(abs(1./sample_period(jj)-100)<0.5)
            [valt,indxt] = max(squeeze(freq_spec(ii,jj,:))./(ft+10).');%max(squeeze(freq_spec(ii,jj,:)));%
        elseif ii==2 && (abs(1./sample_period(jj)-100)<0.5)
            [valt,indxt] = max(squeeze(freq_spec(ii,jj,ft<48)));%max(squeeze(freq_spec(ii,jj,:)));%
        else
            [valt,indxt] = max(squeeze(freq_spec(ii,jj,:)));%max(squeeze(freq_spec(ii,jj,:)));%
        end
        fourier_avg_val(ii,jj) = ft(indxt);
        fourier_avg_std(ii,jj) = 9;
    end
    
end

%%
confidence_intervals = [];
freq_fits_avg = [];
is_shot_good = ones(size(freq_fits,2),length(do_ax));%
ci_tol = [10e-2,5e-2,10e-2];%[50e-6*402,50e-6*250,50e-6*402].*10000;%[50e-6*402,50e-6*402,50e-6*402].*50;%50.0;
ax_indx = 1;
for ax = do_ax
    shot_indx = 1;
    for ii = 1:numfiles
        for jj = 1:length(fitmdl{ii,ax})
            ci = confint(fitmdl{ii,ax}{jj},0.682689492);%0.95
            %         range(ci(1:2,3))./401*1e6/1.96
            %         range(ci(1:2,3))./401*1e6
            confidence_intervals(shot_indx,:,ax) = ci(1:2,3)';
            if range(ci(1:2,3))>ci_tol(ax) || ~fits_done{ii,ax}{jj} || isnan(range(ci(1:2,3)))
                is_shot_good(shot_indx,ax_indx) = 0;
                confidence_intervals(shot_indx,:,ax) = nan.*ci(1:2,3)';
            end
            
            
            shot_indx = shot_indx+1;
        end
        curr_indx = (shot_indx-length(fitmdl{ii,ax})):(shot_indx-1);
        curr_freq=fourier_val(ax_indx,curr_indx);%freq_fits(ax_indx,curr_indx);
        curr_freq=curr_freq(logical(is_shot_good(curr_indx,ax_indx)));
        freq_fits_avg(ax,ii) = nanmean(curr_freq);
        freq_fits_std(ax,ii) = nanmean(range(confidence_intervals(curr_indx,:,ax)));%nanstd(curr_freq);
    end
    ax_indx = ax_indx+1;
end
%%
% shunt_freq_fits = cellfun(@(x) x.c, shunt_fitmdl);
sawfit = cell(size(do_ax,2),1);
f_s = linspace(50,max(sample_freqs).*1.1,15e4);
axcol = {'r','b','k'};
axstyl = {'-','--','-.'};
axshp = {'o','s','^'};
plabel= [];
stfig('Sawtooth plot');
clf
ax_done = do_ax;
freq_vals = fourier_val;%freq_fits;%
is_shot_good = logical(is_shot_good);
fit_ax=2;
weak_ax_fit = fitnlm(sample_freq_vec,fourier_avg_val(fit_ax,:),@nyquist_sawtooth,trap_f_guess(fit_ax));
% trap_f_guess = [420,50,420];
% old_mask = (sample_freqs==[60,70,80,90]);
for plt = 1:length(do_ax)%2%
    ax_c = ax_done(plt);
    data_mask = is_shot_good(:,plt);% & (sample_freqs>144.1).';% &  sample_freqs>80;
    sawfit{plt} = fitnlm(sample_freqs(data_mask),freq_vals(plt,data_mask),@nyquist_sawtooth,trap_f_guess(ax_c));
    %     subplot(2,1,plt)
    hold on
    plabel(plt) = plot(sample_freqs(data_mask),freq_vals(plt,data_mask),'.','Color',axcol{ax_done(plt)},...
        'MarkerSize',15);
    plot(f_s,nyquist_sawtooth(sawfit{plt}.Coefficients.Estimate(1),f_s),...
        'Color',axcol{ax_done(plt)})%401
    errorbar(sample_freqs(data_mask),freq_vals(plt,data_mask),freq_vals(plt,data_mask)-confidence_intervals(data_mask,1,ax_c)',confidence_intervals(data_mask,2,ax_c)'-freq_fits(plt,data_mask),[axcol{ax_done(plt)},'.'])
end

% trap_f_guess = [830,75,830];
% for plt = 2
%     x_sawfit{plt} = fitnlm(x_sample_freq,shunt_freq_fits(:,plt),@nyquist_sawtooth,trap_f_guess(plt));
% %     subplot(2,1,plt)
%     hold on
%     plabel(plt) = plot(x_sample_freq,shunt_freq_fits(:,plt),'.','Color',axcol{plt},...
%         'MarkerSize',15);
%     plot(f_s,nyquist_sawtooth(x_sawfit{plt}.Coefficients.Estimate(1),f_s),...
%         'Color',axcol{do_ax(plt)})
% end

labels = {'Z','X','Y'};
legend(plabel,labels(do_ax),'Location','NorthWest')
% legend(plabel([1,3]),{'Z','Y'},'Location','NorthWest')
% legend(plabel([3]),{'Y'},'Location','NorthWest')
xlabel('Sampling frequency, $f_s$ (Hz)')
ylabel('Apparent frequency, $f_a$ (Hz)')
set(gca,'FontSize',16)
ylim([0 100])
%%
h=stfig('Sawtooth');
clf
subplot(1,2,1)
title('Multiple Zones')
f_s = linspace(45,max(sample_freqs).*1.1,15e4);
avg_plot_vals = fourier_avg_val;%freq_fits_avg
ax_done = [2,3,1];
for plt = 1:length(do_ax)%[3,1,2]%
    ax_c = ax_done(plt);
    data_mask = is_shot_good(:,plt);% & sample_freqs<160 &  sample_freqs>80;
    %     subplot(2,1,plt)
    hold on
    plabel(plt) = plot([1e5],[1],axshp{ax_done(plt)},'Color',axcol{ax_done(plt)},...
        'MarkerSize',8,'MarkerFaceColor',axcol{ax_done(plt)},'LineWidth',1.3);

    plot(sample_freq_vec,avg_plot_vals(ax_c,:),axshp{ax_done(plt)},'Color',axcol{ax_done(plt)},...
        'MarkerSize',4,'MarkerFaceColor',axcol{ax_done(plt)},'LineWidth',1.3);
    if ax_c == 2
        plot(f_s,nyquist_sawtooth(weak_ax_fit.Coefficients.Estimate(1),f_s),axstyl{ax_done(plt)},...
        'Color',axcol{ax_done(plt)},'LineWidth',1.3)%401
    else
    plot(f_s,nyquist_sawtooth(sawfit{ax_c}.Coefficients.Estimate(1),f_s),axstyl{ax_done(plt)},...
        'Color',axcol{ax_c},'LineWidth',1.3)%401
    end
    errorbar(sample_freq_vec,avg_plot_vals(ax_c,:),freq_fits_std(ax_c,:),[axcol{ax_done(plt)},'o'],...
        'CapSize',0,'MarkerSize',1,'LineWidth',1.5)
end

% trap_f_guess = [830,75,830];
% for plt = 2
%     x_sawfit{plt} = fitnlm(x_sample_freq,shunt_freq_fits(:,plt),@nyquist_sawtooth,trap_f_guess(plt));
% %     subplot(2,1,plt)
%     hold on
%     plabel(plt) = plot(x_sample_freq,shunt_freq_fits(:,plt),'.','Color',axcol{plt},...
%         'MarkerSize',15);
%     plot(f_s,nyquist_sawtooth(x_sawfit{plt}.Coefficients.Estimate(1),f_s),...
%         'Color',axcol{do_ax(plt)})
% end
box on
labels = {'X','Y','Z'};%{'Z','X','Y'};
legend(plabel(:),labels(do_ax),'Location','Best')%[2,3,1]
% legend(plabel([1,3]),{'Z','Y'},'Location','NorthWest')
% legend(plabel([3]),{'Y'},'Location','NorthWest')
xlabel('Sampling frequency, $f_s$ (Hz)')
ylabel('Apparent frequency, $f_a$ (Hz)')
set(gca,'FontSize',20)
ylim([0 80])
xlim([48 145])
% %

%
%sawfit = cell(size(do_ax,2),1);
f_s = linspace(50,max(sample_freqs).*1.1,15e4);
axcol = {'r','b','k'};
plabel= [];
%stfig('Sawtooth plot narrow');
subplot(1,2,2)
% clf
ax_done = do_ax;
freq_vals = fourier_val;%freq_fits;%
is_shot_good = logical(is_shot_good);
fit_ax=1;
fitnlm(sample_freq_vec,fourier_avg_val(fit_ax,:),@nyquist_sawtooth,trap_f_guess(fit_ax));
frange = [120 137];
title('Single Zone')
% trap_f_guess = [420,50,420];
% old_mask = (sample_freqs==[60,70,80,90]);
for plt = 1:length(do_ax)%2%
    ax_c = ax_done(plt);
    data_mask = is_shot_good(:,plt) & (sample_freqs>frange(1)).' &  (sample_freqs<frange(2)).';
    sawfit_narrow{plt} = fitnlm(sample_freq_vec,avg_plot_vals(ax_c,:),@nyquist_sawtooth,trap_f_guess(ax_c));
    %     subplot(2,1,plt)
    hold on
    plabel(plt) = plot(sample_freq_vec,avg_plot_vals(ax_c,:),axshp{ax_done(plt)},'Color',axcol{ax_done(plt)},...
        'MarkerSize',6,'MarkerFaceColor',axcol{ax_done(plt)});
    if plt == 2
        plot(f_s,nyquist_sawtooth(weak_ax_fit.Coefficients.Estimate(1),f_s),axstyl{ax_done(plt)},...
        'Color',axcol{ax_done(plt)},'LineWidth',1.3)%401
    else
    plot(f_s,nyquist_sawtooth(sawfit{plt}.Coefficients.Estimate(1),f_s),axstyl{ax_done(plt)},...
        'Color',axcol{ax_done(plt)},'LineWidth',1.3)%401
    end
    errorbar(sample_freq_vec,avg_plot_vals(ax_c,:),freq_fits_std(ax_c,:),[axcol{ax_done(plt)},'o'],...
        'CapSize',0,'MarkerSize',3,'LineWidth',1.5)
end
box on
xlim([frange(1)-1,frange(2)])
xlabel('Sampling frequency, $f_s$ (Hz)')
ylabel('Apparent frequency, $f_a$ (Hz)')
set(gca,'FontSize',20)

%%
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'sawtooth_pdf_v2','-dpdf','-r0')

%%

function Y = damped_sine(p,t)
amp = p(1);
decay = exp(-t/p(2));
freq = p(3);
phase = p(4);
trend = p(5);
Y = amp*decay.*sin(2*pi*freq*t+phase) + trend*t;

end

function f_app = nyquist_sawtooth(f_act,f_samp)
nyquist_zone = ceil(f_act./(.5*f_samp)); %technically offset by 1
zone_even = mod(nyquist_zone,2) == 0;
f_app = 0*f_samp;
f_app(zone_even) =  (nyquist_zone(zone_even)).*f_samp(zone_even)/2-f_act;
f_app(~zone_even) = f_act -  (nyquist_zone(~zone_even)-1).*f_samp(~zone_even)/2;
end


% %
% fsize = 16;
% stfig('PAL fits');
% clf
% h = 4;
% w = 6;
% tiledlayout(4,6)
% ctr = 0;
%%
% for dir_idx = do_files % plot each fit
%     T = linspace(0,.35,100);
%
%     com_disp = 1e3*squeeze(nanmean(bin_data{dir_idx}.vel_zxy.mean-nanmean(bin_data{dir_idx}.vel_zxy.mean,[1,2]),1));
%
%     t_vals = bin_data{dir_idx}.time_cen-min(bin_data{dir_idx}.time_cen);
%     t_fine = linspace(0,t_max,300)';
%     t_vals = t_vals(t_vals<t_max);
%     com_disp = (com_disp(t_vals<t_max,:));
%     ctr = ctr + 1;
%
%     nexttile
%     hold on
%     for ax = 3
%         if fits_done(dir_idx,ax)
%             plot(t_vals,damped_sine(squeeze(p_stash(dir_idx,ax,:))',t_vals),':','Color',axcol{ax})
%             plot(t_vals,com_disp(:,ax),'.','Color',axcol{ax});
%             p=plot(fitmdl{dir_idx,ax});
%             set(p,'Color',axcol{ax})
%             r=plot(fitmdl{dir_idx,ax},t_vals,com_disp(:,ax),'x','residuals');
%             set(r,'Color',axcol{ax})
%             legend('off')
%         else
%             plot(t_fine,mdl_plot(:,ax),':','Color',axcol{ax})
%         end
%     end
%
%     xlim([0,t_max])
%     xlabel('Time (s)')
%     ylabel('Mean velocity (mm/s)')
% end
% cli_header('Done.');

% %

%%

%% Now need to look at the X axis too...
% data_x = '2019_QD\tests\20191107_hf_trap_freq_cal_redo\trap_shunt_method';
% x_import_opts = [];
% x_import_opts.dir = fullfile(data_root,data_x);
% num_x_files = 6;
% x_tdc_data = cell(num_x_files,1);
% for x_samp = 1:num_x_files
%     x_import_opts.shot_num = 5*(x_samp-1)+[1:5];
%     x_tdc_data{x_samp} = import_mcp_tdc_data(x_import_opts);
% end
% % %
%
% x_sample_period = [0.011,0.0112,0.0114,0.01145,0.0115,0.01155];
% x_sample_freq = 1./x_sample_period;
%
% %%
% x_bin_data = cell(num_x_files,1);
% for dir_idx = 1:num_x_files
%     % hist_data = show_txy_raw(tdc_data,'lims',[0.4,.55;-.002,.0045;-.02*1.2,.04*1.2],'num_bins',[300,300,300]);
%     % % Stuff for illustration fig
%     % opts - structure with fields
%     bin_opts.pulses = 150;
%     bin_opts.pulsedt = x_sample_period(dir_idx);
%     bin_opts.t0 = 0.4178;
%     bin_opts.xylim = [-.01,.01;-.02,.02];
%     bin_opts.plot = 1;
%     bin_opts.start_pulse = 1;
%     x_bin_data{dir_idx} = bin_al_pulses(x_tdc_data{dir_idx},bin_opts);
% end
%
% %%
% t_max = 0.25;
% f_x_stash = [20,20,20;
%                 20,20,20;
%                 20,20,20;
%                 35,20,20;
%                 35,20,20;
%                 40,20,20];
%
% axcol = {'r','b','k'};
% fsize = 16;
% stfig('PAL fits, x');
% clf
% h = 4;
% w = 6;
% tiledlayout(2,3)
% ctr = 0;
% do_ax = 1:3;
% shunt_fitmdl = cell(num_x_files,length(do_ax));
% for dir_idx = 1:num_x_files % plot each fit
%     T = linspace(0,.35,100);
%
%     com_disp = 1e3*squeeze(nanmean(x_bin_data{dir_idx}.vel_zxy.mean-nanmean(x_bin_data{dir_idx}.vel_zxy.mean,[1,2]),1));
%
%     t_vals = x_bin_data{dir_idx}.time_cen-min(x_bin_data{dir_idx}.time_cen);
%     t_fine = linspace(0,t_max,300)';
%     t_vals = t_vals(t_vals<t_max);
%     com_disp = (com_disp(t_vals<t_max,:));
%     ctr = ctr + 1;
%
%     nexttile
%     hold on
%
%     for ax = do_ax
%         plot(t_vals,com_disp(:,ax),'.','Color',axcol{ax});
% %         for ax = do_ax
%            p_guess = [max(com_disp(:,ax)),0.2,f_x_stash(dir_idx,ax),0,0];
%            fo = fitoptions('Method','NonlinearLeastSquares',...
%                    'Lower',[0,0,0,-pi,-inf],...
%                    'Upper',[20,2,x_sample_freq(dir_idx)/2,pi,inf],...
%                    'StartPoint',p_guess,...
%                    'MaxFunEvals',1e3);
% %             p_stash(dir_idx,ax,:) = p_guess;
%             ft = fittype('a.*exp(-x./b).*sin(2*pi*c*x+d)+e*x','options',fo);
%             [fit_obj,~] = fit(t_vals,com_disp(:,ax),ft,fo);
%             shunt_fitmdl{dir_idx,ax} = fit_obj;
%
%            fits_done(dir_idx,ax) = 1;
%
% %             f = fitnlm(t_vals,com_disp(:,ax),@damped_sine,p_guess);
% %             [f_pred,pred_ci]= predict(f, t_fine);
% %             plot(t_fine,f_pred)
%             p=plot(shunt_fitmdl{dir_idx,ax});
%             set(p,'Color',axcol{ax})
%             r=plot(shunt_fitmdl{dir_idx,ax},t_vals,com_disp(:,ax),'x','residuals');
%             set(r,'Color',axcol{ax})
%             legend('off')
% %         else
% %             plot(t_fine,mdl_plot(:,ax),':','Color',axcol{ax})
%             plot(t_fine,damped_sine(p_guess,t_fine),':','Color',axcol{ax})
% %         end
%     end
%
%     xlim([0,t_max])
%     xlabel('Time (s)')
%     ylabel('Mean velocity (mm/s)')
% end
% cli_header('Done.');
% function Y = damped_sine_cfit(p,t,)
%     ft = fittype('a*x+b*x.^2','indep','x');me
%     mdl = fit(x,y,ft,'lower',[2 -inf],'upper',[4,inf])
% end