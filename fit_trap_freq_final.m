clear all
% if contains(getenv('COMPUTERNAME'),'RSPE')
%     data_root = '\\AMPLPC29\He_BEC_Archive\EXPERIMENT-DATA\';
% else
%     data_root = 'E:\Data';
% end
data_root = '\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
% data_dir = fullfile(data_root,'2019_QD\tests\20190917_hf_test_freq_cal_trials_1');
data_dir = fullfile(data_root,'20210309_tune_trap_freq_met');
allfiles = dir(data_dir);
fnames = {allfiles.name};
fnames = fnames(~ismember(fnames,{'.','..'}));
numfiles = numel(fnames);

tdc_data = cell(numfiles,1);
sample_period = zeros(numfiles,1);


for dir_idx = 1:numfiles
    this_dir = fnames{dir_idx};
    import_opts.dir = fullfile(data_dir,this_dir);
    % import_opts.shot_num = 6:10;
    tdc_data{dir_idx} = import_mcp_tdc_data(import_opts);
    sample_period(dir_idx) = 1/str2double(this_dir(6:end));%str2double(this_dir(27:end-8))*1e-3;
end
sample_freqs = 1./sample_period;
%%


bin_data = cell(numfiles,1);

for dir_idx = 1:numfiles
    % hist_data = show_txy_raw(tdc_data,'lims',[0.4,.55;-.002,.0045;-.02*1.2,.04*1.2],'num_bins',[300,300,300]);
    % % Stuff for illustration fig
    
    
    % opts - structure with fields
    bin_opts.pulses = 50;
    bin_opts.pulsedt = sample_period(dir_idx);
    bin_opts.t0 = 0.41778;
    bin_opts.xylim = [-.01,.01;-.02,.02];
    bin_opts.plot = 1;
    bin_opts.start_pulse = 1;
    bin_data{dir_idx} = bin_al_pulses(tdc_data{dir_idx},bin_opts);
end

%%

trap_f_guess = [424,40,424];%[75,830,830];
t_max = .3;

% tiledlayout(4,6)

do_ax = [1,2,3];
fitmdl = cell(numfiles,2);
fits_done = zeros(numfiles,3);
f_stash = [   24.0385         0   20.9752
    38.2038         0   42.3726
    16.3271         0   20.2364
    23.5216         0   37.5315
    24.8199         0   31.3806
    10.6151         0   38.3601
    
    40.0000         0   17.7216
    14.4578         0   24.3460
    15.6994         0   30.5773
    23.0387         0   17.4456
    29.9020         0   33.8337
    34.4123         0   33.0226
    
    26.3803         0   23.2230
    20.3224         0   30.1391
    20.9501         0   17.8599
    22.8809         0   24.6410
    18.6935         0   14.9597
    20.7246         0    0.9065
    
    15.6430         0    7.4652
    2.8386         0   12.6144
    2.9807         0    4.1847
    15.7081         0   11.6420
    35.0000         0    9.9874
    11.2864         0   14.8465];
p_stash = zeros(numfiles,3,5);
do_files = 1:numfiles;
for dir_idx = do_files % doing fits
    T = linspace(0,.35,100);
    %
    sample_freq = 1./sample_period(dir_idx);
    if sample_freq == 169
        dum=0;
    end
    nyquist_guess = ceil(trap_f_guess./(.5*sample_freq)); %technically offset by 1
    for jj = 1:3
        if mod(nyquist_guess(jj),2)==0
            f_apparent_guess(jj) =  (nyquist_guess(jj)).*sample_freq/2-trap_f_guess(jj);
            f_recover_guess(jj) =  (nyquist_guess(jj)).*sample_freq/2-f_apparent_guess(jj);
        else
            f_apparent_guess(jj) = trap_f_guess(jj) -  (nyquist_guess(jj)-1).*sample_freq/2;
            f_recover_guess(jj) = f_apparent_guess(jj) +  (nyquist_guess(jj)-1).*sample_freq/2;
        end
    end
    bin_data_curr = bin_data{dir_idx}.vel_zxy.mean-nanmean(bin_data{dir_idx}.vel_zxy.mean,[1,2]);
    
%     com_disp = 1e3*squeeze(nanmean(bin_data{dir_idx}.vel_zxy.mean-nanmean(bin_data{dir_idx}.vel_zxy.mean,[1,2]),1));
    com_disp = 1e3*squeeze(bin_data_curr(end,:,:));
    
    t_vals = bin_data{dir_idx}.time_cen-min(bin_data{dir_idx}.time_cen);
    t_fine = linspace(0,t_max,300)';
    t_vals = t_vals(t_vals<t_max);
    com_disp = (com_disp(t_vals<t_max,:));
    
    A_guess = max(com_disp,[],1);
    time_constant = [.2,.2,.2];
    mdl_plot = A_guess.*exp(-t_fine./time_constant).*sin(2*pi*f_apparent_guess.*t_fine);
    
    % amp tau freq phase trend
    
    for ax = do_ax
        %        p_guess = [A_guess(ax), time_constant(ax),f_stash(dir_idx,ax),0,0];
        p_guess = [A_guess(ax), time_constant(ax),f_apparent_guess(ax),0,0];
        fo = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0,0,-pi,-inf],...
            'Upper',[30,15,sample_freq/2,pi,inf],...
            'StartPoint',p_guess,...
            'MaxFunEvals',3e3);
        p_stash(dir_idx,ax,:) = p_guess;
        ft = fittype('a.*exp(-x./b).*sin(2*pi*c*x+d)+e*x','options',fo);
        com_disp(isnan(com_disp(:,ax)),ax)=0;
        [fit_obj,res] = fit(t_vals,com_disp(:,ax),ft,fo);
        
        ft2 = fittype('a.*sin(2*pi*c*x+d)+b','options',fo);
        p_guess = [A_guess(ax), 0,f_apparent_guess(ax),0];
        fo = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,-15,0,-pi],...
            'Upper',[30,15,sample_freq/2,pi],...
            'StartPoint',p_guess,...
            'MaxFunEvals',3e3);
        [fit_obj_2,res_2] = fit(t_vals,com_disp(:,ax),ft2,fo);
        
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
                end
            end
        end
        fitmdl{dir_idx,ax} = fit_obj;
        
        fits_done(dir_idx,ax) = 1;
    end
    
end
cli_header('Done.');

% %
axcol = {'r','b','k'};
fsize = 16;
stfig('PAL fits');
clf
h = 4;
w = 6;
tiledlayout(4,6)
ctr = 0;
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
%%
confidence_intervals = [];
is_shot_good = ones(size(fitmdl,1),1);
ci_tol = 2.0;
for ax = do_ax
    for ii = 1:numfiles
        ci = confint(fitmdl{ii,ax},0.95);
        if range(ci(1:2,3))>ci_tol
            is_shot_good(ii) = 0;
        end
        confidence_intervals(ii,:,ax) = ci(1:2,3)';
    end
end
%%
freq_fits = cellfun(@(x) x.c, fitmdl(:,[1,2,3]));
% shunt_freq_fits = cellfun(@(x) x.c, shunt_fitmdl);
sawfit = cell(size(do_ax,2),1);
f_s = linspace(60,200,15e4);
plabel= [];
stfig('Sawtooth plot');
clf
ax_done = do_ax;
is_shot_good = logical(is_shot_good);
% old_mask = (sample_freqs==[60,70,80,90]);
data_mask = is_shot_good & sample_freqs<146;% &  sample_freqs>100;
for plt = [1:3]
    ax_c = ax_done(plt);
    sawfit{plt} = fitnlm(sample_freqs(data_mask),freq_fits(data_mask,plt),@nyquist_sawtooth,trap_f_guess(ax_c));
    %     subplot(2,1,plt)
    hold on
    plabel(ax_done(plt)) = plot(sample_freqs(data_mask),freq_fits(data_mask,plt),'.','Color',axcol{ax_done(plt)},...
        'MarkerSize',15);
    plot(f_s,nyquist_sawtooth(sawfit{plt}.Coefficients.Estimate(1),f_s),...
        'Color',axcol{ax_done(plt)})
    errorbar(sample_freqs(data_mask),freq_fits(data_mask,plt),freq_fits(data_mask,plt)-confidence_intervals(data_mask,1,ax_c),confidence_intervals(data_mask,2,ax_c)-freq_fits(data_mask,plt),[axcol{ax_done(plt)},'.'])
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


% legend(plabel,{'Z','X','Y'},'Location','NorthWest')
legend(plabel([1,3]),{'Z','Y'},'Location','NorthWest')
% legend(plabel([3]),{'Y'},'Location','NorthWest')
xlabel('Sampling frequency (Hz)')
ylabel('Apparent frequency (Hz)')
set(gca,'FontSize',16)
ylim([0 100])
%%

% %

%%
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
% function Y = damped_sine_cfit(p,t,)
%     ft = fittype('a*x+b*x.^2','indep','x');me
%     mdl = fit(x,y,ft,'lower',[2 -inf],'upper',[4,inf])
% end