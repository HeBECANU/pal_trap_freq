import_opts.dir = 'Z:\EXPERIMENT-DATA\2020_trap_freq_met\20200911_kick_drop_standard_trap';
data = import_mcp_tdc_data(import_opts);

%% Get log
log_dir = import_opts.dir;
log_file = 'log_LabviewMatlab.txt';
tf = fopen(fullfile(log_dir,log_file));
logdata = textscan(tf,'%s%s%u%s%f%s%f%s%f','Delimiter',',');
wait_times = logdata{9};

%%
num_shots = length(data.shot_num);
bec_cens = zeros(num_shots,3);
for shot_idx = 1:num_shots
    lims = [1.26,1.3;
            -.03,.03;
            -.03,.03];
    this_txy = masktxy_square(data.counts_txy{shot_idx},lims);
    bec_cens(shot_idx,:) = mean(this_txy);
end
bec_cens = bec_cens - mean(bec_cens);
%%
stfig('Plots');
clf
axlabels = ['T','X','Y'];
for axis = 1:3
    subplot(4,1,axis)
    plot(wait_times,bec_cens(:,axis),'*')
    xlabel('Wait time (ms)')
    ylabel(sprintf('%s axis displacement',axlabels(axis)))
    set(gca,'FontSize',14)
    box off
end
% %
subplot(4,1,4)
plot(logdata{9},vecnorm(bec_cens,2,2),'*')
xlabel('Wait time (ms)')
ylabel('Vecnorm displacement')
set(gca,'FontSize',14)
box off
suptitle('Kick-Drop oscillations (standard trap)')
