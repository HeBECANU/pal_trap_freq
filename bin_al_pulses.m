function al_pulses = bin_al_pulses(data, opts)
% adapted from https://github.com/brycehenson/Tune_out_v2/blob/master/lib/bin_al_pulses.m

% Inputs
% data - structure with fields
    % counts_txy
% opts - structure with fields
%     pulses - number of pulses to use
%     pulsedt - gap between pulses
%     t0 - centre of first pulse
%     xylim - rectangular detector region limits
%     plot - visual output enabll
%     c - constants, including gravity g0 and fall_time
%     start_pulse - first pulse index??


iimax = numel(data.counts_txy);
al_pulses = [];
al_pulses.pulsedt = opts.pulsedt;
al_pulses.window = nan(opts.pulses, 3, 2); %initalize
al_pulses.num_counts = nan(iimax, opts.pulses);
al_pulses.pos.mean = nan(iimax, opts.pulses, 3);
al_pulses.pos.std = nan(iimax, opts.pulses, 3);
al_pulses.vel.mean = nan(iimax, opts.pulses, 3);
al_pulses.vel.std = nan(iimax, opts.pulses, 3);

hebec_constants
opts.c = const;
opts.c.tof = sqrt(2*opts.c.fall_distance/opts.c.g0);

fprintf('binning pulses in files %04u:%04u', iimax, 0)
first_good_shot = true;
for shot = 1:iimax
    

    for pulse = 1:opts.pulses
        %set up time window centered arround t0
        t_pulse_cen = opts.t0 + opts.pulsedt ...
                * (opts.start_pulse + pulse - 2);
            trange = t_pulse_cen + opts.pulsedt * [-0.5, 0.5];
        pulse_win_txy = [trange; opts.xylim];
        counts_pulse = masktxy_square(data.counts_txy{shot}, pulse_win_txy);
        if opts.plot>1
            stfig;
            set(gcf, 'Color', [1, 1, 1]);
            subplot(3, 1, 1)
            histogram(counts_pulse(:, 1), 100)
            xlabel('t')
            title('full')
            subplot(3, 1, 2)
            histogram(counts_pulse(:, 2), 100)
            xlabel('x')
            title('full')
            subplot(3, 1, 3)
            histogram(counts_pulse(:, 3), 100)
            xlabel('y')
            title('full')
            pause(0.01)
        end
        if first_good_shot
            %only need to store this on first shot becasue the same for
            %all shots
            al_pulses.window(pulse, :, :) = pulse_win_txy;
            al_pulses.time_cen(pulse, :) = t_pulse_cen;
        end
        al_pulses.num_counts(shot, pulse) = size(counts_pulse(:, 3), 1);
        al_pulses.pos.mean(shot, pulse, :) = mean(counts_pulse, 1);
        al_pulses.pos.std(shot, pulse, :) = std(counts_pulse, 1);
        %TODO: convert to velocity here and then take the mean,std
        vzxy_out = txy_to_vel(counts_pulse, ...
            t_pulse_cen-opts.c.tof, ...
            -opts.c.g0, ...
            opts.c.fall_distance);
        al_pulses.vel_zxy.mean(shot, pulse, :) = mean(vzxy_out, 1);
        al_pulses.vel_zxy.std(shot, pulse, :) = std(vzxy_out, 1);
        al_pulses.vel_zxy.counts{shot,pulse} = vzxy_out;

    end %pulse
    if first_good_shot, first_good_shot = false; end
    %check that the calculated t0 is close enough

    tol_t0_match = 1e-2; %in factors of bin size
    tol_t0_match = tol_t0_match * opts.pulsedt;
    mean_cen_time = mean(al_pulses.pos.mean(shot, :, 1)-al_pulses.time_cen');
    if abs(mean_cen_time) > tol_t0_match
        est_t0 = opts.t0 + mean_cen_time;
        warning('pulses are not centered in time pehaps t0 should be %.5f', est_t0)
    end

    if mod(shot, 10) == 0, fprintf('\b\b\b\b%04u', shot), end
    %to set the pulse t0 right it can be handy to uncomment the next line
    %
end %shots

%check if the average difference between the mean count position (vs the expected center) over all shots is outside a tolerance
tol_t0_match = 1e-3; %in factors of bin size
tol_t0_match = tol_t0_match * opts.pulsedt;
mean_cen_time = nanmean(arrayfun(@(shotnum) mean(al_pulses.pos.mean(shotnum, :, 1)-al_pulses.time_cen'), 1:size(al_pulses.pos.mean, 1)));
if abs(mean_cen_time) > tol_t0_match
    est_t0 = opts.t0 + mean_cen_time;
    fprintf('t0 should be %.6f', est_t0)
end


fprintf('...Done\n')


end