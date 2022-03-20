%try to make aliasing spectrogram


addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path

%%
%osc_function=@(x) (3*sin(x*2*pi*417+1)).^1.4+1.6;
%osc_function=@(x) (3*sin(x*2*pi*417+1))+1.6;
% sawtooth function (has all the harmonics)
%osc_function=@(x) (3*sawtooth(x*2*pi*417+1)).^1.4+1.6;
%trianle function
osc_function=@(x) 3*sawtooth(x*2*pi*417+1,0.5);

samp_freqs=linspace(10,10000,1e4);
samp_intervals=col_vec(linspace(0.001,0.01,1e3));
%samp_intervals=1./samp_freqs;
num_samps=1000; %number of times to sample the oscillation
amplitude_function=@(x) x.^0.1;
%amplitude_function=@(x) x;


iimax=numel(samp_intervals);

spectrogram=[];
spectrogram.fmin=0.5;
spectrogram.fmax=500;
spectrogram.num_f_steps=1e3;
spectrogram.f_sampling=1./ samp_intervals;
spectrogram.f_response=col_vec(linspace(spectrogram.fmin,spectrogram.fmax,spectrogram.num_f_steps));
spectrogram.slice=nan(numel(spectrogram.f_sampling),numel(spectrogram.f_response));

%% loop over the sampling interval
for ii=1:iimax
   samp_tmp.t=(0:(num_samps-1))*samp_intervals(ii);
   samp_tmp.x=osc_function(samp_tmp.t);
   samp_tmp.x=samp_tmp.x-mean(samp_tmp.x);
   freq_amp=fft_tx(samp_tmp.t,samp_tmp.x,'padding',10,'window','hamming');
   %freq_amp= samp_tmp
   spectral_amp=abs(freq_amp(2,:));
   spectral_amp=amplitude_function(spectral_amp);
   interpolated_slice=interp1(freq_amp(1,:),spectral_amp,spectrogram.f_response,'linear',nan);
   spectrogram.slice(ii,:)=interpolated_slice;
   %%
%    clf
%    stfig('single sampling')
%    subplot(2,1,1)
%    plot(spectrogram.f_response,interpolated_slice)
%    xlim([min(spectrogram.f_response),max(spectrogram.f_response)])
%    ylim([0,3])
%    xlabel('freq (Hz)')
%    ylabel('amp (m)')
%   % pause
%    subplot(2,1,2)
%    plot(samp_tmp.t,samp_tmp.x)
%    xlabel('time (s)')
%    ylabel('pos (m)')
end


%%
 
spectrogram_matrix=spectrogram.slice';
spectrogram_matrix(spectrogram_matrix<0)=0;
%spectrogram_matrix=spectrogram_matrix.^0.1;
%spectrogram_matrix=imgaussfilt(spectrogram_matrix,0);

pcolor(spectrogram.f_sampling,spectrogram.f_response,spectrogram_matrix)
shading flat
xlabel('sampling frequency (Hz)')
ylabel('observed frequency (Hz)')
%nlcmap=nonlinear_colormap(viridis,'power',0.1);
%colormap(nlcmap) %plasma
colormap(viridis)

ax=gca;
ax.YDir='normal';