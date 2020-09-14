%try to make aliasing spectrogram


addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path

%%
osc_function=@(x) (3*sin(x*2*pi*417+1)).^1.4+1.6;

samp_freqs=col_vec(linspace(10,1000,1e3));
samp_intervals=1./samp_freqs;
iimax=numel(samp_intervals);
num_samps=repmat(1000,iimax,1);
sample_structures=cell(iimax,1);

spectrogram=[];
spectrogram.fmin=0.5;
spectrogram.fmax=500;
spectrogram.num_f_steps=100000;
spectrogram.f_sampling=1./samp_intervals;
spectrogram.f_response=col_vec(linspace(spectrogram.fmin,spectrogram.fmax,spectrogram.num_f_steps));
spectrogram.slice=nan(numel(spectrogram.f_sampling),numel(spectrogram.f_response));

%%
for ii=1:iimax
   samp_tmp.t=(0:(num_samps-1))*samp_intervals(ii);
   samp_tmp.x=osc_function(samp_tmp.t);
   freq_amp=fft_tx(samp_tmp.t,samp_tmp.x,'padding',10);
   spectral_amp=abs(freq_amp(2,:));
   
   interpolated_slice=interp1(freq_amp(1,:),spectral_amp,spectrogram.f_response,'spline',nan);
   spectrogram.slice(ii,:)=interpolated_slice;
   %%
%    clf
%    plot(spectrogram.f_response,interpolated_slice)
%    xlim([min(spectrogram.f_response),max(spectrogram.f_response)])
%    ylim([0,3])

%    pause(0.001)
end


%%
 
spectrogram_matrix=spectrogram.slice';
spectrogram_matrix(spectrogram_matrix<0)=0;
spectrogram_matrix=spectrogram_matrix.^0.3;
%spectrogram_matrix=imgaussfilt(spectrogram_matrix,0);
imagesc(spectrogram.f_sampling,spectrogram.f_response,spectrogram_matrix)
xlabel('sampling frequency (Hz)')
ylabel('observed frequency (Hz)')
colormap(viridis) %plasma

ax=gca;
ax.YDir='normal';