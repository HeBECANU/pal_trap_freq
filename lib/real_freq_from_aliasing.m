function f_real=real_freq_from_aliasing(app_freq,samp_freq,gradient)
gradient=round(gradient);
f_real=samp_freq*abs(gradient)-sign(gradient)*app_freq;

end

