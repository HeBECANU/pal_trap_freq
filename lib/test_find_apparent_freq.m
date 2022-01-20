%test_find_apparent_freq

f_act = 2;
f_samp = col_vec(linspace(0.5,5,1e6));
f_app1 = nyquist_sawtooth_benchmark(f_act,f_samp);
f_app2 = find_apparent_freq(f_act,f_samp);


if ~all(isalmost(f_app1,f_app2,eps*10))
    error('not eaqul to benchmark')
end


f_app_deriv_anl=find_apparent_freq_deriv(f_act,f_samp);
f_app_deriv_num=num_grad_vecs(@(x) find_apparent_freq(f_act,x),f_samp,1e-6,1);


stfig('results')
clf 
subplot(2,1,1)
hold on
plot(f_samp,f_app1)
plot(f_samp,f_app2)
hold off
legend('benchmark','find_apparent_freq','Location','best')
subplot(2,1,2)


plot(f_samp,find_nyquist_zone(f_act,f_samp))
hold on
plot(f_samp,find_nyquist_grad_zone(f_act,f_samp),'g')
plot(f_samp,abs(f_app_deriv_num),'r')
plot(f_samp,f_app_deriv_num,'c')
plot(f_samp,f_app_deriv_anl,'g')
hold off






function f_app = nyquist_sawtooth_benchmark(f_act,f_samp)
    nyq_zone=find_nyquist_zone(f_act,f_samp);
    zone_even = mod(nyq_zone,2) == 0;
    f_app = 0*f_samp;
    f_app(zone_even) =  -f_act + (nyq_zone(zone_even)).*f_samp(zone_even)/2;
    f_app(~zone_even) = f_act -  (nyq_zone(~zone_even)-1).*f_samp(~zone_even)/2;
end

