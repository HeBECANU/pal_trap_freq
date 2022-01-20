function f_app = find_apparent_freq_deriv(f_act,f_samp)
    nyq_zone=find_nyquist_zone(f_act,f_samp);
    zone_even = mod(nyq_zone,2) == 0;
    f_app = 0*f_samp;
    f_app(zone_even) =   (nyq_zone(zone_even)).*(1/2);
    f_app(~zone_even) =  -(nyq_zone(~zone_even)-1).*(1/2);
end
