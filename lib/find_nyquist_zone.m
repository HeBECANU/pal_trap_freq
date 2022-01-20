function zone = find_nyquist_zone(f_act,f_samp)
    zone= ceil(f_act./(.5*f_samp)); %technically offset by 1
end