function f_app = find_apparent_freq(f_act,f_samp)
    nyq_zone=find_nyquist_zone(f_act,f_samp);
    f_app = ((-1).^is_even(nyq_zone)).*(f_act -  (nyq_zone-is_odd(nyq_zone)).*f_samp/2);
end




%    grad_zone=find_nyquist_grad_zone(f_act,f_samp);
%    f_app = abs(2*grad_zone.*f_samp/2 -f_act );
%end


%some other functions that produce the same value

%     f_app = (-1).^(nyq_zone).*(...
%                 (  nyq_zone+  ((-1+(-1).^nyq_zone)/2)  ).*f_samp/2 ...
%                 -f_act ...
%                 );

% function f_app = nyquist_sawtooth_benchmark(f_act,f_samp)
%     nyq_zone=find_nyquist_zone(f_act,f_samp);
%     zone_even = mod(nyq_zone,2) == 0;
%     f_app = 0*f_samp;
%     f_app(zone_even) =  -f_act + (nyq_zone(zone_even)).*f_samp(zone_even)/2;
%     f_app(~zone_even) = f_act -  (nyq_zone(~zone_even)-1).*f_samp(~zone_even)/2;
% end

