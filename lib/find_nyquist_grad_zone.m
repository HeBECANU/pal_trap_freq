function grad_zone = find_nyquist_grad_zone(f_act,f_samp)
    grad_zone=floor(0.5+f_act./(f_samp));
end

% used this in deriving this function
% function out=shifted_zone(f_act,f_samp)
%     nyq_zone=find_nyquist_zone(f_act,f_samp);
%     out= nyq_zone+  ((-1+(-1).^nyq_zone)/2)  ;
% end