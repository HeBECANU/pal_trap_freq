%%
f_act = 1;
f_samp = linspace(0,3,1e6);

f_plot = nyquist_sawtooth(f_act,f_samp);
stfig('sawtooth example');
clf
plot(f_samp,f_plot,'Color',[100,100,100]./255,'LineWidth',1.5)
hold on
plot(f_samp,0.5.*f_samp,'k--','LineWidth',1.4)
xlim([0 2.5])
ylim([0 1.1])
xlabel('Sampling Frequency, $f_s$/ True Frequency, $f$')
ylabel('Apparent Frequency, $f_a$/ True Frequency, $f$')
set(gca,'FontSize',15)
function f_app = nyquist_sawtooth(f_act,f_samp)
nyquist_zone = ceil(f_act./(.5*f_samp)); %technically offset by 1
zone_even = mod(nyquist_zone,2) == 0;
f_app = 0*f_samp;
f_app(zone_even) =  (nyquist_zone(zone_even)).*f_samp(zone_even)/2-f_act;
f_app(~zone_even) = f_act -  (nyquist_zone(~zone_even)-1).*f_samp(~zone_even)/2;
end