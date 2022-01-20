%%
f_act = linspace(0.1,10,1e3);
f_samp = 2.5;
f_app = find_apparent_freq(f_act,f_samp);
f_pos_branch=freq_pos_branch(f_act,f_samp);


stfig('sawtooth example');
clf
plot(f_act,f_app,'Color',[100,100,100]./255,'LineWidth',1.5)
hold on
plot(f_act,f_pos_branch,'g','LineWidth',1.5)
plot(f_act,-f_pos_branch,'b','LineWidth',1.5)
plot(f_act,0.5.*f_samp,'r--','LineWidth',1.4)
yline( f_samp/2,'r--' )
yline( -f_samp/2,'r--' )
hold off
%xlim([0 2.5])
%ylim([0 1.1])
xlabel('signal freq')
ylabel('Apparent Frequency, $f_a$/ True Frequency, $f$')
%legend('Apparent Frequency','Nyquist Frequency','Location','northwest')
set(gca,'FontSize',15)
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')


%%
f_act = 2;
f_samp = linspace(0.1,10,1e3);
f_app = find_apparent_freq(f_act,f_samp);
f_pos_branch=freq_pos_branch(f_act,f_samp);


stfig('sawtooth example');
clf
plot(f_samp,f_app,'Color',[100,100,100]./255,'LineWidth',1.5)
hold on
plot(f_samp,f_pos_branch,'g','LineWidth',1.5)
plot(f_samp,-f_pos_branch,'b','LineWidth',1.5)
plot(f_samp,0.5.*f_samp,'r--','LineWidth',1.4)
plot(f_samp,-0.5.*f_samp,'r--','LineWidth',1.4)
hold off
%xlim([0 2.5])
ylim([-1.2*f_act 1.2*f_act])
xlabel('signal freq')
ylabel('Apparent Frequency, $f_a$/ True Frequency, $f$')
%legend('Apparent Frequency','Nyquist Frequency','Location','northwest')
set(gca,'FontSize',15)
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')


%%
function f_app = freq_pos_branch(f_act,f_samp)
    f_app = (f_samp/2).*sawtooth(2*pi*(f_act./(f_samp*2) +0.25),1/2);
end
