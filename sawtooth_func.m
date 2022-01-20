%%
f_act = 1;
f_samp = linspace(0.0,3,1e6);
f_app = find_apparent_freq(f_act,f_samp);


stfig('sawtooth example');
clf
plot(f_samp./f_act,f_app./f_act,'Color',[100,100,100]./255,'LineWidth',1.5)
hold on
plot(f_samp./f_act,0.5.*f_samp./f_act,'r--','LineWidth',1.4)
hold off
xlim([0 2.5])
ylim([0 1.1])
xlabel('Normalized Sampling Frequency, $f_s/f$')
ylabel('Normalized Apparent Frequency, $f_a/f$')

% equations 
%txt=text(0.3,0.8,{'$f_{\mathrm{Nyquist}}=\frac{1}{2}f_s$', ...
%    '$N =\lceil f/f_{\mathrm{Nyquist}} \rceil $', ...
%                   '$f_a=\left|f- f_{\mathrm{Nyquist}}(N-isodd(N))\right|$'},...
%   'Color','k','FontSize',20)

%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')


txt=text(0.05,0.65,'$N=$','Color','k','FontSize',20);
vert_pos_N_text=0.65;
f_nyq_edge_lower=2.5;
for ii=1:6
    f_nyq_edge_upper=2*f_act/(ii);
    xline(f_nyq_edge_upper,'-.','Color',[0.5,1,0.5]*0.5)
    hoz_pos=mean([f_nyq_edge_upper,f_nyq_edge_lower]);
    txt=text(hoz_pos,vert_pos_N_text,sprintf('%u',ii),'Color','k','FontSize',20,'HorizontalAlignment','center');
    f_nyq_edge_lower=f_nyq_edge_upper;
end
%set(gca, 'XScale', 'log')

ln=legend('Apparent Frequency, $f_a$',...
    'Nyquist Frequency, $f_{\mathrm{Ny}}$',...
    'Nyquist Zones',...
    'Location','northwest','FontSize',14)
lnpos=get(ln,'Position')
lnpos(1)=0.25
lnpos(2)=0.75
set(ln,'Position',lnpos)
set(ln,'EdgeColor',[1,1,1])
%legend boxoff
set(gca,'FontSize',15)


set(gcf,'Position',[100,100,800,600])
%plot_name='impulse_n_3e5_osc_10mms_in_y';
plot_name=sprintf('alias_sawtooth_th');
out_dir='./figs/bthesis/'
saveas(gcf,[out_dir,plot_name,'.png'])
%saveas(gcf,[out_dir,plot_name,'.pdf'])
saveas(gcf,[out_dir,plot_name,'.svg'])
%saveas(gcf,[out_dir,plot_name,'.fig'])
export_fig([out_dir,plot_name,'.eps'])
export_fig([out_dir,plot_name,'.pdf'])


xlabel('Sampling Frequency, $f_s$ (Hz)')
ylabel('Apparent Frequency, $f_a$ (Hz)')
plot_name=sprintf('alias_sawtooth_th_simple');
saveas(gcf,[out_dir,plot_name,'.png'])
%saveas(gcf,[out_dir,plot_name,'.pdf'])
saveas(gcf,[out_dir,plot_name,'.svg'])
%saveas(gcf,[out_dir,plot_name,'.fig'])
export_fig([out_dir,plot_name,'.eps'])
export_fig([out_dir,plot_name,'.pdf'])

%%
stfig('sawtooth deriv');
f_act = 1;
Nyq_max=15;
f_samp = linspace(2*f_act/Nyq_max,3,1e6);
f_deriv = find_apparent_freq_deriv(f_act,f_samp);

clf
plot(f_samp./f_act,f_deriv,'k','LineWidth',1.5)
xlim([0 2.5])
ylim([-7.5 7.5])
xlabel('Normalized Sampling Frequency, $f_s/f$')
ylabel('Apparent Frequency Deriv. $\partial f_a/\partial f_s$')
set(gca,'FontSize',15)


txt=text(0.4,6.5,'$N$=','Color','k','FontSize',20);
vert_pos_N_text=5.5;
txt=text(0.4,-4,'$\frac{\partial f_a}{\partial f_s}$=','FontSize',25)
vert_pos_d_text=-5.5;
f_nyq_edge_lower=2.5;
for ii=1:6
    f_nyq_edge_upper=2*f_act/(ii);
    xline(f_nyq_edge_upper,'-.','Color',[0.5,1,0.5]*0.5)
    hoz_pos=mean([f_nyq_edge_upper,f_nyq_edge_lower]);
    txt=text(hoz_pos,vert_pos_N_text,sprintf('%u',ii),'Color','k','FontSize',20,'HorizontalAlignment','center');
    f_nyq_edge_lower=f_nyq_edge_upper;
    deriv_here=find_apparent_freq_deriv(f_act,hoz_pos)
    if ii==7
        hoz_pos=hoz_pos-0.025
    end
    deriv_txt=sprintf('%d',deriv_here)
    txt=text(hoz_pos,vert_pos_d_text,deriv_txt,'Color','k','FontSize',20,'HorizontalAlignment','center');
 
end
%set(gca, 'XScale', 'log')

%%
set(gcf,'Position',[100,100,800,600])
%plot_name='impulse_n_3e5_osc_10mms_in_y';
plot_name=sprintf('alias_sawtooth_deriv_th');
out_dir='./figs/bthesis/'
saveas(gcf,[out_dir,plot_name,'.png'])
%saveas(gcf,[out_dir,plot_name,'.pdf'])
saveas(gcf,[out_dir,plot_name,'.svg'])
saveas(gcf,[out_dir,plot_name,'.fig'])
export_fig([out_dir,plot_name,'.eps'])
export_fig([out_dir,plot_name,'.pdf'])



