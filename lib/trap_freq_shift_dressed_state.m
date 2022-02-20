

%addpath('../lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
%set_up_project_path('..')
hebec_constants

%% http://www-lpl.univ-paris13.fr/bec/bec/Teaching/lecture2_2013.pdf
% https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/will_pdf_15083.pdf page 27
% 

oc_time_interval=8e-3;
oc_pulse_time=5e-6;
oc_fraction=0.1e-3;
rabi_pop=@(t,rabi_afreq) sin((1/2)*t*rabi_afreq)^2 ;
oc_rabi_afreq=2*asin(sqrt(oc_fraction))/oc_pulse_time;
%=fzero(@(rabi_afreq) rabi_pop(oc_pulse_time,rabi_afreq) -oc_fraction,100 );
oc_rabi_freq=oc_rabi_afreq/(2*pi);
fprintf('OC rabi freq %.3f kHz \n',oc_rabi_freq*1e-3)
fprintf('expected OC fraction %.3f \n',rabi_pop(oc_pulse_time,oc_rabi_afreq))

%%

trap_afreq=410*2*pi;
x_osc_vmax=10e-3;
x_osc_max=x_osc_vmax/trap_afreq;
trap_bias_afreq=1.0e6*2*pi; %in hz
% find the potneital of the dressed trap
%http://info.phys.unm.edu/~ideutsch/Classes/Phys566F13/Notes%20from%20others/Atom-Light%20Interactions.pdf


hb=const.hb;
m=const.mhe;

omega_trap=trap_afreq;

% if we optimized the drive to get the best outcoupling this would be it
mean_pot_shift_under_osc=(1/4)*x_osc_max^2 *m* omega_trap^2 ;
mean_freq_shift_under_osc=mean_pot_shift_under_osc/hb;
drive_afreq=trap_bias_afreq+ mean_freq_shift_under_osc-0.0e3*2*pi;
fprintf('Drive Freq from trap bottom %.3f kHz \n',1e-3*(drive_afreq-trap_bias_afreq)/(2*pi))
%
pot_e_fun=@(x) hb*trap_bias_afreq+ (1/2) *m *omega_trap^2 * x.^2 ;
detune_fun=@(x) drive_afreq  - pot_e_fun(x)/hb ;
gen_rabi_freq_fun=@(x) sqrt(oc_rabi_afreq.^2 + detune_fun(x).^2);
dressed_pot_e_fun=@(x) (hb/2)*( -detune_fun(x) -gen_rabi_freq_fun(x)  )-drive_afreq*hb;
dressed_pot_g_fun=@(x) (hb/2)*( -detune_fun(x) + gen_rabi_freq_fun(x) );

dressed_e_mag=@(x,mag_g,mag_e) (sqrt(mag_g)*cos((1/2)*acos(-1*detune_fun(x)./gen_rabi_freq_fun(x) )) ...
                               -sqrt(mag_e)*sin((1/2)*acos(-1*detune_fun(x)./gen_rabi_freq_fun(x) )) ).^2;

dressed_g_mag=@(x,mag_g,mag_e) (sqrt(mag_g)*sin((1/2)*acos(-1*detune_fun(x)./gen_rabi_freq_fun(x) )) ...
                               -sqrt(mag_e)*cos((1/2)*acos(-1*detune_fun(x)./gen_rabi_freq_fun(x) )) ).^2;

x_plot_max=x_osc_max*2;

x_samp=linspace(-x_plot_max,x_plot_max,1e3);
dressed_pot_e_samp=dressed_pot_e_fun(x_samp);
dressed_pot_e_samp_subtract=dressed_pot_e_samp-dressed_pot_e_fun(0);
dressed_pot_g_samp=dressed_pot_g_fun(x_samp);
dressed_pot_g_samp_subtract=dressed_pot_g_samp-dressed_pot_g_fun(0);
undressed_pot_e_samp=pot_e_fun(x_samp);
undressed_pot_e_samp_subtract=undressed_pot_e_samp-pot_e_fun(0);
undressed_pot_g_samp=x_samp*0;

% mean graident and potential
dressed_g_mag_sample=dressed_g_mag(x_samp,0,1);
dressed_e_mag_sample=dressed_e_mag(x_samp,0,1);

grad_dressed_pot_e_samp=derivest(dressed_pot_e_fun,x_samp,'MaxStep',mean(diff(x_samp)));
grad_dressed_pot_g_samp=derivest(dressed_pot_g_fun,x_samp,'MaxStep',mean(diff(x_samp)));
grad_undressed_pot_e_samp=derivest(pot_e_fun,x_samp,'MaxStep',mean(diff(x_samp)));

grad_samp_weighted= dressed_g_mag_sample.*grad_dressed_pot_g_samp ...
                    +dressed_e_mag_sample.*grad_dressed_pot_e_samp;

dressed_pot_grad_weighted=cumtrapz(x_samp,grad_samp_weighted);
dressed_pot_grad_weighted_subtract=dressed_pot_grad_weighted-interp1(x_samp,dressed_pot_grad_weighted,0);

dressed_pot_val_weighted_subtract =dressed_g_mag_sample.*dressed_pot_e_samp_subtract ...
                    +dressed_e_mag_sample.*dressed_pot_g_samp_subtract;

curv_samp_weighted=gradient(grad_samp_weighted)./gradient(x_samp);
curv_cen_weighed=interp1(x_samp,curv_samp_weighted,0);
trap_afreq_cen_weighted=sqrt(curv_cen_weighed/m);
fprintf('proj weighted dressed state trap freq %.3f Hz \n',trap_afreq_cen_weighted/(2*pi))

time_w_trap_afreq=(trap_afreq_cen_weighted*(oc_pulse_time/oc_time_interval) +...
                    trap_afreq*(1-oc_pulse_time/oc_time_interval));
time_w_trap_afreq_delta=time_w_trap_afreq-trap_afreq;

fprintf('mean over time trap freq %.3f Hz \n',time_w_trap_afreq/(2*pi))
fprintf('mean over time delta trap freq %.3f Hz \n',time_w_trap_afreq_delta/(2*pi))


time_w2_trap_afreq=sqrt((trap_afreq_cen_weighted^2*(oc_pulse_time/oc_time_interval) +...
                    trap_afreq^2*(1-oc_pulse_time/oc_time_interval)));
time_w2_trap_afreq_delta=time_w2_trap_afreq-trap_afreq;
fprintf('mean over time w^2 weight trap freq %.3f Hz \n',time_w2_trap_afreq/(2*pi))
fprintf('mean over time  w^2 weight delta trap freq %.2f mHz \n',1e3*time_w2_trap_afreq_delta/(2*pi))

%%

font_name='cmr10';
linewidth=1.5;
font_size=12;
colors_main=[[[214,72,154];
[102,181,69];
[151,90,214];
[208,153,44];
[212,73,58]]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2);
colors_shaded=colorspace('HSV->RGB',hsv);

yscale=1e-3/hb;
xscale=1e6;
stfig('dressed potnetial');
clf
%plot(x_samp,undressed_pot_e_samp_subtract)
hold on
plot(x_samp*xscale,dressed_pot_e_samp_subtract*yscale,'-','LineWidth',linewidth,'Color',colors_main(1,:))
plot(x_samp*xscale,dressed_pot_g_samp_subtract*yscale,'-','LineWidth',linewidth,'Color',colors_main(3,:))
plot(x_samp*xscale,dressed_pot_grad_weighted_subtract*yscale,':','LineWidth',linewidth,'Color',colors_main(4,:))
% just seeing how different this one is
%plot(x_samp*xscale,dressed_pot_val_weighted_subtract*yscale,'-','LineWidth',linewidth,'Color',colors_main(4,:))
plot(x_samp*xscale,undressed_pot_e_samp_subtract*yscale,'--','LineWidth',linewidth,'Color',colors_main(1,:))
yline(0,'--','LineWidth',linewidth,'Color',colors_main(3,:))
xline([-1,1]*x_osc_max*xscale,'-.','LineWidth',linewidth,'Color',colors_main(2,:))
hold off
legend('$\left|1^\prime\right\rangle$',...
    '$\left|0^\prime\right\rangle$',...
    'D. Proj. Weighted',...
    '$\left|1\right\rangle$',...
    '$\left|0\right\rangle$',...
    'Osc. Max')
legend('location','north')
ylabel('Potential, $(U(x)-U(0))/\hbar$ (kHz)')
xlabel('Position ($\mu$m)')

yl=ylim;
ylim([yl(1)-diff(yl)*0.05,yl(2)])
box on
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.015,0])

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=600;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[100,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off
%%
fig_name='dressed_state_potentials';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))

%%


stfig('dressed projection');
clf
subplot(1,2,1)
hold on

plot(x_samp*xscale,dressed_e_mag_sample,'-','LineWidth',linewidth,'Color',colors_main(1,:))
plot(x_samp*xscale,dressed_g_mag_sample,'--','LineWidth',linewidth,'Color',colors_main(3,:))
xline([-1,1]*x_osc_max*xscale,'-.','LineWidth',linewidth,'Color',colors_main(2,:))
hold off
legend('$\left|1^\prime\right\rangle$','$\left|0^\prime\right\rangle$','Osc. Max')
legend('location','north')
ylabel('Projection Magnitude')
xlabel('Position ($\mu$m)')
xlim([x_samp(1),x_samp(end)]*xscale)
box on
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.015,0])
sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.01, sp_pos(2)+sp_pos(4)*0.99, 0, 0], 'string', '(A)','FontSize',round(font_size*1.2),'FontName',font_name);



yscale=1e-3*1e-6/hb;
xscale=1e6;
subplot(1,2,2)
hold on
plot(x_samp*xscale,grad_dressed_pot_e_samp*yscale,'-','LineWidth',linewidth,'Color',colors_main(1,:))
plot(x_samp*xscale,grad_dressed_pot_g_samp*yscale,'-','LineWidth',linewidth,'Color',colors_main(3,:))
plot(x_samp*xscale,grad_samp_weighted*yscale,':','LineWidth',linewidth,'Color',colors_main(4,:))
plot(x_samp*xscale,grad_undressed_pot_e_samp*yscale,'--','LineWidth',linewidth,'Color',colors_main(1,:))
%plot(x_samp*xscale,x_samp*m*trap_afreq_cen_weighted^2*yscale,'-.','LineWidth',linewidth,'Color',colors_main(5,:))
xline([-1,1]*x_osc_max*xscale,'-.','LineWidth',linewidth,'Color',colors_main(2,:))
hold off
legend('$\left|1^\prime\right\rangle$',...
    '$\left|0^\prime\right\rangle$',...
        'D. Proj. Weighted',...
        '$\left|1\right\rangle$',...
        'Osc. Max')
legend('location','best')
ylabel('Pot. Grad., $\frac{d}{dx} U(x)/\hbar$ (kHz/$\mu$m)')
xlabel('Position ($\mu$m)')
xlim([x_samp(1),x_samp(end)]*xscale)
ylim(minmax([minmax(grad_dressed_pot_g_samp),minmax(grad_dressed_pot_e_samp)])*yscale)
box on
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.015,0])
sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.01, sp_pos(2)+sp_pos(4)*0.99, 0, 0], 'string', '(B)','FontSize',round(font_size*1.2),'FontName',font_name);


set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=1000;
fig_aspect_ratio=0.3; %0.67;
set(gcf,'Position',[700,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off

%%
fig_name='dressed_pot_fraction_and_grad';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))



