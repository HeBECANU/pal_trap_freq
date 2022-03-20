ft_vec = linspace(11.3032-2,2+12.633,1.5e4);
F_amp = [];
for ii = 1:length(ft_vec)
dt = 8e-3;
tmax = 1.5;
ft = ft_vec(ii);
t = 0:dt:tmax;
y=sin(2.*pi.*ft.*t);
S=y;
%exp(data.mcp_tdc.al_pulses.time_cen.'.*0.2964).*
% size(S)
% Sinterp=smooth(tdat,S,0.1);%interp1(tdat,S,tint,'spline');
% S = Sinterp;
% S = S(1:101);
Fs = 1./dt;   %1./(tint(2)-tint(1));%         % Sampling frequency
T = 1/Fs;             % Sampling period
L = length(S);             % Length of signal
t = (0:L-1)*T;%(1:L);%        % Time vector
Y = fft(S);%fft(S,n);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);%P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;%linspace(0,Fs*0.5,n/2+1);%
F_amp(:,ii) = P1(18:20);
end