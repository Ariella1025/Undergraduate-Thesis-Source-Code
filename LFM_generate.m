clear;clc;

%% LFM信号参数设定
% LFM参数设定
fs = 400e6;                                 % 数据率400MHz
ts = 1/fs;                                  % 数据采样间隔
sig_width_d = 1000;                         % 信号的数字时宽
sig_width_a = sig_width_d*ts;               % 信号的模拟时宽
band_width = 50e6;                          % 信号带宽50MHz
K = band_width/sig_width_a;                 % 调频斜率K
W = linspace(-pi,pi,1000);

%% LFM时域信号生成
t = -sig_width_a/2:ts:sig_width_a/2-ts;     % 信号模拟时间
LFM_sig = exp(1j*pi*K*t.^2);                % LFM信号时域序列

%% LFM信号频谱
LFM_spectrum = fft(LFM_sig,sig_width_d);
LFM_phase_spectrum = angle(LFM_spectrum);   % LFM信号幅度谱
LFM_amplitude_spectrum = abs(LFM_spectrum); % LFM信号相位谱
delta_phase = diff(LFM_phase_spectrum);     % LFM信号相位差分

%% 绘图
figure(1)
plot(t,real(LFM_sig));title('LFM信号实部');
axis([min(t) max(t) min(real(LFM_sig)) max(real(LFM_sig))])
figure(2)
plot(t,imag(LFM_sig));title('LFM信号虚部');
axis([min(t) max(t) min(imag(LFM_sig)) max(imag(LFM_sig))])

figure(3)
plot(W,fftshift(LFM_amplitude_spectrum));title('LFM信号幅度谱');
axis([-pi pi min(fftshift(LFM_amplitude_spectrum)) max(fftshift(LFM_amplitude_spectrum))])
figure(4)
plot(W,fftshift(LFM_phase_spectrum));title('LFM信号相位谱');
axis([-pi pi min(fftshift(LFM_phase_spectrum)) max(fftshift(LFM_phase_spectrum))])
% subplot(313);plot(fftshift(delta_phase));title('LFM信号相位差分');

