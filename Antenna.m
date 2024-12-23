%% 天线方向图
clc
clear all
M=8;
N=8;
dB = 30;
lambda = 0.23; % 信号波长
d = lambda/2; % 阵元间距
Im = chebwin(M,dB);
In = chebwin(N,dB);
theta=meshgrid(eps:pi/180:pi);
phi=meshgrid(eps:pi/360:pi/2)';
phi_0 = 30/180*pi; % 天线主瓣俯仰角
theta_0 = 90/180*pi; % 天线主瓣方位角
cos_psi_0 = cos(phi_0)*cos(theta_0); % 主波束方向的锥角
L = length(phi);
W = length(theta);
for l = 1:L
    for w = 1:W
        cos_psi(l,w) = cos(phi(l))*cos(theta(w));
        F_t(l,w) = sum(In'.*exp(j*2*pi*d/lambda*[0:N-1]*(cos_psi(l,w)-cos_psi_0)))*sum(Im'.*exp(j*2*pi*d/lambda*[0:M-1]*(sin(phi(l))*sin(phi_0))));
    end
    g_n(l) = sum(Im'.*exp(j*2*pi*d/lambda*[0:M-1]*(sin(phi(l))*sin(phi_0))));
end



% (1)确定天线发射方向图
figure(1)
mesh(theta./pi.*180,phi./pi.*180,10*(log10(abs(F_t))));
xlabel('方位角(°)')
ylabel('俯角(°)')
zlabel('F_t(dB)')
xlim([0 180])
ylim([20 75])
title('天线阵列发射方向图')

% (2)确定天线接收方向图
figure(2)
plot(phi./pi.*180,10*(log10(abs(g_n))));
xlabel('俯角(°)')
ylabel('g_n(dB)')
xlim([20 75])
title('天线阵列接收方向图')