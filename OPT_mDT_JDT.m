%空时降维3dt方法
%2012-2-27

clc;clear all;close all;
load clutter_matrix.mat;
j = sqrt(-1);

[NK,L]=size(clutter_matrix);
N=16;
K=10;
CNR=60;
Rc=clutter_matrix*clutter_matrix'/L;
noise=max(max(Rc))/(10^(CNR/10))*eye(N*K);
Rx=Rc+noise;
anoise=max(max(noise));                          %噪声功率
psi0=pi/2; % 目标的锥角余弦位置为0

% 全维STAP
fd=-1:1/50:1;
inv_Rx=pinv(Rx);
for i=1:length(fd)
    Ss=exp(j*pi*(0:N-1)'*cos(psi0));            %目标方向确定时，Ss固定，但目标doppler频率未知，故每一个fd有一个最优权矢量wopt
    St=exp(j*pi*(0:K-1)'*fd(i));
    S=kron(St,Ss); % 直接用该对锥角余弦和归一化多普勒频率生成导向矢量
    wopt=inv_Rx*S/(S'*inv_Rx*S);
    IF(i)=abs(wopt'*S)^2*(10^(CNR/10)+1)*anoise/(wopt'*Rx*wopt);
end
figure
plot(fd,10*log10(abs(IF)))
xlabel('2f_d/f_r');ylabel('IF/dB');
grid on

% 3DT降维STAP
index=0;
for fd=-1:1/50:1
    index=index+1;
    Qt=[exp(j*pi*(0:K-1).'*(fd - 1/K)),exp(j*pi*(0:K-1).'*fd),exp(j*pi*(0:K-1).'*(fd + 1/K))];%时域降维矩阵
    Qs=eye(N); % 空域降维矩阵
    Q=kron(Qt,Qs);
    Ry=Q'*Rx*Q; % 数据降维之后的协方差矩阵，数据降维后被修正成了48*48
    inv_Ry=pinv(Ry);

    Ss=exp(j*pi*(0:N-1).'*cos(psi0));
    St=exp(j*pi*(0:K-1).'*fd);
    S = kron(St,Ss);

    Sy = Q'*S;
    W_3dt=inv_Ry*Sy/(Sy'*inv_Ry*Sy); % 数据矩阵需要降维，使用的导向矩阵需要进行修正，要增加通道的归一化增益
    IF_3dt(index)=abs(W_3dt'*Sy).^2*(10^(CNR/10)+1)*anoise/(W_3dt'*Ry*W_3dt);
end
hold on
plot(-1:1/50:1,10*log10(IF_3dt),'.-');

% JDL降维STAP
index = 0;
for fd=-1:1/50:1
    index = index + 1;
    % 进行数据降维，筛选目标附近的九宫格数据
    Qt = [exp(j*pi*(0:K-1).'*(fd-1/K)),exp(j*pi*(0:K-1).'*fd),exp(j*pi*(0:K-1).'*(fd+1/K))];%时域降维矩阵
    Qs = [exp(j*pi*(0:N-1).'*(cos(psi0)-1/N)),exp(j*pi*(0:N-1).'*cos(psi0)),exp(j*pi*(0:N-1).'*(cos(psi0)+1/N))];%空域降维矩阵

    Q = kron(Qt,Qs);
    Rz = Q'*Rx*Q; % 数据降维后的协方差矩阵，被修正成了9*9
    inv_Rz = pinv(Rz);
    Ss=exp(j*pi*(0:N-1).'*cos(psi0));
    St=exp(j*pi*(0:K-1).'*fd);
    S = kron(St,Ss);

    % 求解新的空间导向矢量，空间导向矢量维度为9
    Sz = Q'*S;
    w_jdl = inv_Rz*Sz/(Sz'*inv_Rz*Sz); 
    IF_jdl(index)=abs(w_jdl'*Sz).^2*(10^(CNR/10)+1)*anoise/(w_jdl'*Rz*w_jdl);
end
hold on
plot(-1:1/50:1,10*log10(IF_jdl),'.-');
legend('最优','3DT','JDL')