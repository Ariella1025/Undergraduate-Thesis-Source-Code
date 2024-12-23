%% 正侧视阵，考虑阵元误差
clc
clear all
j =sqrt(-1);
tic

%% 基础数据
M = 8; % 阵元个数
N = 8; % 阵元个数
K = 16; % 相干脉冲处理个数

lambda = 0.23; % 信号波长
d = lambda/2; % 阵元间距
c = 3e8; % 光速
Re = 849e4; % 地球半径
Lc = 10^(3/10); % 环境损耗

Va = 140; % 载机飞行速度
H = 8000; % 载机飞行高度

CNR = 60; % 杂噪比

dB = 30;
Im = chebwin(M,dB);
In = chebwin(N,dB);

Pt = 200e3; % 发射功率
G = 10^(30/10); % 天线发射增益、天线每一路阵元接收增益

phi_0 = 30/180*pi; % 天线主瓣俯仰角
theta_0 = 90/180*pi; % 天线主瓣方位角
cos_psi_0 = cos(phi_0)*cos(theta_0); % 主波束方向的锥角
alpha = 1/180*pi; % 天线半波束宽度

fr = 6000; % 脉冲重频
tau = 0.1e-5; % 脉宽
Ru = c/(2*fr); % 最大不模糊距离
Rmax = sqrt(2*Re*H); % 最大探测距离
Rmin = 2000; % 最小探测距离
R_cell = c*tau/2; % 距离环宽度

theta = [0:180]/180*pi; % 方位角划分
W = length(theta); % 方位角划分个数

L = floor((Ru)/R_cell); % 距离环划分个数

phi = zeros(1,L); % 每个距离环上的俯仰角
Rg = zeros(1,L); % 每个距离环上的地距
R = zeros(1,L); % 每个距离环的斜距

cos_psi = zeros(L,181); % 每个散射单元对应的锥角

%% 确定每个距离环上的相关数据
for l = 1:L
    Rg(l) = (l-1)*R_cell;
    R(l) = sqrt(H^2 + Rg(l)^2);
    phi(l) = asin(H/R(l));
end

% 确定天线发射方向图和接收
F_t = zeros(L,W); % 天线发射方向图，和方位角和俯仰角有关
g_n = zeros(1,L); % 天线的接受方向图，按列合成，共有N个，忽略阵元误差

for l = 1:L
    for w = 1:W
        cos_psi(l,w) = cos(phi(l))*cos(theta(w));
        F_t(l,w) = sum(In'.*exp(j*2*pi*d/lambda*[0:N-1]*(cos_psi(l,w)-cos_psi_0)))*sum(Im'.*exp(j*2*pi*d/lambda*[0:M-1]*(sin(phi(l))*sin(phi_0))));
    end
    g_n(l) = sum(Im'.*exp(j*2*pi*d/lambda*[0:M-1]*(sin(phi(l))*sin(phi_0))));
end

%% 计算每个距离单元上的相关参数
sigma_0 = zeros(1,L); % 每个距离环上的杂波单元的后向散射系数
Ac = zeros(1,L); % 每个距离环上的杂波单元的面积
sigma_c = sigma_0 .* Ac; % 每个距离环上杂波单元的后向散射面积

eta = zeros(1,L); % 每个距离环上的增益
wt = zeros(L,W);
ws = zeros(L,W);

for l = 1:L
    sigma_0(l) = 0.1*sin(phi(l)) + 8*exp(-(pi/2-phi(l))^2 / (4/180*pi)^2); % 修正等γ模型
    Ac(l) = Rg(l)*alpha*R_cell;
    sigma_c(l) = sigma_0(l) * Ac(l);
    eta(l) = sqrt(Pt*G^2*lambda^2*sigma_c(l)/((4*pi)^3)/Lc); % 每个距离环上的增益
    for w = 1:W
        ws(l,w) = 2*pi*d*cos(phi(l))*cos(theta(w))/lambda;
        wt(l,w) = 4*pi*Va*cos(phi(l))*cos(theta(w))/lambda/fr;
    end
end
flag = zeros(N*K,W);
Clutter = zeros(N*K,1);
Clutter_Matrix = zeros(N*K,L);
Rc_clutter = zeros(N*K,N*K,L); % 用于存储所有距离环的相关矩阵
Rc = zeros(N*K,N*K); % 杂波协方差矩阵

% -----杂波幅度控制----
al=wblrnd(1.5,2.5,L,W); % 杂波服从韦布尔分布

for l = 1:L % 距离循环
    for w = 1:W
        
        % -----杂波频谱控制-----
        sigma_v = Va/(2*sqrt(2*log(2)))*alpha;
        sigma_f=2*sigma_v/lambda;
        delta_f=0.6*(2*sqrt(2*log(2))*sigma_f);
        sk=delta_f*1/sqrt(2*pi*sigma_f^2)*exp(-(([1:5]-3)*delta_f).^2/(2*sigma_f^2));
        ksai=exp(1i*2*pi*rand(1,5));
        xk=ksai.*sqrt(sk);
        for k=1:K
            X(k)=sum(xk.*exp(1i*2*pi*([1:5]-3)*delta_f*k*1/fr));
        end

        for n = 1:N
            for k = 1:K
                s = exp(j*(k-1)*wt(l,w)+j*(n-1)*ws(l,w))/(R(l)^2);
                flag((k-1)*N+n,w) = g_n(l)*F_t(l,w)*al(l,w)*eta(l)*X(k)*s;
            end
        end
    end
    Clutter = sum(flag,2);
    Clutter_Matrix(:,l) = Clutter;
    Rc_clutter(:,:,l) = Clutter_Matrix(:,l)*Clutter_Matrix(:,l)';
    Rc = Rc + Rc_clutter(:,:,l);
end
%% 加入高斯噪声
noise = max(max(Rc))/10^(CNR/10)*eye(size(Rc));         %白噪声只有自相关时不为０
Rcn = Rc+noise;

%% 计算杂波功率谱，并作图
Rcn_v = pinv(Rcn);
Rc_v = pinv(Rc);

psi=acos(linspace(-1,1,361));
fd=linspace(-1/2,1/2,361);
for s=1:length(psi)
    for t=1:length(fd)
        LL=exp(j*2*pi*(0:K-1)'*fd(t));
        PP=exp(j*pi*(0:N-1)'*cos(psi(s)));
        S=kron(LL,PP);
        Q(s,t)=1/(S'*Rcn_v*S);
    end
end

% 作图
C=max(max(abs(Q)));  

z=linspace(-1,1,361);

figure(4);  
mesh(z,z,10*log10(abs(Q)/C));  
ylabel('cos(\psi)');   
xlabel('2fd/fr'); 
zlabel('P/dB'); 
axis([min(z) max(z) min(z) max(z) min(min(10*log10(abs(Q)/C))) max(max(10*log10(abs(Q)/C)))]);

%% 计算杂波特征谱，并作图
[u,s,v] = svd(Rcn);
Eig = eig(Rcn);
Eig = db(sort(abs(Eig),'descend'))/2;

% [u1,s1,v1] = svd(Rc);

% 作图
figure(5)
plot(10*log10(diag(Eig)),'*-');
plot(Eig,'*-')
% hold on
% plot(10*log10(diag(s1)),'o-')
xlabel('特征数目')
ylabel('特征值/dB')

toc










