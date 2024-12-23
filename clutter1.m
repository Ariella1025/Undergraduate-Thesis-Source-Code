    %
% 硕士期间的杂波仿真，准确地说是杂波+干扰
%先产生噪声背景下的干扰，再考虑机载雷达的情况，将杂波考虑进去.主瓣方向为theta=90°，phi=1°
%  统一时域卷积空域
clc;
clear all;
close all;
%% 参数设置
%-- --系统参数------
Pt=230e3;           %峰值发射功率
lamda=0.23;         %工作波长
fr=2434.8;          %重频
B=5e6;              %接收机带宽
F=3;                %噪声系数dB
tao=15e-6;          %脉冲宽度
K=16;               %相参脉冲数
sigma2=1;           %噪声功率   w
k1=1.38e-23;     % 
T=290;              %温度
c=3e8;              %光速
Ls=10^0.4;          %接收机损耗
%-----天线参数-----
%阵元采用余弦方向图
Gel=4;              %阵元增益dB
M=4;                %列向阵元数，
N=16;               %子阵数
d=lamda/2;          %阵元间距
%------载机参数----
H=8000;             %载机高度
v=140;              %载机速度
%------主瓣指向----
theta=90;
phi=1;              %目标俯视角，方位角为theta，锥角为psi   主瓣方向俯仰角1°   列向合成存在相位误差，需要补偿



%------杂波参数----
gamma=-3;                     %归一化散射系数dB
theta_c=(0:359)*(pi/180);     %杂波块个数为360
N_RangeSample = 2000;          %距离采样数   看似是人为计算，实际上是计算出来的，一个不模糊区间的距离门个数Ru/deltaR
Nc=360;
deltaTheta = 2*pi/Nc;          %方位分辨率/rad

open_phasecompensation_c=exp(-1i*2*pi*d/lamda*sin(1*pi/180)*(0:M-1));    %杂波（俯仰）相位补偿算子
open_phasecompensation_c=open_phasecompensation_c/norm(open_phasecompensation_c);
Mtx_synthesisoperator_c=kron(eye(N),open_phasecompensation_c);               %杂波俯仰合成算子

%% ------------目标生成------------
ksi_t=10^(3/10);                                              %信噪比3db     
a=exp(1i*2*pi*(0:N-1)'*0);                                   
b=exp(1i*2*pi*(0:K-1)'*0.3);                                   %2*vt/(lamda*fr)
Xt=sqrt(sigma2*ksi_t/2)*kron(b,a);


%% -----------干扰生成--------------
Sj=1e-3;                 %干扰机的有效辐射功率密度（W/Hz）  1000
R_j=150e3;              %干扰源所在斜距,考虑干扰源在地面   150km
phi_j=asin(H/R_j);      %干扰源的俯仰角
Fr_j=sum(exp(1i*2*pi*d/lamda*(0:M-1).'*(sin(phi_j)-sin(1*pi/180))));       %干扰的接收方向图  .*Fe_j?
Fr_j=Fr_j/max(max(Fr_j));
Gr_j=M*(10^(Gel/10))*(abs(Fr_j)).^2;                                       %干扰的接收增益
S_j=exp(1i*2*pi*(0:N-1)'*0.05);                                            %16*1
J0=Sj*Gr_j*lamda^2/((4*pi)^2*R_j^2*Ls);           %干扰功率
ksi_j=J0/(k1*T*10^(F/10));                        %干噪比   分母原本乘了一个B，陈师兄
%% ------------杂波生成----------------
Ru=c/(2*fr);                    %最大不模糊距离约为61.6公里
d0=4.12*(sqrt(H)+0)*1000;       %雷达视距 368.5km
fs = 2.5e6;
deltaR = c/(2*B);              % 距离分辨率m     30m
Re = 8.49e6;                   %地球曲率半径km
Ft_c= zeros(1,Nc);             %发射方向图
Fr_c= zeros(1,Nc);             %接收方向图
Xcjn=zeros(N*K,N_RangeSample); %杂波+干扰接收数据
Xcn=zeros(N*K,N_RangeSample);  %杂波接收数据
S_c=zeros(N*K,Nc);
Gt_c=zeros(1,Nc);              %发射增益
Gr_c=zeros(1,Nc);              %接收增益
rcjn=zeros(N*K);
rcn=zeros(N*K);
Xc=zeros(N*K,N_RangeSample);
Rc_real=zeros(N*K);
Rcjn_real=zeros(N*K);
for l = 1:N_RangeSample                             %距离采样遍历
    Rcell = (l-1)*deltaR;                           
    for m=3                                         %对于第三个模糊区间     150km
        Rc =(m-1)*Ru + Rcell;                       %不同模糊区间的杂波距离  (m-1)*Ru + Rcell
        if Rc>=H &&  Rc<= d0
            sin_psi = (2*H*Re+H^2-Rc^2)/2/Rc/Re;       %擦地角           若不严格考虑擦地角而将俯仰角视为擦地角，最终无法得出结果！！！！
            if sin_psi>0 && (1-sin_psi)>(1e-10)
                phi_c=asin(-(2*H*Re+H^2+Rc^2)/2/Rc/(Re+H));         %杂波俯仰角
                psi_c=acos(cos(theta_c)*cos(phi_c));                %计算锥角
                cos_psi = sqrt(1-sin_psi^2);                        %擦地角余弦
                Fe_c=abs(cos(theta_c+90*pi/180));                   %阵元方向图采用余弦方向图
                Ft_c=sum(exp(1i*2*pi*d/lamda*(0:N-1).'*(cos(psi_c)-cos(90*pi/180)*cos(1*pi/180))))*sum(exp(1i*2*pi*d/lamda*(0:M-1)'*(sin(phi_c)-sin(1*pi/180)))).*Fe_c;    %整个阵面总的发射方向图
                Fr_c=sum(exp(1i*2*pi*d/lamda*(0:M-1).'*(sin(phi_c)-sin(1*pi/180)))).*Fe_c;
                Ft_c=Ft_c/max(max(Ft_c));
                Fr_c=Fr_c/max(max(Fr_c));
                Gt_c=M*N*(10^(Gel/10))*(abs(Ft_c)).^2;
                Gr_c=M*(10^(Gel/10))*(abs(Fr_c)).^2;
                sigma0=10^(gamma/10)*sin_psi;                                   %散射系数,等gama模型
                % ClutterPatchArea = Rc*deltaTheta*deltaR/cos_psi;
                % sigma = sigma0*ClutterPatchArea;
                Ag=c*tao/2*1/cos_psi*Rc*sin(1/360*2*pi);                   %杂波块的地面可分辨面积 长*宽
                ksi_c=sqrt((Pt*B*tao*Gt_c.*Gr_c*lamda^2*Ag*sigma0)/((4*pi)^3*Rc^4*k1*T*B*10^(F/10)*Ls));    %杂噪比
                fd_c=2*cos(phi_c)*cos(theta_c)*v/(lamda*fr);
                St_c=exp(1i*(0:K-1)'*2*pi*fd_c);                                %时域导向矢量
                fs_cn=cos(phi_c)*cos(theta_c)*2*pi*d/lamda;
                fs_cm=sin(phi_c)*2*pi*d/lamda;
                Ss_c=(ones(N,1)*ksi_c).*(Mtx_synthesisoperator_c*kron(exp(1i*(0:N-1)'*fs_cn),exp(1i*(0:M-1)'*fs_cm)));    %空域导向矢量,面阵的导向矢量与线阵不同
                for i=1:Nc
                    S_c(:,i)=kron(St_c(:,i),Ss_c(:,i));                             %时域卷积空域     NK*360
                end
                if l == 1000
                    Rc_real = S_c*S_c';                                             %每一个距离门上真实(理论)的Rc
                end
                Xc(:,l)=S_c*(rand(1,Nc)+1i*randn(1,Nc))';                       %每一个距离采样点上的杂波数据
                Xn=sqrt(sigma2/2)*(randn(N*K,1)+1i*randn(N*K,1));               %Xn=wgn(N*K,1,0)  0db是统计意义上的值，在非大样本情况下实际功率不是0db
                % Xjt=sqrt(sigma2*ksi_j/2)*(randn(K,1)+1i*randn(K,1));          %干扰的时域数据，全多普勒域
                % if l>=1000&&l<1200
                Xj=sqrt(sigma2*ksi_j/2).*kron((randn(K,1)+1i*randn(K,1)),S_j);   % 时域卷积空域
                % else
                %     Xj=zeros(N*K,1);
                % end
                Xcjn(:,l) = Xcjn(:,l)+Xc(:,l)+Xn+Xj;                     %第l个距离单元的空时快拍数据  256*1
                Xcn(:,l) = Xcn(:,l)+Xc(:,l)+Xn;  
            end
        end
    end
end

%% 画图
Rcjn_real=Rc_real+Xj*Xj'+eye(N*K);
Rcn_real=Rc_real+eye(N*K);
for p=1:3*N*K                 %用3NK个样本来进行估计
    rcjn=rcjn+Xcjn(:,p)*Xcjn(:,p)';    %+eye(N*K)
    rcn=rcn+Xcn(:,p)*Xcn(:,p)';
end
Rcjn=rcjn/(3*N*K);                  %估计值
Rcn=rcn/(3*N*K); 

%计算杂波功率谱P并绘图
fs_num=101;
fd_num=101;
fs_n=linspace(-0.5,0.5,fs_num);
fd_n=linspace(-0.5,0.5,fd_num);
Rcjn_inv=inv(Rcjn);
P=zeros(fd_num,fs_num);
for i=1:fd_num
    St=exp(1i*(0:K-1)'*2*pi*fd_n(i));
    for j=1:fs_num
        Ss=exp(1i*(0:N-1)'*2*pi*fs_n(j));     %此处虚数原先用j，混淆了。。
        S=kron(St,Ss);                        %先时域再空域！！！！
        S=S/norm(S);                          %归一化使噪声为功率在图中显示为零
        P(j,i)=10*log10(abs(1/(S'*Rcjn_inv*S)));
    end
end

figure(1);
colormap(jet);
surf(fd_n,fs_n,P);
shading interp;lighting gouraud;colorbar;
title('杂波功率谱');xlabel('fd');ylabel('fs');zlabel('P/dB');hold on;  %绘制杂波功率谱图

%特征谱图
EigVal = eig(Rcjn);
EigVal = db(sort(abs(EigVal),'descend'))/2;
EigVal2 = eig(Rcn);
EigVal2 = db(sort(abs(EigVal2),'descend'))/2;

figure(2);
plot(EigVal','*-');
title('特征谱');
grid on;
hold on 
plot(EigVal2','+-');
legend('杂波+干扰','杂波','Location','NorthWest'); 

%距离多普勒图
%先对每个阵元接收数据做DFT，后空域合成
clutter_echo_matrix = zeros(N,K,N_RangeSample);
Vec_Qs = chebwin(N,60).';                               %
Vec_Qs = Vec_Qs/(sqrt(Vec_Qs*Vec_Qs'));                %归一化的空域
u_s = Vec_Qs;                                          %1*N的空域合成算子
Vec_Qt = chebwin(K,60).';                              %1*K
Vec_Qt = Vec_Qt/(sqrt(Vec_Qt*Vec_Qt'));
fd_t=-0.5:0.01:0.49;                                  %(-K/2:K/2-1)/K;
fdata = zeros(N,100,N_RangeSample);
fdata2 = zeros(100,N_RangeSample);

for l=1:N_RangeSample
    clutter_echo_matrix(:,:,l)=reshape(Xcjn(:,l),N,K); %NK*L矩阵变维为N*K*L
    for n = 1:N                                        %先对每个阵元数据进行DFT,得到每个阵元数据仍为1*K(频域)
        fdata(n,:,l)=(Vec_Qt.*clutter_echo_matrix(n,:,l)*exp(1i*2*pi*(0:K-1).'*fd_t));         %(1*K).*(1*K)*(K*K)第l个距离门中第n个阵元进行DFT处理以后的回波数据，exp(j*2*pi*(0:K-1).'*fd))为DFT权值
    end
    fdata2(:,l) = u_s*fdata(:,:,l);  %第l个距离门中对每个阵元数据进行DFT后再用空域合成算子进行空域合成
end

x=1:length(fd_t);                       %1:size(fdata2,1)
y=1:N_RangeSample;                       %1:size(fdata2,2)
figure(2);colormap(jet);
mesh(x,y,20*log10(abs(fdata2.')));
shading interp;
lighting gouraud;colorbar;
title('处理前距离-多普勒图');
% axis([1 length(fd_t) 1 N_RangeSample]);
xlabel('多普勒通道');ylabel('距离门');
