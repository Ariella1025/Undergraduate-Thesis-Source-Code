%                杂波建模     正侧阵
%  ======================================================================
clc;    close all;  clear all;
j=sqrt(-1);
alpha =0*pi/180;      %改这个 改变杂波形状 0*pi/180   45*pi/180    90*pi/180
M=8; %阵元个数
N=10;%阵元个数
K=16;                   % 采样快拍数
N2=N;
N1=M;
lamda=0.32;
% lambda=0.12;

Va=120;                 % 飞机速度120m/s
H=1000;                 % 飞机高度
Vc=0;                   % 目标方位向速度
Vr=100;                 % 目标距离向速度
prf=6000;               % PRF：脉冲重复频率
c=3*10^8;
Ru=c/(2*prf);           % 最大不模糊距离,100km
Rmax=400000;            % 雷达远地斜距 100km
Rmin=1000;              % 近地距离  飞机正下方开始
Re = 8490000;
theta0=90/180*pi;       % 天线指向的方位角
phi0=30/180*pi;         % 天线指向的俯仰角
theta=[0:180]/180*pi;   % 方位角方位（0~180°）
distance=15;            % 距离环宽度
dis = distance;
L=200;                  % 距离环个数：L>2MN
dd = lamda/2;


R_top=0.203;%0.203;
h=0.182;



S=1; %不考虑距离模糊
for i=1:L
     Rl=Rmin+(S-1)*Ru+(i-1)*dis;
     phi(i)=asin(H/Rl+(Rl^2-H^2)/(2*Rl*(Re+H)));   %每个距离环对应的俯仰角  单位 度
end

arry_col=M;                                    %将选中阵元合成 任意列数 (偶数可被N2sel整除)
arry_row=N;                                     %将选中阵元合成 任意行数 偶数 可被N1sel整除
T8=kron(eye(arry_row),kron(ones(N/arry_row,1),kron(eye(arry_col),ones(M/arry_col,1)))); %生成聚焦矩阵T


% 
% px = linspace(0,lamda/2*(N-1),N);      %等效面阵横坐标
% pz = linspace(0,lamda/2*(M-1),M);      %等效面阵纵坐标


px = linspace(-R_top,R_top,N);      %等效面阵横坐标
pz = linspace(0,h,M);             %等效面阵纵坐标

k=1;
for i=1:length(px)
    p(1,1+(i-1)*M:i*M) = px(i);
    for m=1:length(pz)
        p(2,k) = pz(m);
        k=k+1;
    end
end


py = zeros(1,N*M);
plot3(p(1,:),py,p(2,:),'ko');
hold on;
xlabel('\it x');
ylabel('\it y');
zlabel('\it z');
title('等效阵元位置')

d=2*R_top/N;

SCR=30;                %信杂比
CNR=60;                 %杂噪比
%%%%%%-----------天线方向图-------------
dB=20;
IM=chebwin(M,dB);   %天线阵列电流强度，加窗                          ??
IN=chebwin(N,dB);   %天线阵行电流强度，加窗

F_arr=zeros(N,1);                     
Frow = zeros(1,length(theta));

theta3db=1*pi/180;                          %半功率波束宽度
r=zeros(1,N*K);

%%%%%%%%%%------------杂波反射------------%%%目标所在的方位角在天线指向的位置        并假定目标位于中间距离环上 
clutter1=zeros(N*K,1);          %用于存储每一次快拍（距离环）的列矢量（杂波＋信号）　   临时存放
clutter2=zeros(N*K,1);          %用于存储每一次快拍（距离环）的列矢量（杂波）
signal=zeros(N*K,1);            %目标信号矩阵
clutter1_matrix=zeros(N*K,L);   %用于存储所有距离环的杂波＋信号 矩阵
clutter2_matrix=zeros(N*K,L);   %用于存储除去动目标的所有距离环杂波矩阵 （包含动目标所在距离环的杂波）
Rcs_clutter_signal=zeros(N*K,N*K,L);    %三维矩阵 用于存储每一个距离环一次快拍的相关矩阵（未平均）（包含信号）
Rc_clutter=zeros(N*K,N*K,L);            %三维矩阵 用于存储每一个距离环一次快拍的相关矩阵（未平均）（不包含信号）
Rcs=zeros(N*K);                         %包含动目标距离环估计的协方差矩阵
Rc=zeros(N*K);                          %包含动目标距离环但不包含动目标信号估计的协方差矩阵 （杂波＋噪声 协方差矩阵）
for i=1:L
    i;
    Vr_temp = Vr;
    
    %% 海情
    ss=5; %海情等级
    A=4*10e-7*10^(0.6*(ss+1));
    B= pi/2;
    beita0 = (2.44*(ss+1)^1.08)/57.29;
    he = 0.025+0.046*ss^1.72;
    yitac = asin(lamda/(4*pi*he));
    k = 1.9;
    yita = phi(i);
    if yita<yitac
        seigemac0 = (yita/yitac)^k;
    else
        seigemac0 = 1;
	end
    seigema0 = A*seigemac0*sin(yita)/lamda+cot(beita0)^2*exp(-tan(B-yita)^2/tan(beita0)^2);
    RR = Rmin+(S-1)*Ru+(i-1)*dis;
    sr=RR*theta3db*dis/sqrt(2);
    seigema0 = seigema0*sr;
    
    
    
    B_coe=seigema0*exp(1j*2*pi*randn(1,length(theta)));        %第i个距离环上的杂波反射系数,幅度服从正态分布，相位服从正态分布
    B_coe_target=exp(1j*2*pi*randn(1,1));
    
    uu = cos(phi(i))*sin(theta);
    vv = cos(phi(i))*cos(theta);
    hh = sin(phi(i));
     
    
    % 面阵按列合成接收方向图过程
    w4 = repmat(exp(1j*2*pi*dd/lamda*(sin(phi(i))-sin(phi0)).*p(2,:)'),1,length(theta));
       
    F_arr = T8'*(w4);
    
    
    % 面阵  计算发射方向图
    for m=1:length(theta)
        ww = exp(1j*2*pi*dd/lamda*(sin(phi(i))-sin(phi0)).*p(2,:)').*exp(1j*2*pi*dd/lamda*(cos(theta(m))*cos(phi(i))-cos(theta0)*cos(phi0)).*p(1,:)');
        Frow(1,m) = sum(ww);
    end
    
    for k=1:K        %每个距离环的快拍信号（脉冲个数）
        for n=1:N    %阵元信号
            a=4*pi*Va/lamda/prf;
            Wt=4*pi*Va/lamda/prf*cos(theta +alpha)*cos(phi(i));        %第k次快拍时的时域角频率   i距离环的多普勒角频率        
            Ws=4*pi*d/lamda*cos(theta)*cos(phi(i));                    %第n个天线的空域角频率
            s=exp(j*(k-1)*Wt+j*(n-1)*Ws)/(Rmin+(i-1)*dis)^2;                        %信号强度衰减 1/R^2
            clutter1((k-1)*N+n,1)=sum(F_arr(n).*Frow.*B_coe.*s);       %ff为能量?　f^2*B 为回波信号ｓ（ｔ）
            clutter2((k-1)*N+n,1)=sum(F_arr(n).*Frow.*B_coe.*s);
            r(1,(k-1)*N+n)=sum(F_arr(n).*(Frow).*exp(j*((n-1)*Ws+(k-1)*Wt))/(Rmin+(i-1)*dis)^2); 
            %    -----     加入信号    ------
            if i==L/2   
                Vr_temp = Vr;
                for mm = 1: 2
                    Vr_temp = Vr_temp+100;
                    Wt=4*pi*(Va-Vc)/lamda/prf*cos(theta0+alpha)*cos(phi(i))+4*pi*Vr_temp/lamda/prf*  sin(theta0+alpha)*cos(phi(i));      %目标的时域角频率   两个速度VcVr
                    Ws=4*pi*d/lamda*cos(theta0)*cos(phi(i));                                                               %目标的空域角频率
                    signal((k-1)*N+n,1)=sum(sqrt(10^(SCR/10))*F_arr(n).*Frow.*B_coe_target.*exp(j*(k-1)*Wt+j*(n-1)*Ws))/(Rmin+(i-1)*dis)^2;%(k-1)*N+n 第k-1个快拍 的第n个空间导向信号
                    clutter1((k-1)*N+n,1)=clutter1((k-1)*N+n,1)+signal((k-1)*N+n,1);
                end
            end
            
        end
    end
    clutter1_matrix(:,i)=clutter1;
    clutter2_matrix(:,i)=clutter2;
    Rcs_clutter_signal(:,:,i)=clutter1_matrix(:,i)*clutter1_matrix(:,i)';    %自相关矩阵?　ｔ＝ｋ　ｓ关于脉冲ｋ　（时间）的自相关矩阵
    Rc_clutter(:,:,i)=clutter2_matrix(:,i)*clutter2_matrix(:,i)';
end


%%%%%%%%估计在距离环L/2处杂波协方差矩阵（动目标位于L/2距离环处）%%%%%%%%%%%%%%%%
for i=1:L         %用2NK个样本估计协方差矩阵
    Rcs=Rcs+Rcs_clutter_signal(:,:,i);
end
Rcs=Rcs/L;             %Rcs协方差矩阵中包含目标信号

noise=max(max(Rcs))/10^(CNR/10)*eye(size(Rcs));         %白噪声只有自相关时不为０

Rcs=Rcs+noise;

Rcs_v=pinv(Rcs);

%%%% capon成像
%%%%--杂波功率谱，处理器二维频响图，杂波加信号的功率谱和频响图%其实就是二维搜索，取一个锥角，然后在此锥角固定时搜索所有的多普勒频率%%%%%%%%%%
psi=acos(linspace(-1,1,361));
fd=linspace(-1/2,1/2,361);
for s=1:length(psi)
    for t=1:length(fd)
        LL=exp(j*2*pi*(0:K-1)'*fd(t));
        PP=exp(j*pi*(0:N-1)'*cos(psi(s)));
        S=kron(LL,PP);
%         Q1(s,t)=1/(S'*Rc_v*S);   %蓝
        temp = Rcs_v*S;
        Q3(s,t)=1/(S'*Rcs_v*S);  %蓝
%         Q3(s,t)=S'*Rcs*S;  %蓝
%         W = Rcs_v * S;
%         y = W*clutter1_matrix;
    end
end


% %% CST
% % I=chebwin(N,40);
% I = ones(N,1);
% Ws0=I/(I'*I);
% % Ss_fai0=exp(j*2*pi/lambda*d*[0:P-1]'*cos(fai0)*cos(theta0));
% WSs_fai0=Ws0.*exp(1j*2*pi*[0:N-1]'/lamda*sin(phi0));
% % for k=1:K
% %     St_fdj=exp(j*2*pi/K*(k-1)*[1:K]');
% %     WSt_fdj=St_fdj;
% %     WSts=kron(WSt_fdj,WSs_fai0);
% %     for m=1:L
% %         X3=clutter1_matrix(:,m);
% %         Y1(m,k)=WSts'*X3;
% %     end
% % 
% % end
% St_fdj=exp(j*2*pi/K*(k-1)*[1:K]');
% 
% for k = 1:N*K
%     St_fdj=exp(1j*2*pi/K*(k-1)*[1:K]');
%     for m = 1:L
%         y_CST(m,k)=St_fdj'*clutter1_matrix(:,m);
%     end
% end
% A=max(max(abs(Q1))); 
C=max(max(abs(Q3)));  

x=cos(psi);
y=4*Va/lamda*x/prf;
z=linspace(-1,1,361);

% figure();  mesh(x,y,10*log10(abs(Q1)/A));  xlabel('cos(\psi)');   ylabel('2fd/fr'); zlabel('P/dB');axis([min(x) max(x) min(y) max(y) min(min(10*log10(abs(Q1)/A))) max(max(10*log10(abs(Q1)/A)))]);
% figure();  mesh(x,y,10*log10(abs(Q3)/C));  xlabel('cos(\psi)');   ylabel('2fd/fr'); zlabel('P/dB'); axis([min(x) max(x) min(y) max(y) min(min(10*log10(abs(Q3)/C))) max(max(10*log10(abs(Q3)/C)))]);

figure();  mesh(z,y,10*log10(abs(Q3)/C));  xlabel('cos(\psi)');   ylabel('2fd/fr'); zlabel('P/dB'); axis([min(z) max(z) min(y) max(y) min(min(10*log10(abs(Q3)/C))) max(max(10*log10(abs(Q3)/C)))]);


