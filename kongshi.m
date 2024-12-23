%                �Ӳ���ģ     ������
%  ======================================================================
clc;    close all;  clear all;
j=sqrt(-1);
alpha =0*pi/180;      %����� �ı��Ӳ���״ 0*pi/180   45*pi/180    90*pi/180
M=8; %��Ԫ����
N=10;%��Ԫ����
K=16;                   % ����������
N2=N;
N1=M;
lamda=0.32;
% lambda=0.12;

Va=120;                 % �ɻ��ٶ�120m/s
H=1000;                 % �ɻ��߶�
Vc=0;                   % Ŀ�귽λ���ٶ�
Vr=100;                 % Ŀ��������ٶ�
prf=6000;               % PRF�������ظ�Ƶ��
c=3*10^8;
Ru=c/(2*prf);           % ���ģ������,100km
Rmax=400000;            % �״�Զ��б�� 100km
Rmin=1000;              % ���ؾ���  �ɻ����·���ʼ
Re = 8490000;
theta0=90/180*pi;       % ����ָ��ķ�λ��
phi0=30/180*pi;         % ����ָ��ĸ�����
theta=[0:180]/180*pi;   % ��λ�Ƿ�λ��0~180�㣩
distance=15;            % ���뻷���
dis = distance;
L=200;                  % ���뻷������L>2MN
dd = lamda/2;


R_top=0.203;%0.203;
h=0.182;



S=1; %�����Ǿ���ģ��
for i=1:L
     Rl=Rmin+(S-1)*Ru+(i-1)*dis;
     phi(i)=asin(H/Rl+(Rl^2-H^2)/(2*Rl*(Re+H)));   %ÿ�����뻷��Ӧ�ĸ�����  ��λ ��
end

arry_col=M;                                    %��ѡ����Ԫ�ϳ� �������� (ż���ɱ�N2sel����)
arry_row=N;                                     %��ѡ����Ԫ�ϳ� �������� ż�� �ɱ�N1sel����
T8=kron(eye(arry_row),kron(ones(N/arry_row,1),kron(eye(arry_col),ones(M/arry_col,1)))); %���ɾ۽�����T


% 
% px = linspace(0,lamda/2*(N-1),N);      %��Ч���������
% pz = linspace(0,lamda/2*(M-1),M);      %��Ч����������


px = linspace(-R_top,R_top,N);      %��Ч���������
pz = linspace(0,h,M);             %��Ч����������

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
title('��Ч��Ԫλ��')

d=2*R_top/N;

SCR=30;                %���ӱ�
CNR=60;                 %�����
%%%%%%-----------���߷���ͼ-------------
dB=20;
IM=chebwin(M,dB);   %�������е���ǿ�ȣ��Ӵ�                          ??
IN=chebwin(N,dB);   %�������е���ǿ�ȣ��Ӵ�

F_arr=zeros(N,1);                     
Frow = zeros(1,length(theta));

theta3db=1*pi/180;                          %�빦�ʲ������
r=zeros(1,N*K);

%%%%%%%%%%------------�Ӳ�����------------%%%Ŀ�����ڵķ�λ��������ָ���λ��        ���ٶ�Ŀ��λ���м���뻷�� 
clutter1=zeros(N*K,1);          %���ڴ洢ÿһ�ο��ģ����뻷������ʸ�����Ӳ����źţ���   ��ʱ���
clutter2=zeros(N*K,1);          %���ڴ洢ÿһ�ο��ģ����뻷������ʸ�����Ӳ���
signal=zeros(N*K,1);            %Ŀ���źž���
clutter1_matrix=zeros(N*K,L);   %���ڴ洢���о��뻷���Ӳ����ź� ����
clutter2_matrix=zeros(N*K,L);   %���ڴ洢��ȥ��Ŀ������о��뻷�Ӳ����� ��������Ŀ�����ھ��뻷���Ӳ���
Rcs_clutter_signal=zeros(N*K,N*K,L);    %��ά���� ���ڴ洢ÿһ�����뻷һ�ο��ĵ���ؾ���δƽ�����������źţ�
Rc_clutter=zeros(N*K,N*K,L);            %��ά���� ���ڴ洢ÿһ�����뻷һ�ο��ĵ���ؾ���δƽ�������������źţ�
Rcs=zeros(N*K);                         %������Ŀ����뻷���Ƶ�Э�������
Rc=zeros(N*K);                          %������Ŀ����뻷����������Ŀ���źŹ��Ƶ�Э������� ���Ӳ������� Э�������
for i=1:L
    i;
    Vr_temp = Vr;
    
    %% ����
    ss=5; %����ȼ�
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
    
    
    
    B_coe=seigema0*exp(1j*2*pi*randn(1,length(theta)));        %��i�����뻷�ϵ��Ӳ�����ϵ��,���ȷ�����̬�ֲ�����λ������̬�ֲ�
    B_coe_target=exp(1j*2*pi*randn(1,1));
    
    uu = cos(phi(i))*sin(theta);
    vv = cos(phi(i))*cos(theta);
    hh = sin(phi(i));
     
    
    % �����кϳɽ��շ���ͼ����
    w4 = repmat(exp(1j*2*pi*dd/lamda*(sin(phi(i))-sin(phi0)).*p(2,:)'),1,length(theta));
       
    F_arr = T8'*(w4);
    
    
    % ����  ���㷢�䷽��ͼ
    for m=1:length(theta)
        ww = exp(1j*2*pi*dd/lamda*(sin(phi(i))-sin(phi0)).*p(2,:)').*exp(1j*2*pi*dd/lamda*(cos(theta(m))*cos(phi(i))-cos(theta0)*cos(phi0)).*p(1,:)');
        Frow(1,m) = sum(ww);
    end
    
    for k=1:K        %ÿ�����뻷�Ŀ����źţ����������
        for n=1:N    %��Ԫ�ź�
            a=4*pi*Va/lamda/prf;
            Wt=4*pi*Va/lamda/prf*cos(theta +alpha)*cos(phi(i));        %��k�ο���ʱ��ʱ���Ƶ��   i���뻷�Ķ����ս�Ƶ��        
            Ws=4*pi*d/lamda*cos(theta)*cos(phi(i));                    %��n�����ߵĿ����Ƶ��
            s=exp(j*(k-1)*Wt+j*(n-1)*Ws)/(Rmin+(i-1)*dis)^2;                        %�ź�ǿ��˥�� 1/R^2
            clutter1((k-1)*N+n,1)=sum(F_arr(n).*Frow.*B_coe.*s);       %ffΪ����?��f^2*B Ϊ�ز��źţ󣨣���
            clutter2((k-1)*N+n,1)=sum(F_arr(n).*Frow.*B_coe.*s);
            r(1,(k-1)*N+n)=sum(F_arr(n).*(Frow).*exp(j*((n-1)*Ws+(k-1)*Wt))/(Rmin+(i-1)*dis)^2); 
            %    -----     �����ź�    ------
            if i==L/2   
                Vr_temp = Vr;
                for mm = 1: 2
                    Vr_temp = Vr_temp+100;
                    Wt=4*pi*(Va-Vc)/lamda/prf*cos(theta0+alpha)*cos(phi(i))+4*pi*Vr_temp/lamda/prf*  sin(theta0+alpha)*cos(phi(i));      %Ŀ���ʱ���Ƶ��   �����ٶ�VcVr
                    Ws=4*pi*d/lamda*cos(theta0)*cos(phi(i));                                                               %Ŀ��Ŀ����Ƶ��
                    signal((k-1)*N+n,1)=sum(sqrt(10^(SCR/10))*F_arr(n).*Frow.*B_coe_target.*exp(j*(k-1)*Wt+j*(n-1)*Ws))/(Rmin+(i-1)*dis)^2;%(k-1)*N+n ��k-1������ �ĵ�n���ռ䵼���ź�
                    clutter1((k-1)*N+n,1)=clutter1((k-1)*N+n,1)+signal((k-1)*N+n,1);
                end
            end
            
        end
    end
    clutter1_matrix(:,i)=clutter1;
    clutter2_matrix(:,i)=clutter2;
    Rcs_clutter_signal(:,:,i)=clutter1_matrix(:,i)*clutter1_matrix(:,i)';    %����ؾ���?�������롡���������롡��ʱ�䣩������ؾ���
    Rc_clutter(:,:,i)=clutter2_matrix(:,i)*clutter2_matrix(:,i)';
end


%%%%%%%%�����ھ��뻷L/2���Ӳ�Э������󣨶�Ŀ��λ��L/2���뻷����%%%%%%%%%%%%%%%%
for i=1:L         %��2NK����������Э�������
    Rcs=Rcs+Rcs_clutter_signal(:,:,i);
end
Rcs=Rcs/L;             %RcsЭ��������а���Ŀ���ź�

noise=max(max(Rcs))/10^(CNR/10)*eye(size(Rcs));         %������ֻ�������ʱ��Ϊ��

Rcs=Rcs+noise;

Rcs_v=pinv(Rcs);

%%%% capon����
%%%%--�Ӳ������ף���������άƵ��ͼ���Ӳ����źŵĹ����׺�Ƶ��ͼ%��ʵ���Ƕ�ά������ȡһ��׶�ǣ�Ȼ���ڴ�׶�ǹ̶�ʱ�������еĶ�����Ƶ��%%%%%%%%%%
psi=acos(linspace(-1,1,361));
fd=linspace(-1/2,1/2,361);
for s=1:length(psi)
    for t=1:length(fd)
        LL=exp(j*2*pi*(0:K-1)'*fd(t));
        PP=exp(j*pi*(0:N-1)'*cos(psi(s)));
        S=kron(LL,PP);
%         Q1(s,t)=1/(S'*Rc_v*S);   %��
        temp = Rcs_v*S;
        Q3(s,t)=1/(S'*Rcs_v*S);  %��
%         Q3(s,t)=S'*Rcs*S;  %��
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


