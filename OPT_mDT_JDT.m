%��ʱ��ά3dt����
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
anoise=max(max(noise));                          %��������
psi0=pi/2; % Ŀ���׶������λ��Ϊ0

% ȫάSTAP
fd=-1:1/50:1;
inv_Rx=pinv(Rx);
for i=1:length(fd)
    Ss=exp(j*pi*(0:N-1)'*cos(psi0));            %Ŀ�귽��ȷ��ʱ��Ss�̶�����Ŀ��dopplerƵ��δ֪����ÿһ��fd��һ������Ȩʸ��wopt
    St=exp(j*pi*(0:K-1)'*fd(i));
    S=kron(St,Ss); % ֱ���øö�׶�����Һ͹�һ��������Ƶ�����ɵ���ʸ��
    wopt=inv_Rx*S/(S'*inv_Rx*S);
    IF(i)=abs(wopt'*S)^2*(10^(CNR/10)+1)*anoise/(wopt'*Rx*wopt);
end
figure
plot(fd,10*log10(abs(IF)))
xlabel('2f_d/f_r');ylabel('IF/dB');
grid on

% 3DT��άSTAP
index=0;
for fd=-1:1/50:1
    index=index+1;
    Qt=[exp(j*pi*(0:K-1).'*(fd - 1/K)),exp(j*pi*(0:K-1).'*fd),exp(j*pi*(0:K-1).'*(fd + 1/K))];%ʱ��ά����
    Qs=eye(N); % ����ά����
    Q=kron(Qt,Qs);
    Ry=Q'*Rx*Q; % ���ݽ�ά֮���Э����������ݽ�ά����������48*48
    inv_Ry=pinv(Ry);

    Ss=exp(j*pi*(0:N-1).'*cos(psi0));
    St=exp(j*pi*(0:K-1).'*fd);
    S = kron(St,Ss);

    Sy = Q'*S;
    W_3dt=inv_Ry*Sy/(Sy'*inv_Ry*Sy); % ���ݾ�����Ҫ��ά��ʹ�õĵ��������Ҫ����������Ҫ����ͨ���Ĺ�һ������
    IF_3dt(index)=abs(W_3dt'*Sy).^2*(10^(CNR/10)+1)*anoise/(W_3dt'*Ry*W_3dt);
end
hold on
plot(-1:1/50:1,10*log10(IF_3dt),'.-');

% JDL��άSTAP
index = 0;
for fd=-1:1/50:1
    index = index + 1;
    % �������ݽ�ά��ɸѡĿ�긽���ľŹ�������
    Qt = [exp(j*pi*(0:K-1).'*(fd-1/K)),exp(j*pi*(0:K-1).'*fd),exp(j*pi*(0:K-1).'*(fd+1/K))];%ʱ��ά����
    Qs = [exp(j*pi*(0:N-1).'*(cos(psi0)-1/N)),exp(j*pi*(0:N-1).'*cos(psi0)),exp(j*pi*(0:N-1).'*(cos(psi0)+1/N))];%����ά����

    Q = kron(Qt,Qs);
    Rz = Q'*Rx*Q; % ���ݽ�ά���Э������󣬱���������9*9
    inv_Rz = pinv(Rz);
    Ss=exp(j*pi*(0:N-1).'*cos(psi0));
    St=exp(j*pi*(0:K-1).'*fd);
    S = kron(St,Ss);

    % ����µĿռ䵼��ʸ�����ռ䵼��ʸ��ά��Ϊ9
    Sz = Q'*S;
    w_jdl = inv_Rz*Sz/(Sz'*inv_Rz*Sz); 
    IF_jdl(index)=abs(w_jdl'*Sz).^2*(10^(CNR/10)+1)*anoise/(w_jdl'*Rz*w_jdl);
end
hold on
plot(-1:1/50:1,10*log10(IF_jdl),'.-');
legend('����','3DT','JDL')