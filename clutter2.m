%*****�������Ӳ�ģ��******
clc
clear
tic
lambda=0.23;d=lambda/2;fr=2434.8;
M=8;N=6;K=8;
V=140;H=12000;V1=10;%V1Ϊ����
R0=50e3;    %Ŀ�����ھ��뻷б��
theta0=60*pi/180;%��������λ��
phi0=asin(H/R0);%������������
thetap=60*pi/180;  %����������������з���н�
theta_a=1*pi/180;    %ˮƽ��λ���ڲ������
pt=200e3;
Im=chebwin(M,20);In=chebwin(N,20);
CNR=60;
%***************
%****�״�ϵͳ����****
GG=pt*lambda^2/(4*pi)^3;%����
tau=0.1e-6;%�����ȣ�
R_cell=150;%����ֱ���
yangben=2*N*K;%
%************
Re=8490e3;c=3e8;
k=1.38e-23;T0=290;Bn=70e6;Fn=3.5;Fn=10^(3.5/10);
Pn=k*T0*Bn*Fn;
%******
theta=[0:pi/179:pi];
W=length(theta);
R2=zeros(N*K,N*K);
%***************
earth_leixing=2;%�����������Ƶر�����
clutter_leixing=2;%�����Ӳ���������
ops=0;
for l=1:yangben %�����Ų���ѵ������
    ops=ops+1;
    cr=zeros(N,K);
    r(l)=R0+(l-yangben/2)*R_cell;
    phi(l)=asin(H/r(l));
%     aaa=unifrnd(0,2*pi,1,1);
%      gaosi_noise=Pn*randn(N*K,1).*exp(j*2*pi*rand(N*K,1));
    for w=1:W   %��λ��ѭ��
        cos_psi(w)=cos(theta(w))*cos(phi(l));
        %*****�ռ䡢ʱ���Ƶ��
        ws(l,w)=2*pi*d/lambda*cos(theta(w))*cos(phi(l));
        wt(l,w)=4*pi*V/(lambda*fr)*cos(theta(w)+thetap)*cos(phi(l));
        fd(l,w)=2*V/lambda*cos(theta(w)+thetap)*cos(phi(l))/fr;
        %*****���߷���ͼ
        F(l,w)=sum(In'.*exp(j*2*pi*d/lambda*[0:N-1]*(cos(theta(w))-cos(theta0))))*sum(Im'.*exp(j*2*pi*d/lambda*[0:M-1]*(cos(theta(w))*cos(phi(l))-cos(theta0)*cos(phi0))));
        for n=1:N
            G(n)=sum(Im'.*exp(j*2*pi*d/lambda*[0:M-1]*(sin(phi(l))-sin(phi0))));
        end
        %****�Ӳ�ɢ�䵥Ԫ���**
        Ac=r(l)*theta_a*c*tau/2;    %�󸩽�������Ӳ�ɢ�䵥Ԫ�����theta_a����λ��ƽ���ڲ�����ȡ�
        %***�Ӳ�ɢ��ϵ��sigma0*****
        if earth_leixing==1 %ũ��
            A=0.004;B=pi/2;beta0=0.2;sigma_c0=1;
        elseif earth_leixing==2 %����
            A=0.0126;B=pi/2;beta0=0.4;sigma_c0=1;
        elseif earth_leixing==3 %��ɽ
            A=0.04;B=1.24;beta0=0.5;sigma_c0=1;
        end
        u=sqrt(c/lambda*1e-9)/4.7;
        theta_g=asin(H/r(l)+H^2/(2*Re*r(l))-r(l)/(2*Re));
        sigma0(l)=A*sigma_c0*sin(theta_g)/lambda+u*cot(beta0)^2*exp(-tan(B-theta_g)^2/tan(beta0)^2);
        %***�Ӳ��ֱ浥Ԫ�״�ɢ������
        sigma_c=sigma0(l)*Ac;
        %***�Ӳ����ȿ���
        if clutter_leixing==1   %�����ֲ�
            al=raylrnd(1.5,1,K);
        elseif clutter_leixing==2   %������̬�ֲ�
            al=lognrnd(0.8,0.2,1,K);
        elseif clutter_leixing==3
            al=wblrnd(1.5,2.5,1,K);
        end
        %****�Ӳ�������Ƶ�׿���**��˹��
        sigma_v1=0.0066*V1;sigma_v2=V/(2*sqrt(2*log(2)))*theta_a;
        sigma_v=sqrt(sigma_v1^2+sigma_v2^2);
        sigma_f=2*sigma_v/lambda;
        delta_f=0.6*(2*sqrt(2*log(2))*sigma_f);
        sk=delta_f*1/sqrt(2*pi*sigma_f^2)*exp(-(([1:5]-3)*delta_f).^2/(2*sigma_f^2));
        kesai=exp(1i*2*pi*rand(1,5));
        xk=kesai.*sqrt(sk);
        for k=1:K
            X(k)=sum(xk.*exp(1i*2*pi*([1:5]-3)*delta_f*k*1/fr));
        end       
        %***ϵͳ������**
        gaosi_noise=sqrt(10)*randn(1,K);%��ֵΪ0������Ϊ1�ĸ�˹�ֲ�������
%         Z=sqrt(-2*log(randn(1,K)));
%         gaosi_noise=Pn*Z.*exp(j*unifrnd(0,2*pi,1,K));%������˹����
       aaa=unifrnd(0,2*pi,1,1);
        %***�Ӳ��ز��ź�����
        for n=1:N       %
%             cr(n,:)=cr(n,:)+sqrt(GG)*sqrt(sigma_c)*F(l,w)*G(n).*al.*exp(j*aaa).*X.*exp(j*(n-1)*ws(l,w)+j*[0:K-1]*wt(l,w))/r(l)^2+gaosi_noise;
            cr(n,:)=cr(n,:)+(sqrt(sigma_c)*F(l,w)*G(n)*al.*exp(1i*unifrnd(-pi,pi,1,K)).*X.*exp(1i*(n-1)*ws(l,w)+1i*((1:K)-1)*wt(l,w)))/r(l)^2+gaosi_noise;
        end %��Ԫѭ��
    end %��λ��ѭ��
    crr=reshape(cr,N*K,1);%+gaosi_noise;
    R1=crr*crr';R2=R2+R1;        
end %������ѭ��
figure(1)   %�Ƕȡ��������ռ�����λ-�����ռ���
plot(fd(1,:),cos_psi(:),'*-')
hold on
plot(fd(40,:),cos_psi(:),'*-')
hold on
grid on
hold off
R2=R2/yangben;
R=R2/(max(max(abs(R2))));
R=10^(CNR/10)*R+eye(N*K);
R_inv=inv(R);
[u,s,v]=svd(R);
figure(2)   %�Ӳ�����ֵ�����ɶȣ�
plot(10*log10(diag(s)),'*-');   
grid on
figure(3)   %���߷���ͼ
mesh(10*log10(abs(F)))
plot(theta/pi*180,10*log10(abs(F(48,:)))) 

figure(4)  %�Ӳ�������
mesh(abs(R))


wss=[-1:2/80:1];
wtt=[-1:2/80:1];
for ss=1:length(wss)
    Ss=conj(exp(1i*[0:N-1]*ws(ss))');
    for tt=1:length(wtt)
        St=conj(exp(1i*[0:K-1]*wtt(tt))');
        SS=kron(Ss,St);
        mu=inv(SS'*R_inv*SS);
        P(ss,tt)=1/abs((SS'*R_inv*SS));
%         Wopt(ss,tt)=mu*R_inv*SS;
        IF(ss,tt)=SS'*R_inv*SS*trace(R)/(SS'*SS);
    end
end

P=P/max(max(abs(P)));
figure(5)   %�Ӳ�������
mesh(wss,wtt,10*log10(P))
figure(6)   %��������
mesh(wss,wtt,10*log10(IF))

toc
