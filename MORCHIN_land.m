clear
clc
C = 299792458;
Fs = 1e9;
Lambda = C/Fs;

AA = [0.00126 0.004 0.0126 0.04];
BB = [pi/2 pi/2 pi/2 1.24];
beita00 = [0.14 0.2 0.4 0.5];
u = sqrt(Fs)/4.7;
yita = linspace(0,90,400)/180*pi;
for j=1:length(AA)
    
    A = AA(j);
    B = BB(j);
    beita0 = beita00(j);
    he = 9.3*beita0^2.2;
    yitac = asin(Lambda/(4*pi*he));
    k = 1.9;
    for i = 1:length(yita)
        if j==1
            if yita(i)<yitac
                seigemac0 = yita(i)/yitac;
            else
                seigemac0 = 1;
            end
        else
            seigemac0 = 1;
        end
        
        seigema0(i) = A*seigemac0*sin(yita(i))/Lambda+u*cot(beita0)^2*exp(-tan(B-yita(i))^2/tan(beita0)^2);
    end
    
    figure(1),
    plot(yita/pi*180,seigema0,'LineWidth',3)
    hold on
end
legend('沙漠','农田','丘陵','高山')
xlabel('入射角')
ylabel('地面散射系数（dB）')
title('地面杂波散射系数')
