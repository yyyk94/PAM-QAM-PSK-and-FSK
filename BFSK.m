clc
clear

%generating transmitter symbols
fc=1e2;%载波频率100HZ
symbolRate = 1e3;
nSamplesPerSymbol = fc;
symbolLength=1e4;
Fs=symbolRate*sqrt(nSamplesPerSymbol);%时间轴频率步进/采样速率
t_end = 1e3;
t=0:1/Fs:(t_end-1/Fs);
g=randi([0,1],1,symbolLength);

tz1=g(ceil(10*t+(1/Fs))).*cos(2*pi*fc*t);%载波频率为100HZ
tz2=~g(ceil(10*t+(1/Fs))).*cos(2*pi*2*fc*t);%载波频率为200HZ
txSamples=tz1+tz2;

% t1 = (1:1:1e3)./Fs;
% y1 = txSamples(1:1e3);
% subplot(2, 2, 1), plot(t1,abs(y1));
% xlabel('加噪前时域')
% subplot(2, 2, 2), plot(10*log10(abs(fftshift(fft(txSamples,2^15)))));
% ylabel('10\times log_{10}A');
% xlabel('加噪前频域')



SER = [];BER = []; 
for snr_dB = 0:2:30
%Add AWGN
%snr_dB = 9.8;
snr = 10^(snr_dB./10);
Psignal = mean(abs(txSamples.^2));%信号功率 
Pnoise = Psignal./snr*nSamplesPerSymbol;
noise = sqrt(Pnoise).*(randn(1,length(txSamples))+1j*randn(1,length(txSamples)))./sqrt(2);%复数所以/sqrt2
%注意功率
rxSamples = txSamples + noise;

% y2 = rxSamples(1:1e3);
% subplot(2, 2, 3), plot(t1,abs(y2));
% xlabel('加噪后时域')
% subplot(2, 2, 4), plot(10*log10(abs(fftshift(fft(rxSamples,2^15)))));
% ylabel('10\times log_{10}A');
% xlabel('加噪后频域')

%coherent reception 
[b1,a1]=butter(2,[80,120]*2/Fs);%设计巴特沃斯带通滤波器，100HZ信号通过，2阶，系数为a1,b1
[b2,a2]=butter(2,[180,220]*2/Fs);%设计巴特沃斯带通滤波器，200HZ信号通过，2阶，系数为a2,b2

sg1=filter(b1,a1,rxSamples);%100HZ信号通过该BPF
sag1=filter(b2,a2,rxSamples);%200HZ信号通过该BPF

sg2=2*sg1.*cos(2*pi*fc*t);%100HZ信号通过相乘器
sag2=2*sag1.*cos(2*pi*2*fc*t);%200HZ信号通过相乘器

[b3,a3]=butter(2,10*2/Fs);%设计巴特沃斯低通滤波器，100HZ信号通过 
[b4,a4]=butter(2,10*2/Fs);%设计巴特沃斯低通滤波器，200HZ信号通过
sg3= filter(b3,a3,sg2);%100HZ信号通过该LPF
sag3= filter(b4, a4, sag2);%200HZ信号通过该LPF


%判决
LL=t_end/2;
for i=1:symbolLength
    if sg3((i-1)*t_end+LL)>=sag3((i-1)*t_end+LL)%取中间一点作为判决点比较两路信号输出值
        u(i)=1;
    else 
        u(i)=0;
    end
end


%误码率
[numbers,pe] = symerr(g,u);
BER_0 = pe;
BER = [BER BER_0];
end
figure;
semilogy(0:2:30,BER());
grid on;
legend('BFSK');
xlabel('SNR(dB)');
ylabel('误码率BER');
title('BER-SNR曲线');





