clc
clear

%generating transmitter symbols
fc=1e2;%�ز�Ƶ��100HZ
symbolRate = 1e3;
nSamplesPerSymbol = fc;
symbolLength=1e4;
Fs=symbolRate*sqrt(nSamplesPerSymbol);%ʱ����Ƶ�ʲ���/��������
t_end = 1e3;
t=0:1/Fs:(t_end-1/Fs);
g=randi([0,1],1,symbolLength);

tz1=g(ceil(10*t+(1/Fs))).*cos(2*pi*fc*t);%�ز�Ƶ��Ϊ100HZ
tz2=~g(ceil(10*t+(1/Fs))).*cos(2*pi*2*fc*t);%�ز�Ƶ��Ϊ200HZ
txSamples=tz1+tz2;

% t1 = (1:1:1e3)./Fs;
% y1 = txSamples(1:1e3);
% subplot(2, 2, 1), plot(t1,abs(y1));
% xlabel('����ǰʱ��')
% subplot(2, 2, 2), plot(10*log10(abs(fftshift(fft(txSamples,2^15)))));
% ylabel('10\times log_{10}A');
% xlabel('����ǰƵ��')



SER = [];BER = []; 
for snr_dB = 0:2:30
%Add AWGN
%snr_dB = 9.8;
snr = 10^(snr_dB./10);
Psignal = mean(abs(txSamples.^2));%�źŹ��� 
Pnoise = Psignal./snr*nSamplesPerSymbol;
noise = sqrt(Pnoise).*(randn(1,length(txSamples))+1j*randn(1,length(txSamples)))./sqrt(2);%��������/sqrt2
%ע�⹦��
rxSamples = txSamples + noise;

% y2 = rxSamples(1:1e3);
% subplot(2, 2, 3), plot(t1,abs(y2));
% xlabel('�����ʱ��')
% subplot(2, 2, 4), plot(10*log10(abs(fftshift(fft(rxSamples,2^15)))));
% ylabel('10\times log_{10}A');
% xlabel('�����Ƶ��')

%coherent reception 
[b1,a1]=butter(2,[80,120]*2/Fs);%��ư�����˹��ͨ�˲�����100HZ�ź�ͨ����2�ף�ϵ��Ϊa1,b1
[b2,a2]=butter(2,[180,220]*2/Fs);%��ư�����˹��ͨ�˲�����200HZ�ź�ͨ����2�ף�ϵ��Ϊa2,b2

sg1=filter(b1,a1,rxSamples);%100HZ�ź�ͨ����BPF
sag1=filter(b2,a2,rxSamples);%200HZ�ź�ͨ����BPF

sg2=2*sg1.*cos(2*pi*fc*t);%100HZ�ź�ͨ�������
sag2=2*sag1.*cos(2*pi*2*fc*t);%200HZ�ź�ͨ�������

[b3,a3]=butter(2,10*2/Fs);%��ư�����˹��ͨ�˲�����100HZ�ź�ͨ�� 
[b4,a4]=butter(2,10*2/Fs);%��ư�����˹��ͨ�˲�����200HZ�ź�ͨ��
sg3= filter(b3,a3,sg2);%100HZ�ź�ͨ����LPF
sag3= filter(b4, a4, sag2);%200HZ�ź�ͨ����LPF


%�о�
LL=t_end/2;
for i=1:symbolLength
    if sg3((i-1)*t_end+LL)>=sag3((i-1)*t_end+LL)%ȡ�м�һ����Ϊ�о���Ƚ���·�ź����ֵ
        u(i)=1;
    else 
        u(i)=0;
    end
end


%������
[numbers,pe] = symerr(g,u);
BER_0 = pe;
BER = [BER BER_0];
end
figure;
semilogy(0:2:30,BER());
grid on;
legend('BFSK');
xlabel('SNR(dB)');
ylabel('������BER');
title('BER-SNR����');





