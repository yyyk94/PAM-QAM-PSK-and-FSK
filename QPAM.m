clc
clear

%generating transmitter symbols
bitLength = 2^16;
M = 4;  %4PAM
symbolLength = bitLength./log2(M);

txBits = randi([0,1],1,bitLength);

modTable = [-3 -1 1 3];

txBitsIndex = reshape(txBits,2 ,symbolLength);  %分组
txSymbols = modTable([1 2]*txBitsIndex+1);%转为十进制

P_Symbols = mean(abs(txSymbols).^2);

%pulse shaping
nSamplesPerSymbol = 2 ;
span = 60; %单边长度
roll = 0.1; %滚降系数
%画图
color = ['r','g','b','y','m'];
BER_f = [];

symbolRate = 100e6;
Fs = symbolRate*nSamplesPerSymbol;
t_end = 2^10;
t = (1:t_end)/Fs;
rrcFilter = rcosdesign(roll, span ,nSamplesPerSymbol);
txSamples  = upfirdn(txSymbols, rrcFilter, nSamplesPerSymbol, 1);
txSamples = txSamples./sqrt( mean(abs(txSamples).^2));
% txSamples = upsample(txSymbols,nSamplesPerSymbol);
% txSamples = conv(txSamples,rrcFilter,'same');


% y1 = txSamples(1:t_end);
% subplot(2, 2, 1), plot(t,abs(y1));
% xlabel('加噪前时域')
% subplot(2, 2, 2), plot(10*log10(abs(fftshift(fft(txSamples)))));
% ylabel('10\times log_{10}A');
% xlabel('加噪前频域')





SER = [];BER = []; 
for snr_dB = 5:2:30

%Add AWGN
% snr_dB = 15;
snr = 10^(snr_dB./10);
Psignal = mean(abs(txSamples.^2));%信号功率 
Pnoise = Psignal./snr*nSamplesPerSymbol;
noise = sqrt(Pnoise).*(randn(1,length(txSamples))+1j*randn(1,length(txSamples)))./sqrt(2);%复数所以/sqrt2
%注意功率
rxSamples = txSamples + noise;


% y2 = txSamples(1:t_end);
% subplot(2, 2, 3), plot(t,abs(y2));
% xlabel('加噪后时域')
% subplot(2, 2, 4), plot(10*log10(abs(fftshift(fft(rxSamples)))));
% ylabel('10\times log_{10}A');
% xlabel('加噪后频域')
% scatterplot(rxSamples);

%matched filtering
rxSymbols = upfirdn(rxSamples,rrcFilter,1,nSamplesPerSymbol);
rxSymbols = rxSymbols(span+1:end-span);
rxSymbols = rxSymbols./sqrt(mean(abs(rxSymbols).^2));
rxSymbols = rxSymbols.*sqrt(P_Symbols);
% scatterplot(rxSymbols);
        

%demapping
        
for i = 1:M
    dist(i,:) = abs(rxSymbols-modTable(i)).^2;      
end

[~,ind]=min(dist);
decSymbols = modTable(ind);

SER = sum(abs(decSymbols-txSymbols)>0.01)./symbolLength;

rxBits = zeros(1,bitLength);     
rxBits(1:2:end) = bitget(ind-1,1);     
rxBits(2:2:end) = bitget(ind-1,2);     

%BER =sum(abs(rxBits-txBits)>0.01)./bitLength;
BER = [BER sum(abs(rxBits-txBits)>0.01)./bitLength];
end
figure;
semilogy(5:2:30,BER(),'m');
grid on;
legend('QPAM');
xlabel('信噪比SNR/dB');
ylabel('误码率BER');
title('BER-SNR曲线');
