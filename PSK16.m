clc
clear

%generating transmitter symbols
bitLength = 2^16;
M = 16;  %16PSK
symbolLength = bitLength./log2(M);
txBits = randi([0,1],1,bitLength);

%modTable = [exp(i*pi/4) exp(-i*pi/4) exp(i*3/4*pi) exp(-i*3/4*pi)];
modTable = [exp(i*pi/16) exp(-i*pi/16) exp(i*3/16*pi) exp(-i*3/16*pi) exp(i*5/16*pi) exp(-i*5/16*pi) exp(i*7/16*pi) exp(-i*7/16*pi)...... 
    exp(i*9/16*pi) exp(-i*9/16*pi) exp(i*11/16*pi) exp(-i*11/16*pi)  exp(i*13/16*pi) exp(-i*13/16*pi)  exp(i*15/16*pi) exp(-i*15/16*pi)];

txBitsIndex = reshape(txBits,4,symbolLength);  %分组

txSymbols = modTable([1 2 4 8]*txBitsIndex +1);

%pulse shaping
nSamplesPerSymbol = 4 ;
span = 60; %单边长度
roll = 0.1; %滚降系数
%画图
color = ['r','g','b','y','m'];
BER_f = [];

t_end = 2^10;
symbolRate = 100e5;
Fs = symbolRate*nSamplesPerSymbol;
t = (1:t_end)/Fs;
rrcFilter = rcosdesign(roll, span ,nSamplesPerSymbol);
% txSamples = upsample(txSymbols,nSamplesPerSymbol);
% txSamples = conv(txSamples,rrcFilter,'same');
txSamples  = upfirdn(txSymbols, rrcFilter, nSamplesPerSymbol, 1);
txSamples = txSamples./ sqrt(mean(abs(txSamples).^2));

% y1 = txSamples(1:t_end);
% subplot(2, 2, 1), plot(t,abs(y1));
% xlabel('加噪前时域')
% subplot(2, 2, 2), plot(10*log10(abs(fftshift(fft(txSamples)))));
% ylabel('10\times log_{10}A');
% xlabel('加噪前频域')

SER = [];BER = []; 
for snr_dB = 5:2:30
%Add AWGN
%snr_dB = 20;
snr = 10^(snr_dB./10);
Psignal = mean(abs(txSymbols.^2));%信号功率 
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
% rxSymbols = downsample(rxSamples,nSamplesPerSymbol);
% rxSymbols = conv(rxSymbols,rrcFilter,'same');
rxSymbols = rxSymbols./sqrt(mean(abs(rxSymbols).^2));
%scatterplot(rxSymbols);
        

%demapping
        
for i = 1:M
    dist(i,:) = abs(rxSymbols-modTable(i)).^2;      
end

[~,ind]=min(dist);
decSymbols = modTable(ind);

SER = sum(abs(decSymbols-txSymbols)>0.01)./symbolLength;

rxBits = zeros(1,bitLength);     
rxBits(1:4:end)= mod((ind-1),2); %十进制低位
rxBits(2:4:end)=mod(((ind-1)-rxBits(1:4:end))./2,2); 
rxBits(3:4:end)=mod((((ind-1)-rxBits(1:4:end))./2-rxBits(2:4:end))./2,2);
rxBits(4:4:end)=((((ind-1)-rxBits(1:4:end))./2-rxBits(2:4:end))./2-rxBits(3:4:end))./2;


%BER =sum(abs(rxBits-txBits)>0.01)./bitLength;
BER = [BER sum(abs(rxBits-txBits)>0.01)./bitLength];
end
figure;
semilogy(5:2:30,BER(),'g');
grid on;
legend('16PSK');
xlabel('信噪比SNR/dB');
ylabel('误码率BER');
title('BER-SNR曲线');