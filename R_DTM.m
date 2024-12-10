clc;
% close all;
dbstop if error;
fileName = 'backward1T4R.bin';
%% 基本参数设置
numOfADCSamples = 128;
numOfChirpsInFrame = 128;
adcSampleByte = 2;
IQ = 2;
numOfADCBits = 16;
TX = 2;
RX = 4;
adcBinInfo = dir(fileName);
adcBinSize = adcBinInfo.bytes;
numOfFrames = adcBinSize/numOfADCSamples/numOfChirpsInFrame/TX/RX/adcSampleByte/IQ;
numOfAllChirps = numOfChirpsInFrame*numOfFrames;
fs = 4000e3;%ksps
c = 3e8;
ts = numOfADCSamples/fs;
freqSlope = 99.987e12;
b = ts*freqSlope;
timePerFrame = 50e-3;%50e-3;
% frequencyPerChirp = 1/timePerChirp;
startFrequency = 77e9;
deltaR = c/(2*b);
deltaV = (c/startFrequency)/(2*numOfFrames*timePerFrame);
%% bin文件转化为adc数据
radarData = readDCA1000(fileName,numOfADCSamples,TX*RX);%TI官方bin文件处理
radarDataFinal = reshape(radarData(1,:),numOfADCSamples,numOfAllChirps);
samplesInFastTime = numOfADCSamples;
samplesInSlowTime = numOfFrames;
time = (0:samplesInSlowTime-1)*timePerFrame;
%% RDM
clc;

hamWinR = hamming(samplesInFastTime);
hamWinD = hamming(numOfChirpsInFrame);

radarRTM = zeros(samplesInFastTime,samplesInSlowTime);
radarDTM = zeros(samplesInFastTime,samplesInSlowTime);

temp = zeros(samplesInFastTime,numOfChirpsInFrame);
fft_1d = zeros(samplesInFastTime,numOfChirpsInFrame);
fft_2d = zeros(samplesInFastTime,numOfChirpsInFrame);
for i = 1:numOfFrames
    fft_1d = fft(radarDataFinal(:,(i-1)*numOfChirpsInFrame+1:i*numOfChirpsInFrame));%帧内速度维fft
    radarRTM(:,i) = fft_1d(:,1) - mean(fft_1d,2);%静物去除
    fft_1d = fft_1d - mean(fft_1d,2);%零速去除
    for k = 1:samplesInFastTime
        fft_2d(k,:) = fftshift(hamWinD'.*fft(fft_1d(k,:)));%二维fft doppler
    end
    radarDTM(:,i) = mean(fft_2d,1);%average
end

figure;
imagesc(time,deltaR*(0:samplesInFastTime-1),abs(radarRTM));
xlabel('Time(s)');ylabel('Range(m)');title('RTM');
colorbar;

figure;
imagesc(time,deltaV*(-samplesInFastTime/2:samplesInFastTime/2-1),abs(radarDTM));
xlabel('Time(s)');ylabel('Velocity(m/s)');title('DTM');
colorbar;
