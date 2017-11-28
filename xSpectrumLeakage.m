% xSpectrumLeakage.m
% 徐文宇，20171108
% 分析频谱泄漏
% ver:---
close all;clear;

% plot画线类型
plotLineType='';        % '' 实线
plotDottedLineType=':'; % ':'虚线

%% @_@{整数周期信号}*****************************************************
% 采样周期(时间间隔)
timeSample = 0.000625;
% 采样频率，采样周期的倒数
frequencySample = 1/timeSample;
% 原始信号的频率
frequencySignal = 80;
% 每个周期的采样点数
samplingPointsNumEachPeriod = frequencySample / frequencySignal;

% 采样周期数目
samplingPeriodNum = 8;

% 采样起止时间
timeStart = 0;
timeEnd = timeStart+samplingPeriodNum*samplingPointsNumEachPeriod *timeSample-timeSample;
% 采样时间序列
timeSequence = timeStart:timeSample:timeEnd;
% 采样时长
duration = timeEnd-timeStart+timeSample;
% 采样点数目
samplingPointsNum = duration/timeSample+1;
% 采样频率序列
frequencySequency = 0:frequencySample/(samplingPointsNum-1):frequencySample-frequencySample/(samplingPointsNum-1);

% 对原始信号进行采样
dataSample = sin(2*pi*frequencySignal*timeSequence);%+2*sin(2*pi*40*t);
% 对采样信号进行Fourier变换
fftDataSample = fft(dataSample);

% 绘制原始信号
figure('name','Integate Periods of Signal','NumberTitle','off');
plot(timeSequence,dataSample,'*-');
title('Integate Periods of Signal');
grid on;
xlim([timeStart,timeStart+ceil(samplingPeriodNum)*samplingPointsNumEachPeriod *timeSample]);
ylim([min(dataSample)-0.1,max(dataSample)+0.1]);
% xTick & xTickLabel
xTick = timeStart:timeSample*samplingPointsNumEachPeriod:timeStart+ceil(samplingPeriodNum)*samplingPointsNumEachPeriod *timeSample;
xTickLabel = cell(1,ceil(samplingPeriodNum));
for k = 1:length(xTick)
    if k==1
        xTickLabel{k} = sprintf('%.4f',xTick(k));
    else
        xTickLabel{k} = sprintf('[%d]%.4f',k-1,xTick(k));
    end
end
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
xlabel('Time (/$s$)','Interpreter','latex');
ylabel('Amplitude');

% 绘制频谱
figure('name','Fourier Spectrum (Integate Periods of Signal)','NumberTitle','off');
stem(frequencySequency-frequencySample/2,abs(fftshift(fftDataSample)));
title('Fourier Spectrum (Integate Periods of Signal)');
grid on;
xlim([min(frequencySequency-frequencySample/2),max(frequencySequency-frequencySample/2)]);
% xTick & xTickLabel
xTick2 = [-ceil(frequencySample/2):frequencySignal*2:-frequencySignal*2,-frequencySignal,...
    0,frequencySignal,frequencySignal*2:frequencySignal*2:ceil(frequencySample/2)];
xTickLabel2 = cell(1,length(xTick2));
for k = 1:length(xTick2)
  xTickLabel2{k} = num2str(xTick2(k));
end
set(gca, 'XTick', xTick2);set(gca, 'XTickLabel',xTickLabel2);
xlabel('Frequency (/$Hz$)','Interpreter','latex');
ylabel('Magnitude');

% 显示整数周期信号及其Hilbert变换
figure('name','Integate Periods of Signal and its Hilbert Transform','NumberTitle','off');
% 整数周期信号
plot(timeSequence,dataSample,                       plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 期望Hilbert变换结果
expectedHilbert=-cos(2*pi*frequencySignal*timeSequence);
plot(timeSequence,expectedHilbert,plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 实际Hilbert变换结果
plot(timeSequence,imag(hilbert(dataSample)),        plotLineType,'Color','m','LineWidth',0.5,'MarkerSize',2);
title('Integate Periods of Signal and its Hilbert Transform');
legend('Signal','Expected HT','Actual HT','Location','southoutside','Orientation','horizontal');
xlim([timeStart,timeStart+ceil(samplingPeriodNum)*samplingPointsNumEachPeriod *timeSample]);
ylim([min(dataSample)-0.1,max(dataSample)+0.1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
grid on;


%% @_@{非整数周期信号}**************************************************
% 采样周期(时间间隔)
timeSample = 0.000625;
% 采样频率，采样周期的倒数
frequencySample = 1/timeSample;
% 原始信号的频率
frequencySignal = 80;
% 每个周期的采样点数
samplingPointsNumEachPeriod = frequencySample / frequencySignal;

% 采样周期数目
samplingPeriodNum = 8.5;

% 采样起止时间
timeStart = 0;
timeEnd = timeStart+samplingPeriodNum*samplingPointsNumEachPeriod *timeSample-timeSample;
% 采样时间序列
timeSequence = timeStart:timeSample:timeEnd;
% 采样时长
duration = timeEnd-timeStart+timeSample;
% 采样点数目
samplingPointsNum = duration/timeSample+1;
% 采样频率序列
frequencySequency = 0:frequencySample/(samplingPointsNum-1):frequencySample-frequencySample/(samplingPointsNum-1);

% 对原始信号进行采样
dataSample = sin(2*pi*frequencySignal*timeSequence);%+2*sin(2*pi*40*t);
% 对采样信号进行Fourier变换
fftDataSample = fft(dataSample);

% 绘制原始信号
figure('name','Non-Integate Periods of Signal','NumberTitle','off');
plot(timeSequence,dataSample,'*-');
title('Non-Integate Periods of Signal');
grid on;
xlim([timeStart,timeStart+ceil(samplingPeriodNum)*samplingPointsNumEachPeriod *timeSample]);
ylim([min(dataSample)-0.1,max(dataSample)+0.1]);
% xTick & xTickLabel
xTick = timeStart:timeSample*samplingPointsNumEachPeriod:timeStart+ceil(samplingPeriodNum)*samplingPointsNumEachPeriod *timeSample;
xTickLabel = cell(1,ceil(samplingPeriodNum));
for k = 1:length(xTick)
    if k==1
        xTickLabel{k} = sprintf('%.4f',xTick(k));
    else
        xTickLabel{k} = sprintf('[%d]%.4f',k-1,xTick(k));
    end
end
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
xlabel('Time (/$s$)','Interpreter','latex');
ylabel('Amplitude');

% 绘制频谱
figure('name','Fourier Spectrum (Non-Integate Periods of Signal)','NumberTitle','off');
stem(frequencySequency-frequencySample/2,abs(fftshift(fftDataSample)));
title('Fourier Spectrum (Non-Integate Periods of Signal)');
grid on;
xlim([min(frequencySequency-frequencySample/2),max(frequencySequency-frequencySample/2)]);
% xTick & xTickLabel
xTick2 = [-ceil(frequencySample/2):frequencySignal*2:-frequencySignal*2,-frequencySignal,...
    0,frequencySignal,frequencySignal*2:frequencySignal*2:ceil(frequencySample/2)];
xTickLabel2 = cell(1,length(xTick2));
for k = 1:length(xTick2)
  xTickLabel2{k} = num2str(xTick2(k));
end
set(gca, 'XTick', xTick2);set(gca, 'XTickLabel',xTickLabel2);
xlabel('Frequency (/$Hz$)','Interpreter','latex');
ylabel('Magnitude');

% 显示非整数周期信号及其Hilbert变换
figure('name','Non-Integate Periods of Signal and its Hilbert Transform','NumberTitle','off');
% 非整数周期信号
plot(timeSequence,dataSample,                       plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 期望Hilbert变换结果
expectedHilbert=-cos(2*pi*frequencySignal*timeSequence);
plot(timeSequence,expectedHilbert,plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 实际Hilbert变换结果
plot(timeSequence,imag(hilbert(dataSample)),        plotLineType,'Color','m','LineWidth',0.5,'MarkerSize',2);
title('Non-Integate Periods of Signal and its Hilbert Transform');
legend('Signal','Expected HT','Actual HT','Location','southoutside','Orientation','horizontal');
xlim([timeStart,timeStart+ceil(samplingPeriodNum)*samplingPointsNumEachPeriod *timeSample]);
ylim([min(imag(hilbert(dataSample)))-0.1,max(imag(hilbert(dataSample)))+0.1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
grid on;

