% xAboutHilbertTransformNonIntegateFrequency.m
% 徐文宇，20171106，周一
% 分析非整数周期采样对Hilbert变换以及Hilbert域相位的影响
% ver：---
close all;clear;

%% @_@{设置基本参数}****************************************************
% 图像条纹参数
width=1024; height=800; period=128;

%% --频率调制
% 整数周期标记
integatePeriodsFlag=1;
% 非整数周期标记
nonIntegatePeriodsFlag=1;

% 信号范围
numOfPeriods=8;
startOfSignal=1;
endOfSignal=startOfSignal+period*numOfPeriods-1;
lengthOfSignal=endOfSignal-startOfSignal+1;
% xTick & xTickLabel
xTick=zeros(1,numOfPeriods+2);
xTickLabel=cell(1,numOfPeriods+2);
for xt=0:numOfPeriods
    xTick(xt+1)=floor(xt*period); xTickLabel{xt+1}=num2str(xTick(xt+1));
end
xTick(end-1)=lengthOfSignal-1-period/2; xTickLabel{end-1}=num2str(lengthOfSignal-1-period/2);
xTick(end)=lengthOfSignal-1; xTickLabel{end}=num2str(lengthOfSignal-1);
% yTick & yTickLabel
yTickNum=8;
yTick=zeros(1,yTickNum+1);
yTickLabel=cell(1,yTickNum+1);
yTick(yTickNum/2+1)=0;
yTickLabel{yTickNum/2+1}='0';
for xt=1:yTickNum/2
    yTick(yTickNum/2+1+xt)=floor( xt*256/(yTickNum/2)); yTickLabel{yTickNum/2+1+xt}=num2str(yTick(yTickNum/2+1+xt)); 
    yTick(yTickNum/2+1-xt)=floor(-xt*256/(yTickNum/2)); yTickLabel{yTickNum/2+1-xt}=num2str(yTick(yTickNum/2+1-xt));
end
% xTickPart & xTickLabelPart for Fourier Spectrum
partNum=30;
xTickFourierSpectrum=zeros(1,partNum+1);
xTickLabelFourierSpectrum=cell(1,partNum+1);
for xt=0:partNum/2
    xTickFourierSpectrum(partNum/2+xt+1)=lengthOfSignal/2+1+xt+xt; xTickLabelFourierSpectrum{partNum/2+xt+1}=num2str( xt+xt);
    xTickFourierSpectrum(partNum/2-xt+1)=lengthOfSignal/2+1-xt-xt; xTickLabelFourierSpectrum{partNum/2-xt+1}=num2str(-xt-xt);
end
xTickFourierSpectrum(partNum/2+1)=lengthOfSignal/2+1; xTickLabelFourierSpectrum{partNum/2+1}=num2str(0);

% 相位误差显示有效区间
upPhaseErrorBound=2; bottomPhaseErrorBound=-2;

% plot画线类型
plotLineType='';        % '' 实线
plotDottedLineType=':'; % ':'虚线

%% @_@{整数周期频率}****************************************************
% 整数周期标记
if integatePeriodsFlag==1
%% --生成整数周期单位半幅度条纹信号
moveNumAll=24;
fringeListAllIntegatePeriods=cell(moveNumAll,1);
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllIntegatePeriods{k} = floor(255*0.5*(cos(((0:lengthOfSignal-1)-sf)/period*2*pi)+1)/2);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalIntegatePeriods=SelectNStepFring(fringeListAllIntegatePeriods,moveNumPart);
expectedHilbertOfFringeListFractionalIntegatePeriods=cell(size(fringeListFractionalIntegatePeriods));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalIntegatePeriods{k}=255.0/2*(sin(((1:lengthOfSignal)-sf)/period*2*pi))/2;
end

%% --计算理想空域相位
wrappedPhaseAllIntegatePeriods=GetWrapPhase(fringeListAllIntegatePeriods,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalIntegatePeriods=GetWrapPhase(fringeListFractionalIntegatePeriods,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertIntegatePeriods=HilbertPerRow(fringeListFractionalIntegatePeriods,moveNumPart);
wrappedPhaseFractionalHilbertIntegatePeriods=GetWrapPhaseWithHilbert(fringeListFractionalHilbertIntegatePeriods,moveNumPart);

%% --显示整数周期数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示整数周期信号及其Hilbert变换
figure('name','Integate Periods of Signal and its Hilbert Transform','NumberTitle','off');
% 整数周期信号
plot(fringeListFractionalIntegatePeriods{1},                       plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 期望Hilbert变换结果
plot(expectedHilbertOfFringeListFractionalIntegatePeriods{1},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 实际Hilbert变换结果
plot(imag(hilbert(fringeListFractionalIntegatePeriods{1})),        plotLineType,'Color','m','LineWidth',0.5,'MarkerSize',2);
title('Integate Periods of Signal and its Hilbert Transform');
legend('Signal','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 显示整数周期信号的Fourier频谱
figure('name','Fourier Spectrum of Integate Periods of Signal','NumberTitle','off');
plot(abs(fftshift(fft(fringeListFractionalIntegatePeriods{1}))),plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Fourier Spectrum of Integate Periods of Signal');
% xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);grid on;
xlim([0,lengthOfSignal-1]);grid on;
% set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
figure('name','Integate Periods of Original Fringe','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(fringeListFractionalIntegatePeriods{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(expectedHilbertOfFringeListFractionalIntegatePeriods{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 实际Hilbert变换结果
    plot(imag(hilbert(fringeListFractionalIntegatePeriods{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    if k==moveNumPart
        legend('Original Fringe','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);ylim([-96 160]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Integate Periods)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalIntegatePeriods-wrappedPhaseAllIntegatePeriods,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertIntegatePeriods-wrappedPhaseAllIntegatePeriods,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'Color','g','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Frequency Changed by Step Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------integatePeriods-------------\n');
fprintf('        Mean of   Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('Max positive of   Space Phase Error: %+f\n',max( wrappedErrorSpace));
fprintf('Max negative of   Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('        RMSE of   Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);

end


%% @_@{非整数周期}******************************************************
% 非整数周期标记
if nonIntegatePeriodsFlag==1
%% --生成非整数周期单位半幅度条纹信号
moveNumAll=24;
fringeListAllNonIntegatePeriods=cell(moveNumAll,1);
for k=1:moveNumAll
    fringeListAllNonIntegatePeriods{k}=fringeListAllIntegatePeriods{k}(startOfSignal:endOfSignal-period/2);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalNonIntegatePeriods=SelectNStepFring(fringeListAllNonIntegatePeriods,moveNumPart);
expectedHilbertOfFringeListFractionalNonIntegatePeriods=cell(size(fringeListFractionalIntegatePeriods));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalNonIntegatePeriods{k}=255.0/2*(sin(((1:lengthOfSignal-period/2)-sf)/period*2*pi))/2;
end

%% --计算理想空域相位
wrappedPhaseAllNonIntegatePeriods=GetWrapPhase(fringeListAllNonIntegatePeriods,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalNonIntegatePeriods=GetWrapPhase(fringeListFractionalNonIntegatePeriods,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertNonIntegatePeriods=HilbertPerRow(fringeListFractionalNonIntegatePeriods,moveNumPart);
wrappedPhaseFractionalHilbertNonIntegatePeriods=GetWrapPhaseWithHilbert(fringeListFractionalHilbertNonIntegatePeriods,moveNumPart);

%% --显示非整数周期数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示非整数周期信号及其Hilbert变换
figure('name','Non-Integate Periods of Signal and its Hilbert Transform','NumberTitle','off');
% 非整数周期信号
plot(fringeListFractionalNonIntegatePeriods{1},                       plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 期望Hilbert变换结果
plot(expectedHilbertOfFringeListFractionalNonIntegatePeriods{1},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 实际Hilbert变换结果
plot(imag(hilbert(fringeListFractionalNonIntegatePeriods{1})),        plotLineType,'Color','m','LineWidth',0.5,'MarkerSize',2);
title('Non-Integate Periods of Signal and its Hilbert Transform');
legend('Signal','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 显示整数周期信号的Fourier频谱
figure('name','Fourier Spectrum of Non-Integate Periods of Signal','NumberTitle','off');
plot(abs(fftshift(fft(fringeListFractionalNonIntegatePeriods{1}))),plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Fourier Spectrum of Non-Integate Periods of Signal');
% xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);
grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
figure('name','Non-Integate Periods of Original Fringe','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(fringeListFractionalNonIntegatePeriods{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(expectedHilbertOfFringeListFractionalNonIntegatePeriods{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 实际Hilbert变换结果
    plot(imag(hilbert(fringeListFractionalNonIntegatePeriods{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    if k==moveNumPart
        legend('Original Fringe','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);ylim([-96 160]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示Hilbert变换期望值与实际值的条纹误差曲线
figure('name','Fringe Error between Expected HT and Actual HT','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    plot(imag(hilbert(fringeListFractionalNonIntegatePeriods{k}))-expectedHilbertOfFringeListFractionalNonIntegatePeriods{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Non-Integate Periods)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalNonIntegatePeriods    -wrappedPhaseAllNonIntegatePeriods,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 非整数周期的Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertNonIntegatePeriods-wrappedPhaseAllNonIntegatePeriods,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Non-Integate Periods)');
legend('Space Phase Error','HT Phase Error','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------nonIntegatePeriods-------------\n');
fprintf('        Mean of   Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('Max positive of   Space Phase Error: %+f\n',max( wrappedErrorSpace));
fprintf('Max negative of   Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('        RMSE of   Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);

end

