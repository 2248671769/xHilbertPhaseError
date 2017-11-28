% xAboutHilbertTransformDifferentAmplitude.m
% 徐文宇，20171102，周四
% 改变条纹信号的调制度，分析整数周期采样与非整数周期采样对Hilbert变换以及Hilbert域相位的影响
% ver：---
close all;clear;

%% @_@{设置基本参数}*************************************************
% 图像条纹参数
width=1280; height=800; period=128;

%% --幅度调制
% 阶跃调制幅度标记
stepFunctionModulateFlag=1;
% 对称连续弧段调制幅度标记
symmetricalArcModulateFlag=1;
% 非对称连续弧段调制幅度标记
asymmetricalArcModulateFlag=1;

% 信号范围
numOfPeriods=8;
startOfSignal=1;
endOfSignal=startOfSignal+period*numOfPeriods-1;
lengthOfSignal=endOfSignal-startOfSignal+1;
% xTick & xTickLabel
xTick=zeros(1,numOfPeriods+1);
xTickLabel=cell(1,numOfPeriods+1);
for xt=0:numOfPeriods
    xTick(xt+1)=floor(xt*period); xTickLabel{xt+1}=num2str(xTick(xt+1));
end
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
partNum=16;
xTickFourierSpectrum=zeros(1,partNum+1);
xTickLabelFourierSpectrum=cell(1,partNum+1);
for xt=0:partNum/2
    xTickFourierSpectrum(partNum/2+xt+1)=lengthOfSignal/2+1+xt+xt; xTickLabelFourierSpectrum{partNum/2+xt+1}=num2str(( xt+xt)*period/8);
    xTickFourierSpectrum(partNum/2-xt+1)=lengthOfSignal/2+1-xt-xt; xTickLabelFourierSpectrum{partNum/2-xt+1}=num2str((-xt-xt)*period/8);
end
xTickFourierSpectrum(partNum/2+1)=lengthOfSignal/2+1; xTickLabelFourierSpectrum{partNum/2+1}=num2str(0);

% 相位误差显示有效区间
upPhaseErrorBound=2; bottomPhaseErrorBound=-2;

% plot画线类型
plotLineType='';        % '' 实线
plotDottedLineType=':'; % ':'虚线

%% @_@{生成24幅全部条纹图像(单一幅度)}*******************************
%% --生成24幅全部条纹图像并从中抽取出数步相移条纹图像
% -单位半幅度条纹信号
moveNumAll=24;
fringeListAllUnit=cell(moveNumAll,1);
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllUnit{k} = floor(255*0.5*(cos(((0:lengthOfSignal-1)-sf)/period*2*pi)+1)/2);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractional=SelectNStepFring(fringeListAllUnit,moveNumPart);

%% --计算理想空域相位
wrappedPhaseAll=GetWrapPhase(fringeListAllUnit,moveNumAll);
AB=GetABNormal(fringeListAllUnit,wrappedPhaseAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractional=GetWrapPhase(fringeListFractional,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbert=HilbertPerRow(fringeListFractional,moveNumPart);
wrappedPhaseFractionalHilbert=GetWrapPhaseWithHilbert(fringeListFractionalHilbert,moveNumPart);

%% --显示图表
% 显示信号及其Hilbert变换
figure('name','Original Fringe','NumberTitle','off');
plot(fringeListFractional{2});hold on;
plot(imag(hilbert(fringeListFractional{1})),plotLineType,'MarkerSize',2);
title('Original Fringe and its Hilbert Transform');
legend('Original Fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);


%% @_@{阶跃式调制幅度}***********************************************
if stepFunctionModulateFlag==1
%% --调整所有信号幅度为左全右半阶跃状态
% 阶跃函数
filterStepAmplitude=[ones(1,lengthOfSignal/2) 2*ones(1,lengthOfSignal/2)];
fringeListAllStepAmplitude=cell(size(fringeListAllUnit));
fringeListFractionalStepAmplitude=cell(size(fringeListFractional));
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalStepAmplitude=cell(size(fringeListFractional));
for k=1:moveNumAll
    fringeListAllStepAmplitude{k}=floor((fringeListAllUnit{k}-255.0/4).*filterStepAmplitude+255.0/2);
end
for k=1:moveNumPart
    fringeListFractionalStepAmplitude{k}=floor((fringeListFractional{k}-255.0/4).*filterStepAmplitude+255.0/2);
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalStepAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterStepAmplitude;
end

%% --计算理想空域相位
wrappedPhaseAllStepAmplitude=GetWrapPhase(fringeListAllStepAmplitude,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalStepAmplitude=GetWrapPhase(fringeListFractionalStepAmplitude,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertStepAmplitude=HilbertPerRow(fringeListFractionalStepAmplitude,moveNumPart);
wrappedPhaseFractionalHilbertStepAmplitude=GetWrapPhaseWithHilbert(fringeListFractionalHilbertStepAmplitude,moveNumPart);

%% --显示阶跃函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示阶跃函数及其Hilbert变换
htStepAmplitudeFilter=imag(hilbert(filterStepAmplitude));
figure('name','Step Function and its Hilbert Transform','NumberTitle','off');
plot(filterStepAmplitude,  plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
plot(htStepAmplitudeFilter,plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Step Function and its Hilbert Transform');
legend('Step Function','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
% 显示阶跃函数的Fourier频谱
figure('name','Fourier Spectrum of Step Function','NumberTitle','off');
plot(abs(fftshift(fft(filterStepAmplitude))),plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Fourier Spectrum of Step Function');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);

% 显示信号及其Hilbert变换
figure('name','Amplitude Modulated by Step Function','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(fringeListFractionalStepAmplitude{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(expectedHilbertOfFringeListFractionalStepAmplitude{k}
    % 实际Hilbert变换结果
    plot(imag(hilbert(fringeListFractionalStepAmplitude{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    if k==moveNumPart
        legend('Original Fringe','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示Hilbert变换期望值与实际值的条纹误差曲线
figure('name','Fringe Error between Expected HT and Actual HT','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    plot(imag(hilbert(fringeListFractionalStepAmplitude{k}))-expectedHilbertOfFringeListFractionalStepAmplitude{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Modulated by Step Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalStepAmplitude-wrappedPhaseAllStepAmplitude,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertStepAmplitude-wrappedPhaseAllStepAmplitude,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Amplitude Modulated by Step Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/阶跃式Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------stepAmplitudeModulate-------------\n');
fprintf('        Mean of   Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('Max positive of   Space Phase Error: %+f\n',max( wrappedErrorSpace));
fprintf('Max negative of   Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('        RMSE of   Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);

end

%% @_@{对称连续弧段调制幅度}****************************************
if symmetricalArcModulateFlag==1
%% --以对称连续弧段函数调制所有信号
% 对称连续弧段函数
arcCenter=0.5;
filterSymmetricalArcAmplitude=1-2.5*((0:lengthOfSignal-1)/lengthOfSignal-arcCenter).^2;
fringeListAllSymmetricalArcAmplitude=cell(size(fringeListAllUnit));
fringeListFractionalSymmetricalArcAmplitude=cell(size(fringeListFractional));
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalSymmetricalArcAmplitude=cell(size(fringeListFractional));
for k=1:moveNumAll
    fringeListAllSymmetricalArcAmplitude{k}=floor(fringeListAllUnit{k}.*filterSymmetricalArcAmplitude);
end
for k=1:moveNumPart
    fringeListFractionalSymmetricalArcAmplitude{k}=floor(fringeListFractional{k}.*filterSymmetricalArcAmplitude);
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalSymmetricalArcAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterSymmetricalArcAmplitude;
end

%% --计算理想空域相位
wrappedPhaseAllSymmetricalArcAmplitude=GetWrapPhase(fringeListAllSymmetricalArcAmplitude,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalSymmetricalArcAmplitude=GetWrapPhase(fringeListFractionalSymmetricalArcAmplitude,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertSymmetricalArcAmplitude=HilbertPerRow(fringeListFractionalSymmetricalArcAmplitude,moveNumPart);
wrappedPhaseFractionalHilbertSymmetricalArcAmplitude=GetWrapPhaseWithHilbert(fringeListFractionalHilbertSymmetricalArcAmplitude,moveNumPart);

%% --显示对称连续弧段函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 组合显示原始信号、对称连续弧段函数、调制后的组合函数
figure('name','Symmetrical Arc Function','NumberTitle','off');
plot(filterSymmetricalArcAmplitude, plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Symmetrical Arc Function');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 组合显示原始信号、调制后的组合函数
figure('name','1/4 Step of Original Fringe & Modulated Signal by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListFractional{1},                       plotLineType,'LineWidth',1.5,'MarkerSize',2);hold on;
plot(fringeListFractionalSymmetricalArcAmplitude{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title('1/4 Step of Original Fringe & Modulated Signal by Symmetrical Arc Function');
legend('Fringe Signal','Modulated Signal','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 分别显示原始信号、调制后的组合函数、对称连续弧段函数的Fourier变换
figure('name','Fourier Transform of Fringe Signal, Modulated Signal, Symmetrical Arc Function','NumberTitle','off');
% 原始信号
subplot(3,1,1)
plot(abs(fftshift(fft(fringeListFractional{1}))),                       plotLineType,'LineWidth',1,'MarkerSize',2);
title('Fourier Spectrum of Fringe Signal');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
% 调制后的组合函数
subplot(3,1,2)
plot(abs(fftshift(fft(fringeListFractionalSymmetricalArcAmplitude{1}))),plotLineType,'LineWidth',1,'MarkerSize',2);
title('Fourier Spectrum of Modulated Signal');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
% 对称连续弧段函数
subplot(3,1,3)
plot(abs(fftshift(fft(filterSymmetricalArcAmplitude))),              plotLineType,'LineWidth',1,'MarkerSize',2);
title('Fourier Spectrum of Symmetrical Arc Function');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);

% 显示对称连续弧段函数及其Hilbert变换
htSymmetricalArcAmplitudeFilter=imag(hilbert(filterSymmetricalArcAmplitude));
figure('name','Symmetrical Arc Function and its Hilbert Transform','NumberTitle','off');
plot(filterSymmetricalArcAmplitude,                 plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
plot(htSymmetricalArcAmplitudeFilter,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
title('Symmetrical Arc Function and its Hilbert Transform');
legend('Symmetrical Arc Function ','HT','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
figure('name','Amplitude Modulated by Symmetrical Arc Function','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(fringeListFractionalSymmetricalArcAmplitude{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(expectedHilbertOfFringeListFractionalSymmetricalArcAmplitude{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 实际Hilbert变换结果
    plot(imag(hilbert(fringeListFractionalSymmetricalArcAmplitude{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    if k==moveNumPart
        legend('Original Fringe','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);ylim([-64-8 128]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示Hilbert变换期望值与实际值的条纹误差曲线
figure('name','Fringe Error between Expected HT and Actual HT','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    plot(imag(hilbert(fringeListFractionalSymmetricalArcAmplitude{k}))-expectedHilbertOfFringeListFractionalSymmetricalArcAmplitude{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Modulated by Symmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalSymmetricalArcAmplitude    -wrappedPhaseAllSymmetricalArcAmplitude,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertSymmetricalArcAmplitude-wrappedPhaseAllSymmetricalArcAmplitude,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Amplitude Modulated by Symmetrical Arc Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/阶跃式Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------symmetricalArcModulate-------------\n');
fprintf('        Mean of   Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('Max positive of   Space Phase Error: %+f\n',max( wrappedErrorSpace));
fprintf('Max negative of   Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('        RMSE of   Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);
end

%% @_@{非对称连续弧段调制幅度}**************************************
if asymmetricalArcModulateFlag==1
%% --以非对称连续弧段函数调制所有信号
% 非对称连续弧段函数
arcCenter=0.625;
filterAsymmetricalArcAmplitude=1-1.0*((0:lengthOfSignal-1)/lengthOfSignal-arcCenter).^2;
fringeListAllAsymmetricalArcAmplitude=cell(size(fringeListAllUnit));
fringeListFractionalAsymmetricalArcAmplitude=cell(size(fringeListFractional));
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalAsymmetricalArcAmplitude=cell(size(fringeListFractional));
for k=1:moveNumAll
    fringeListAllAsymmetricalArcAmplitude{k}=floor(fringeListAllUnit{k}.*filterAsymmetricalArcAmplitude);
end
for k=1:moveNumPart
    fringeListFractionalAsymmetricalArcAmplitude{k}=floor(fringeListFractional{k}.*filterAsymmetricalArcAmplitude);
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalAsymmetricalArcAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterAsymmetricalArcAmplitude;
end

%% --计算理想空域相位
wrappedPhaseAllAsymmetricalArcAmplitude=GetWrapPhase(fringeListAllAsymmetricalArcAmplitude,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalAsymmetricalArcAmplitude=GetWrapPhase(fringeListFractionalAsymmetricalArcAmplitude,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAsymmetricalArcAmplitude=HilbertPerRow(fringeListFractionalAsymmetricalArcAmplitude,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcAmplitude=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcAmplitude,moveNumPart);

%% --显示非对称连续弧段函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 组合显示原始信号、对称非连续弧段函数、调制后的组合函数
figure('name','Asymmetrical Arc Function','NumberTitle','off');
plot(filterAsymmetricalArcAmplitude, plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Asymmetrical Arc Function');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 组合显示原始信号、调制后的组合函数
figure('name','1/4 Step of Original Fringe & Modulated Signal by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListFractional{1},                        plotLineType,'LineWidth',1.5,'MarkerSize',2);hold on;
plot(fringeListFractionalAsymmetricalArcAmplitude{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title('1/4 Step of Original Fringe & Modulated Signal by Asymmetrical Arc Function');
legend('Fringe Signal','Modulated Signal','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
 
% 分别显示原始信号、调制后的组合函数、非对称连续弧段函数的Fourier变换
figure('name','Fourier Transform of Fringe Signal, Modulated Signal, Asymmetrical Arc Function','NumberTitle','off');
% 原始信号
subplot(3,1,1)
plot(abs(fftshift(fft(fringeListFractional{1}))),                        plotLineType,'LineWidth',1,'MarkerSize',2);
title('Fourier Spectrum of Fringe Signal');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
% 调制后的组合函数
subplot(3,1,2)
plot(abs(fftshift(fft(fringeListFractionalAsymmetricalArcAmplitude{1}))),plotLineType,'LineWidth',1,'MarkerSize',2);
title('Fourier Spectrum of Modulated Signal');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
% 非对称连续弧段函数
subplot(3,1,3)
plot(abs(fftshift(fft(filterAsymmetricalArcAmplitude))),              plotLineType,'LineWidth',1,'MarkerSize',2);
title('Fourier Spectrum of Asymmetrical Arc Function');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);

% 显示非对称连续弧段函数及其Hilbert变换
htAsymmetricalArcAmplitudeFilter=imag(hilbert(filterAsymmetricalArcAmplitude));
figure('name','Asymmetrical Arc Function and its Hilbert Transform','NumberTitle','off');
plot(filterAsymmetricalArcAmplitude,  plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
plot(htAsymmetricalArcAmplitudeFilter,plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Asymmetrical Arc Function and its Hilbert Transform');
legend('Asymmetrical Arc Function ','HT','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(fringeListFractionalAsymmetricalArcAmplitude{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(expectedHilbertOfFringeListFractionalAsymmetricalArcAmplitude{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 实际Hilbert变换结果
    plot(imag(hilbert(fringeListFractionalAsymmetricalArcAmplitude{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    if k==moveNumPart
        legend('Original Fringe','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);ylim([-64-8 128]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示Hilbert变换期望值与实际值的条纹误差曲线
figure('name','Fringe Error between Expected HT and Actual HT','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    plot(imag(hilbert(fringeListFractionalAsymmetricalArcAmplitude{k}))-expectedHilbertOfFringeListFractionalAsymmetricalArcAmplitude{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcAmplitude-wrappedPhaseAllAsymmetricalArcAmplitude,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertAsymmetricalArcAmplitude-wrappedPhaseAllAsymmetricalArcAmplitude,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Amplitude Modulated by Asymmetrical Arc Function)');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Location','NorthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/阶跃式Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------asymmetricalArcModulate-------------\n');
fprintf('        Mean of   Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('Max positive of   Space Phase Error: %+f\n',max( wrappedErrorSpace));
fprintf('Max negative of   Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('        RMSE of   Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);
end




