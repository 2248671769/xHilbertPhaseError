% xGammaPhaseError.m
% 徐文宇，20171114，周四
% 模拟gamma效应，改变条纹信号的调制度，分析条纹调制度归一化样对Hilbert变换以及Hilbert域相位的影响
% ver：---
close all;clear;

%% @_@{设置基本参数}-☉☉☉⊿⊿⊿⊿
% 图像条纹参数
width=1024; height=800; period=128;

%% --调制标记
% 对非对称连续调制频率信号[完整周期]标记,归一化: S3/B
asymmetricalArcFrequencyModulateNormWithDividingBFlag=1;
% % 对称连续弧段调制幅度标记
symmetricalArcModulateFlag=0;
% 非对称连续弧段调制幅度Gamma非线性响应条纹标记,归一化：(S2-A)/B
asymmetricalArcModulateNormWithABFlag=0;
% 非对称连续弧段调制幅度Gamma非线性响应条纹标记,归一化：S2/B
asymmetricalArcModulateNormWithDividingBFlag=0;
% 均值为四分幅度的半幅度Gamma非线性响应条纹标记,归一化：(S1-A)/B}
halfAmplitudeNormWithABFlag=0;
% 均值为四分幅度的半幅度Gamma非线性响应条纹标记,归一化：S1-A
halfAmplitudeNormWithSubtractingAFlag=0;
% 均值为四分幅度的半幅度Gamma非线性响应条纹标记,归一化：S1/B
halfAmplitudeNormWithDividingBFlag=0;
% 均值为0的Gamma非线性响应条纹条纹标记,归一化：S0/B
unit0NormWithDividingBFlag=0;
 

% 信号范围
numOfPeriods=8;
startOfSignal=1;
endOfSignal=startOfSignal+period*numOfPeriods-1;
lengthOfSignal=endOfSignal-startOfSignal+1;

% 横坐标序列
horizontalIndex=0:lengthOfSignal-1;

% xTick & xTickLabel
if mod(lengthOfSignal,period)==0
    xTick=zeros(1,numOfPeriods+1);
    xTickLabel=cell(1,numOfPeriods+1);
    for xt=0:numOfPeriods
        xTick(xt+1)=floor(xt*period); xTickLabel{xt+1}=num2str(xTick(xt+1));
    end
    xTick(end)=lengthOfSignal-1; xTickLabel{end}=num2str(lengthOfSignal-1);
else
    xTick=zeros(1,numOfPeriods+2);
    xTickLabel=cell(1,numOfPeriods+2);
    for xt=0:numOfPeriods
        xTick(xt+1)=floor(xt*period); xTickLabel{xt+1}=num2str(xTick(xt+1));
    end
    xTick(end)=lengthOfSignal-1; xTickLabel{end}=num2str(lengthOfSignal-1);
end
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

%% --Gamma曲线
gamma=3.0;
filterGamma=255*power((0:255)/255,gamma);
figure('name','Gamma Curve','NumberTitle','off');
plot(filterGamma,plotLineType,'LineWidth',1,'MarkerSize',2);hold on;
plot((0:256-1),'g-.','LineWidth',0.5);
title(sprintf('Gamma Curve ($\\gamma$=%1.2f)',gamma),'Interpreter','latex');
xlim([0,256]);ylim([0,256]);grid on;
set(gca, 'XTick', yTick);set(gca, 'XTickLabel',yTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);


%% @_@{生成24幅全部条纹图像(单一幅度)}☉☉☉☉☉☉☉☉☉☉
%% --生成24幅全部条纹图像并从中抽取出数步相移条纹图像
% -单位半幅度条纹信号
moveNumAll=24;
fringeListAllUnit=cell(moveNumAll,1);
fringeListAllUnitGamma=cell(moveNumAll,1);
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllUnit{k} = floor(255*0.5*(cos(((0:lengthOfSignal-1)-sf)/period*2*pi)+1)/2);
    fringeListAllUnitGamma{k} =floor(filterGamma(floor(fringeListAllUnit{k})+1));
%     255*0.5*power(floor((255*0.5*(cos(((0:lengthOfSignal-1)-sf)/period*2*pi)+1)/2))/255,gamma);   
end

% -显示gamma效应前后的条纹信号
figure('name','Gamma Effect for Unit Fringes','NumberTitle','off');
plot(fringeListAllUnit{1});hold on;
plot(fringeListAllUnitGamma{1},plotLineType,'MarkerSize',2);
title(sprintf('Gamma Effect ($\\gamma$=%1.2f) for Unit Fringes',gamma),'Interpreter','latex');
legend('Original Fringe','Gamma Effect','Location','NorthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% -计算理想空域相位
wrappedPhaseAllUnit=GetWrapPhase(fringeListAllUnit,moveNumAll);

% -计算理想空域相位[Gamma效应)
wrappedPhaseAllGamma=GetWrapPhase(fringeListAllUnitGamma,moveNumAll);

% 抽取出数步相移条纹图像
moveNumPart=3;
fringeListFractional=SelectNStepFring(fringeListAllUnit,moveNumPart);
fringeListGammaFractional=SelectNStepFring(fringeListAllUnitGamma,moveNumPart);

%% --计算理想空域相位
% wrappedPhaseAll0=GetWrapPhase(fringeListAllUnit,moveNumAll);
% AB=GetABNormal(fringeListAllUnit,wrappedPhaseAllGamma);

%% --计算数步相移条纹的空域相位
wrappedPhaseGammaFractional=GetWrapPhase(fringeListGammaFractional,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListGammaFractionalHilbert=HilbertPerRow(fringeListGammaFractional,moveNumPart);
wrappedPhaseGammaFractionalHilbert=GetWrapPhaseWithHilbert(fringeListGammaFractionalHilbert,moveNumPart);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma Effect) for Unit Fringes','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseGammaFractional    -wrappedPhaseAllUnit,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseGammaFractionalHilbert-wrappedPhaseAllUnit,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,plotDottedLineType,'LineWidth',1.5,'MarkerSize',2);hold on;
% 均值相位误差
plot((wrappedErrorSpace+wrappedErrorHT)/2,plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Gamma Phase Error ($\\gamma$=%1.2f) for Unit Fringes',gamma),'Interpreter','latex');
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);


%% @_@{对称连续弧段调制幅度}☉☉☉☉☉☉☉☉☉☉☉☉☉
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
    fringeListAllSymmetricalArcAmplitude{k}=floor((fringeListAllUnit{k}-255.0/4).*filterSymmetricalArcAmplitude+255.0/2);
end
for k=1:moveNumPart
    fringeListFractionalSymmetricalArcAmplitude{k}=floor(fringeListFractional{k}.*filterSymmetricalArcAmplitude);
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalSymmetricalArcAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterSymmetricalArcAmplitude;
end

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
    xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
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
plot(extractValidPhaseErrorWithBounds(wrappedPhaseFractional                              -wrappedPhaseAllGamma,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertSymmetricalArcAmplitude-wrappedPhaseAllGamma,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Amplitude Modulated by Symmetrical Arc Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/阶跃式Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------symmetricalArcModulate-------------\n');
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractional-wrappedPhaseAllGamma,upPhaseErrorBound,bottomPhaseErrorBound);
fprintf('          Mean of Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('  Max positive of Space Phase Error: %+f\n',max( phasepedErrorSpace));
fprintf('  Max negative of Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('          RMSE of Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertSymmetricalArcAmplitude-wrappedPhaseAllGamma,upPhaseErrorBound,bottomPhaseErrorBound);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( phasepedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);
end

%% @_@{非对称连续弧段调制幅度Gamma非线性响应条纹标记,归一化：(S2-A)/B}☉☉☉☉☉☉☉☉☉☉☉
if asymmetricalArcModulateNormWithABFlag==1
%% --以非对称连续弧段函数调制所有信号
% 非对称连续弧段函数
arcCenter=0.625;
filterAsymmetricalArcAmplitude=1-1.0*((0:lengthOfSignal-1)/lengthOfSignal-arcCenter).^2;
fringeListAllAsymmetricalArcAmplitude=cell(moveNumAll,1);
fringeListAllAsymmetricalArcAmplitudeGamma=cell(moveNumAll,1);
% 期望的Hilbert变换
% expectedHilbertOfFringeListFractionalAsymmetricalArcAmplitude=cell(moveNumPart,1);
for k=1:moveNumAll
    fringeListAllAsymmetricalArcAmplitude{k}=floor(fringeListAllUnit{k}.*filterAsymmetricalArcAmplitude);
    fringeListAllAsymmetricalArcAmplitudeGamma{k}=floor(filterGamma(floor(fringeListAllAsymmetricalArcAmplitude{k})+1));
end
% for k=1:moveNumPart
%     fringeListFractionalAsymmetricalArcAmplitude{k}=floor(fringeListFractional{k}.*filterAsymmetricalArcAmplitude);
%     sf=-period*(k-1)/moveNumPart;
% %     expectedHilbertOfFringeListFractionalAsymmetricalArcAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterAsymmetricalArcAmplitude;
% end

% 抽取出数步相移条纹图像
fringeListFractionalAsymmetricalArcAmplitude     =SelectNStepFring(fringeListAllAsymmetricalArcAmplitude     ,moveNumPart);
fringeListFractionalAsymmetricalArcAmplitudeGamma=SelectNStepFring(fringeListAllAsymmetricalArcAmplitudeGamma,moveNumPart);

%% --计算理想空域相位[Gamma]
wrappedPhaseAllAsymmetricalArcAmplitudeGamma=GetWrapPhase(fringeListAllAsymmetricalArcAmplitudeGamma,moveNumAll);

%% --计算数步相移条纹的空域相位[Gamma]
wrappedPhaseFractionalAsymmetricalArcAmplitudeGamma=GetWrapPhase(fringeListFractionalAsymmetricalArcAmplitudeGamma,moveNumPart);

%% --计算数步相移条纹[Gamma]的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAsymmetricalArcAmplitudeGamma=HilbertPerRow(fringeListFractionalAsymmetricalArcAmplitudeGamma,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcAmplitudeGamma=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcAmplitudeGamma,moveNumPart);

%% --显示非对称连续弧段函数、数步相移条纹信号[Gamma]及其Hilbert变换、相位误差图表
% 组合显示原始信号、对称非连续弧段函数、调制后的组合函数
figure('name','Asymmetrical Arc Function','NumberTitle','off');
plot(filterAsymmetricalArcAmplitude, plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Asymmetrical Arc Function');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 组合显示原始信号、调制后的组合函数
figure('name',sprintf('1/%d Step of Original Fringe & Modulated Signal by Asymmetrical Arc Function [Gamma]',moveNumPart),'NumberTitle','off');
plot(fringeListFractional{1},                             plotLineType,'LineWidth',1.5,'MarkerSize',2);hold on;
plot(fringeListFractionalAsymmetricalArcAmplitude{1},     plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
plot(fringeListFractionalAsymmetricalArcAmplitudeGamma{1},plotLineType,'Color','m','LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringes [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
legend('Fringe Signal','Modulated Signal','Modulated Signal [Gamma]','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% [Gamma without Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma] (Amplitude Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcAmplitudeGamma-wrappedPhaseAllAsymmetricalArcAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertAsymmetricalArcAmplitudeGamma-wrappedPhaseAllAsymmetricalArcAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
% 平均相位误差
wrappedErrorMean=(wrappedErrorSpace+wrappedErrorHT)/2;
plot(wrappedErrorMean,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','NorthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% --对信号[Gamma]归一化[Normalization]处理(S-A)/B及归一化后的理想空域相位、Hilbert域相位、平均相位
% 理想信号[Gamma]
ABAllAmplitudeModulatedGamma=GetABNormal(fringeListAllAsymmetricalArcAmplitudeGamma,wrappedPhaseAllAsymmetricalArcAmplitudeGamma);
fringeListAllAsymmetricalArcAmplitudeGammaNorm=NormalizeFringe(fringeListAllAsymmetricalArcAmplitudeGamma,moveNumAll,ABAllAmplitudeModulatedGamma);
wrappedPhaseAllAsymmetricalArcAmplitudeGammaNorm=GetWrapPhase(fringeListAllAsymmetricalArcAmplitudeGammaNorm,moveNumAll);
% 数步信号[Gamma]
ABFractionalAmplitudeModulatedGamma=GetABNormal(fringeListFractionalAsymmetricalArcAmplitudeGamma,wrappedPhaseFractionalAsymmetricalArcAmplitudeGamma);
fringeListFractionalAsymmetricalArcAmplitudeGammaNorm=NormalizeFringe(fringeListFractionalAsymmetricalArcAmplitudeGamma,moveNumPart,ABFractionalAmplitudeModulatedGamma);
wrappedPhaseFractionalAsymmetricalArcAmplitudeGammaNorm=GetWrapPhase(fringeListFractionalAsymmetricalArcAmplitudeGammaNorm,moveNumPart);
% 数步信号Hilbert[Gamma]
fringeListFractionalHilbertAsymmetricalArcAmplitudeGammaNorm=HilbertPerRow(fringeListFractionalAsymmetricalArcAmplitudeGammaNorm,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcAmplitudeGammaNorm=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcAmplitudeGammaNorm,moveNumPart);

% 显示归一化[Normalization]的信号
figure('name',sprintf('1/%d Step of Fringe Modulated by Asymmetrical Arc Function [Gamma/Normalization: (S-A)/B]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalAsymmetricalArcAmplitudeGammaNorm{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringe Modulated by Asymmetrical Arc Function [Gamma/Normalization: (S-A)/B]',moveNumPart));
xlim([0,lengthOfSignal-1]);ylim([-1,1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma with Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma/Normalization: (S-A)/B] (Amplitude Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpaceNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcAmplitudeGammaNorm-wrappedPhaseAllAsymmetricalArcAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHTNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertAsymmetricalArcAmplitudeGammaNorm-wrappedPhaseAllAsymmetricalArcAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 平均相位误差
wrappedErrorMeanNorm=(wrappedErrorSpaceNorm+wrappedErrorHTNorm)/2;
plot(wrappedErrorMeanNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error [Gamma ($\\gamma$=%1.2f)/Normalization: (S-A)/B]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','NorthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

end

%% @_@{非对称连续弧段调制幅度Gamma非线性响应条纹标记,归一化：S2/B}☉☉☉☉☉☉☉☉☉☉☉
if asymmetricalArcModulateNormWithDividingBFlag==1
%% --以非对称连续弧段函数调制所有信号
% 非对称连续弧段函数
arcCenter=0.625;
filterAsymmetricalArcAmplitude=1-1.0*((0:lengthOfSignal-1)/lengthOfSignal-arcCenter).^2;
fringeListAllAsymmetricalArcAmplitude=cell(moveNumAll,1);
fringeListAllAsymmetricalArcAmplitudeGamma=cell(moveNumAll,1);
% 期望的Hilbert变换
% expectedHilbertOfFringeListFractionalAsymmetricalArcAmplitude=cell(moveNumPart,1);
for k=1:moveNumAll
    fringeListAllAsymmetricalArcAmplitude{k}=floor(fringeListAllUnit{k}.*filterAsymmetricalArcAmplitude);
    fringeListAllAsymmetricalArcAmplitudeGamma{k}=floor(filterGamma(floor(fringeListAllAsymmetricalArcAmplitude{k})+1));
end
% for k=1:moveNumPart
%     fringeListFractionalAsymmetricalArcAmplitude{k}=floor(fringeListFractional{k}.*filterAsymmetricalArcAmplitude);
%     sf=-period*(k-1)/moveNumPart;
% %     expectedHilbertOfFringeListFractionalAsymmetricalArcAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterAsymmetricalArcAmplitude;
% end

% 抽取出数步相移条纹图像
fringeListFractionalAsymmetricalArcAmplitude     =SelectNStepFring(fringeListAllAsymmetricalArcAmplitude     ,moveNumPart);
fringeListFractionalAsymmetricalArcAmplitudeGamma=SelectNStepFring(fringeListAllAsymmetricalArcAmplitudeGamma,moveNumPart);

%% --计算理想空域相位[Gamma]
wrappedPhaseAllAsymmetricalArcAmplitudeGamma=GetWrapPhase(fringeListAllAsymmetricalArcAmplitudeGamma,moveNumAll);

%% --计算数步相移条纹的空域相位[Gamma]
wrappedPhaseFractionalAsymmetricalArcAmplitudeGamma=GetWrapPhase(fringeListFractionalAsymmetricalArcAmplitudeGamma,moveNumPart);

%% --计算数步相移条纹[Gamma]的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAsymmetricalArcAmplitudeGamma=HilbertPerRow(fringeListFractionalAsymmetricalArcAmplitudeGamma,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcAmplitudeGamma=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcAmplitudeGamma,moveNumPart);

%% --显示非对称连续弧段函数、数步相移条纹信号[Gamma]及其Hilbert变换、相位误差图表
% 组合显示原始信号、对称非连续弧段函数、调制后的组合函数
figure('name','Asymmetrical Arc Function','NumberTitle','off');
plot(filterAsymmetricalArcAmplitude, plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Asymmetrical Arc Function');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 组合显示原始信号、调制后的组合函数
figure('name',sprintf('1/%d Step of Original Fringe & Modulated Signal by Asymmetrical Arc Function [Gamma]',moveNumPart),'NumberTitle','off');
plot(fringeListFractional{1},                             plotLineType,'LineWidth',1.5,'MarkerSize',2);hold on;
plot(fringeListFractionalAsymmetricalArcAmplitude{1},     plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
plot(fringeListFractionalAsymmetricalArcAmplitudeGamma{1},plotLineType,'Color','m','LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringes [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
legend('Fringe Signal','Modulated Signal','Modulated Signal [Gamma]','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% [Gamma without Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma] (Amplitude Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcAmplitudeGamma-wrappedPhaseAllAsymmetricalArcAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertAsymmetricalArcAmplitudeGamma-wrappedPhaseAllAsymmetricalArcAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
% 平均相位误差
wrappedErrorMean=(wrappedErrorSpace+wrappedErrorHT)/2;
plot(wrappedErrorMean,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','NorthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% --对信号[Gamma]归一化[Normalization]处理S/B及归一化后的理想空域相位、Hilbert域相位、平均相位
% 理想信号[Gamma]
ABAllAmplitudeModulatedGamma=GetABNormal(fringeListAllAsymmetricalArcAmplitudeGamma,wrappedPhaseAllAsymmetricalArcAmplitudeGamma);
% fringeListAllAsymmetricalArcAmplitudeGammaNorm=NormalizeFringe(fringeListAllAsymmetricalArcAmplitudeGamma,moveNumAll,ABAllAmplitudeModulatedGamma);
for k=1:moveNumAll
    fringeListAllAsymmetricalArcAmplitudeGammaNorm{k}=fringeListAllAsymmetricalArcAmplitudeGamma{k}./ABAllAmplitudeModulatedGamma{2};
end
wrappedPhaseAllAsymmetricalArcAmplitudeGammaNorm=GetWrapPhase(fringeListAllAsymmetricalArcAmplitudeGammaNorm,moveNumAll);
% 数步信号[Gamma]
ABFractionalAmplitudeModulatedGamma=GetABNormal(fringeListFractionalAsymmetricalArcAmplitudeGamma,wrappedPhaseFractionalAsymmetricalArcAmplitudeGamma);
% fringeListFractionalAsymmetricalArcAmplitudeGammaNorm=NormalizeFringe(fringeListFractionalAsymmetricalArcAmplitudeGamma,moveNumPart,ABFractionalAmplitudeModulatedGamma);
for k=1:moveNumPart
    fringeListFractionalAsymmetricalArcAmplitudeGammaNorm{k}=fringeListFractionalAsymmetricalArcAmplitudeGammaNorm{k}./ABFractionalAmplitudeModulatedGamma{2};
end
wrappedPhaseFractionalAsymmetricalArcAmplitudeGammaNorm=GetWrapPhase(fringeListFractionalAsymmetricalArcAmplitudeGammaNorm,moveNumPart);
% 数步信号Hilbert[Gamma]
fringeListFractionalHilbertAsymmetricalArcAmplitudeGammaNorm=HilbertPerRow(fringeListFractionalAsymmetricalArcAmplitudeGammaNorm,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcAmplitudeGammaNorm=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcAmplitudeGammaNorm,moveNumPart);

% 显示归一化[Normalization]的信号
figure('name',sprintf('1/%d Step of Fringe Modulated by Asymmetrical Arc Function [Gamma/Normalization: S/B]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalAsymmetricalArcAmplitudeGammaNorm{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringe Modulated by Asymmetrical Arc Function [Gamma/Normalization: S/B]',moveNumPart));
xlim([0,lengthOfSignal-1]);ylim([-1,1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma with Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma/Normalization: S/B] (Amplitude Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpaceNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcAmplitudeGammaNorm-wrappedPhaseAllAsymmetricalArcAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHTNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertAsymmetricalArcAmplitudeGammaNorm-wrappedPhaseAllAsymmetricalArcAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 平均相位误差
wrappedErrorMeanNorm=(wrappedErrorSpaceNorm+wrappedErrorHTNorm)/2;
plot(wrappedErrorMeanNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error [Gamma ($\\gamma$=%1.2f)/Normalization: S/B]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','NorthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

end

%% @_@{均值为四分幅度的半幅度Gamma非线性响应条纹标记,归一化：(S1-A)/B}☉☉☉☉☉☉☉☉☉☉☉☉
if halfAmplitudeNormWithABFlag==1
%% --以非对称连续弧段函数调制所有信号
fringeListAllHalfAmplitude=cell(moveNumAll,1);
fringeListAllHalfAmplitudeGamma=cell(moveNumAll,1);
% 期望的Hilbert变换
% expectedHilbertOfFringeListFractionalHalfAmplitude=cell(moveNumPart,1);
for k=1:moveNumAll
    fringeListAllHalfAmplitude{k}=floor(fringeListAllUnit{k});
    fringeListAllHalfAmplitudeGamma{k}=floor(filterGamma(floor(fringeListAllHalfAmplitude{k})+1));
end
% for k=1:moveNumPart
%     fringeListFractionalHalfAmplitude{k}=floor(fringeListFractional{k}.*filterHalfAmplitude);
%     sf=-period*(k-1)/moveNumPart;
% %     expectedHilbertOfFringeListFractionalHalfAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterHalfAmplitude;
% end

% 抽取出数步相移条纹图像
fringeListFractionalHalfAmplitude     =SelectNStepFring(fringeListAllHalfAmplitude     ,moveNumPart);
fringeListFractionalHalfAmplitudeGamma=SelectNStepFring(fringeListAllHalfAmplitudeGamma,moveNumPart);

%% --计算理想空域相位[Gamma]
wrappedPhaseAllHalfAmplitudeGamma=GetWrapPhase(fringeListAllHalfAmplitudeGamma,moveNumAll);

%% --计算数步相移条纹的空域相位[Gamma]
wrappedPhaseFractionalHalfAmplitudeGamma=GetWrapPhase(fringeListFractionalHalfAmplitudeGamma,moveNumPart);

%% --计算数步相移条纹[Gamma]的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertHalfAmplitudeGamma=HilbertPerRow(fringeListFractionalHalfAmplitudeGamma,moveNumPart);
wrappedPhaseFractionalHilbertHalfAmplitudeGamma=GetWrapPhaseWithHilbert(fringeListFractionalHilbertHalfAmplitudeGamma,moveNumPart);

%% --显示数步相移条纹信号[Gamma]及其Hilbert变换、相位误差图表
% 组合显示原始信号、调制后的组合函数
figure('name',sprintf('1/%d Step of Original Fringe (Half Amplitude) [Gamma]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalHalfAmplitude{1},     plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
plot(fringeListFractionalHalfAmplitudeGamma{1},plotLineType,'Color','m','LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringes (Half Amplitude) [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
legend('Fringe Signal','Fringe Signal [Gamma]','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示理想折叠相位
figure('name','Wrapped Phase (Half Amplitude) [Gamma/Normalization: (S-A)/B]','NumberTitle','off');
plot(wrappedPhaseAllHalfAmplitudeGamma,plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Wrapped Phase (Unit0 Amplitude) [Gamma/Normalization: (S-A)/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma without Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma]','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHalfAmplitudeGamma-wrappedPhaseAllHalfAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertHalfAmplitudeGamma-wrappedPhaseAllHalfAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
% 平均相位误差
wrappedErrorMean=(wrappedErrorSpace+wrappedErrorHT)/2;
plot(wrappedErrorMean,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Half Amplitude) [Gamma ($\\gamma$=%1.2f)]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% 对信号[Gamma]归一化[Normalization]处理(S-A)/B及归一化后的理想空域相位、Hilbert域相位、平均相位
% 理想信号[Gamma]
ABAllAmplitudeModulatedGamma=GetABNormal(fringeListAllHalfAmplitudeGamma,wrappedPhaseAllHalfAmplitudeGamma);
% 显示背景光强A
figure('name','Background Illumination A (Unit0 Amplitude) [Gamma/Normalization: (S-A)/B]','NumberTitle','off');
plot(ABAllAmplitudeModulatedGamma{1},plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Background Illumination A (Unit0 Amplitude) [Gamma/Normalization: (S-A)/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
% 显示调制度B
figure('name','Amplitude B (Unit0 Amplitude) [Gamma/Normalization: (S-A)/B]','NumberTitle','off');
plot(ABAllAmplitudeModulatedGamma{2},plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Amplitude B (Unit0 Amplitude) [Gamma/Normalization: (S-A)/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
fringeListAllHalfAmplitudeGammaNorm=cell(size(fringeListAllHalfAmplitudeGamma));
for k=1:moveNumAll
    fringeListAllHalfAmplitudeGammaNorm{k}=(fringeListAllHalfAmplitudeGamma{k}-ABAllAmplitudeModulatedGamma{1})./ABAllAmplitudeModulatedGamma{2};
end
% fringeListAllHalfAmplitudeGammaNorm=NormalizeFringe(fringeListAllHalfAmplitudeGamma,moveNumAll,ABAllAmplitudeModulatedGamma);
wrappedPhaseAllHalfAmplitudeGammaNorm=GetWrapPhase(fringeListAllHalfAmplitudeGammaNorm,moveNumAll);
% 数步信号[Gamma]
ABFractionalAmplitudeModulatedGamma=GetABNormal(fringeListFractionalHalfAmplitudeGamma,wrappedPhaseFractionalHalfAmplitudeGamma);
fringeListFractionalHalfAmplitudeGammaNorm=cell(size(fringeListFractionalHalfAmplitudeGamma));
for k=1:moveNumPart
    fringeListFractionalHalfAmplitudeGammaNorm{k}=(fringeListFractionalHalfAmplitudeGamma{k}-ABFractionalAmplitudeModulatedGamma{1})./ABFractionalAmplitudeModulatedGamma{2};
end
% fringeListFractionalHalfAmplitudeGammaNorm=NormalizeFringe(fringeListFractionalHalfAmplitudeGamma,moveNumPart,ABFractionalAmplitudeModulatedGamma);
wrappedPhaseFractionalHalfAmplitudeGammaNorm=GetWrapPhase(fringeListFractionalHalfAmplitudeGammaNorm,moveNumPart);
% 数步信号Hilbert[Gamma]
fringeListFractionalHilbertHalfAmplitudeGammaNorm=HilbertPerRow(fringeListFractionalHalfAmplitudeGammaNorm,moveNumPart);
wrappedPhaseFractionalHilbertHalfAmplitudeGammaNorm=GetWrapPhaseWithHilbert(fringeListFractionalHilbertHalfAmplitudeGammaNorm,moveNumPart);

% 显示归一化[Normalization]的信号
figure('name',sprintf('1/%d Step of Fringe (Half Amplitude) [Gamma/Normalization: (S-A)/B]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalHalfAmplitudeGammaNorm{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringe (Half Amplitude) [Gamma/Normalization: (S-A)/B]',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma with Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma/Normalization: (S-A)/B] (Half Amplitude)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpaceNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHalfAmplitudeGammaNorm-wrappedPhaseAllHalfAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHTNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertHalfAmplitudeGammaNorm-wrappedPhaseAllHalfAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 平均相位误差
wrappedErrorMeanNorm=(wrappedErrorSpaceNorm+wrappedErrorHTNorm)/2;
plot(wrappedErrorMeanNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Half Amplitude) [Gamma ($\\gamma$=%1.2f)/Normalization: (S-A)/B]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

end

%% @_@{均值为四分幅度的半幅度Gamma非线性响应条纹标记,归一化：S1-A}☉☉☉☉☉☉☉☉☉☉☉☉
if halfAmplitudeNormWithSubtractingAFlag==1
%% --以非对称连续弧段函数调制所有信号
fringeListAllHalfAmplitude=cell(moveNumAll,1);
fringeListAllHalfAmplitudeGamma=cell(moveNumAll,1);
% 期望的Hilbert变换
% expectedHilbertOfFringeListFractionalHalfAmplitude=cell(moveNumPart,1);
for k=1:moveNumAll
    fringeListAllHalfAmplitude{k}=floor(fringeListAllUnit{k});
    fringeListAllHalfAmplitudeGamma{k}=floor(filterGamma(floor(fringeListAllHalfAmplitude{k})+1));
end
% for k=1:moveNumPart
%     fringeListFractionalHalfAmplitude{k}=floor(fringeListFractional{k}.*filterHalfAmplitude);
%     sf=-period*(k-1)/moveNumPart;
% %     expectedHilbertOfFringeListFractionalHalfAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterHalfAmplitude;
% end

% 抽取出数步相移条纹图像
fringeListFractionalHalfAmplitude     =SelectNStepFring(fringeListAllHalfAmplitude     ,moveNumPart);
fringeListFractionalHalfAmplitudeGamma=SelectNStepFring(fringeListAllHalfAmplitudeGamma,moveNumPart);

%% --计算理想空域相位[Gamma]
wrappedPhaseAllHalfAmplitudeGamma=GetWrapPhase(fringeListAllHalfAmplitudeGamma,moveNumAll);

%% --计算数步相移条纹的空域相位[Gamma]
wrappedPhaseFractionalHalfAmplitudeGamma=GetWrapPhase(fringeListFractionalHalfAmplitudeGamma,moveNumPart);

%% --计算数步相移条纹[Gamma]的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertHalfAmplitudeGamma=HilbertPerRow(fringeListFractionalHalfAmplitudeGamma,moveNumPart);
wrappedPhaseFractionalHilbertHalfAmplitudeGamma=GetWrapPhaseWithHilbert(fringeListFractionalHilbertHalfAmplitudeGamma,moveNumPart);

%% --显示数步相移条纹信号[Gamma]及其Hilbert变换、相位误差图表
% 组合显示原始信号、调制后的组合函数
figure('name',sprintf('1/%d Step of Original Fringe (Half Amplitude) [Gamma]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalHalfAmplitude{1},     plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
plot(fringeListFractionalHalfAmplitudeGamma{1},plotLineType,'Color','m','LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringes (Half Amplitude) [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
legend('Fringe Signal','Fringe Signal [Gamma]','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示理想折叠相位
figure('name','Wrapped Phase (Half Amplitude) [Gamma/Normalization: S-A]','NumberTitle','off');
plot(wrappedPhaseAllHalfAmplitudeGamma,plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Wrapped Phase (Unit0 Amplitude) [Gamma/Normalization: S-A]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma without Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma]','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHalfAmplitudeGamma-wrappedPhaseAllHalfAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertHalfAmplitudeGamma-wrappedPhaseAllHalfAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
% 平均相位误差
wrappedErrorMean=(wrappedErrorSpace+wrappedErrorHT)/2;
plot(wrappedErrorMean,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Half Amplitude) [Gamma ($\\gamma$=%1.2f)]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% 对信号[Gamma]归一化[Normalization]处理(S-A)及归一化后的理想空域相位、Hilbert域相位、平均相位
% 理想信号[Gamma]
ABAllAmplitudeModulatedGamma=GetABNormal(fringeListAllHalfAmplitudeGamma,wrappedPhaseAllHalfAmplitudeGamma);
% 显示背景光强A
figure('name','Background Illumination A (Unit0 Amplitude) [Gamma/Normalization: S-A]','NumberTitle','off');
plot(ABAllAmplitudeModulatedGamma{1},plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Background Illumination A (Unit0 Amplitude) [Gamma/Normalization: S-A]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
% 显示调制度B
figure('name','Amplitude B (Unit0 Amplitude) [Gamma/Normalization: S-A]','NumberTitle','off');
plot(ABAllAmplitudeModulatedGamma{2},plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Amplitude B (Unit0 Amplitude) [Gamma/Normalization: S-A]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
fringeListAllHalfAmplitudeGammaNorm=cell(size(fringeListAllHalfAmplitudeGamma));
for k=1:moveNumAll
    fringeListAllHalfAmplitudeGammaNorm{k}=fringeListAllHalfAmplitudeGamma{k}-ABAllAmplitudeModulatedGamma{1};
end
% fringeListAllHalfAmplitudeGammaNorm=NormalizeFringe(fringeListAllHalfAmplitudeGamma,moveNumAll,ABAllAmplitudeModulatedGamma);
wrappedPhaseAllHalfAmplitudeGammaNorm=GetWrapPhase(fringeListAllHalfAmplitudeGammaNorm,moveNumAll);
% 数步信号[Gamma]
ABFractionalAmplitudeModulatedGamma=GetABNormal(fringeListFractionalHalfAmplitudeGamma,wrappedPhaseFractionalHalfAmplitudeGamma);
fringeListFractionalHalfAmplitudeGammaNorm=cell(size(fringeListFractionalHalfAmplitudeGamma));
for k=1:moveNumPart
    fringeListFractionalHalfAmplitudeGammaNorm{k}=fringeListFractionalHalfAmplitudeGamma{k}-ABFractionalAmplitudeModulatedGamma{1};
end
% fringeListFractionalHalfAmplitudeGammaNorm=NormalizeFringe(fringeListFractionalHalfAmplitudeGamma,moveNumPart,ABFractionalAmplitudeModulatedGamma);
wrappedPhaseFractionalHalfAmplitudeGammaNorm=GetWrapPhase(fringeListFractionalHalfAmplitudeGammaNorm,moveNumPart);
% 数步信号Hilbert[Gamma]
fringeListFractionalHilbertHalfAmplitudeGammaNorm=HilbertPerRow(fringeListFractionalHalfAmplitudeGammaNorm,moveNumPart);
wrappedPhaseFractionalHilbertHalfAmplitudeGammaNorm=GetWrapPhaseWithHilbert(fringeListFractionalHilbertHalfAmplitudeGammaNorm,moveNumPart);

% 显示归一化[Normalization]的信号
figure('name',sprintf('1/%d Step of Fringe (Half Amplitude) [Gamma/Normalization: S-A]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalHalfAmplitudeGammaNorm{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringe (Half Amplitude) [Gamma/Normalization: S-A]',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma with Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma/Normalization: S-A] (Half Amplitude)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpaceNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHalfAmplitudeGammaNorm-wrappedPhaseAllHalfAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHTNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertHalfAmplitudeGammaNorm-wrappedPhaseAllHalfAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 平均相位误差
wrappedErrorMeanNorm=(wrappedErrorSpaceNorm+wrappedErrorHTNorm)/2;
plot(wrappedErrorMeanNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Half Amplitude) [Gamma ($\\gamma$=%1.2f)/Normalization: S-A]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

end

%% @_@{均值为四分幅度的半幅度Gamma非线性响应条纹标记,归一化：S1/B}☉☉☉☉☉☉☉☉☉☉☉☉
if halfAmplitudeNormWithDividingBFlag==1
%% --以非对称连续弧段函数调制所有信号
fringeListAllHalfAmplitude=cell(moveNumAll,1);
fringeListAllHalfAmplitudeGamma=cell(moveNumAll,1);
% 期望的Hilbert变换
% expectedHilbertOfFringeListFractionalHalfAmplitude=cell(moveNumPart,1);
for k=1:moveNumAll
    fringeListAllHalfAmplitude{k}=floor(fringeListAllUnit{k});
    fringeListAllHalfAmplitudeGamma{k}=floor(filterGamma(floor(fringeListAllHalfAmplitude{k})+1));
end
% for k=1:moveNumPart
%     fringeListFractionalHalfAmplitude{k}=floor(fringeListFractional{k}.*filterHalfAmplitude);
%     sf=-period*(k-1)/moveNumPart;
% %     expectedHilbertOfFringeListFractionalHalfAmplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterHalfAmplitude;
% end

% 抽取出数步相移条纹图像
fringeListFractionalHalfAmplitude     =SelectNStepFring(fringeListAllHalfAmplitude     ,moveNumPart);
fringeListFractionalHalfAmplitudeGamma=SelectNStepFring(fringeListAllHalfAmplitudeGamma,moveNumPart);

%% --计算理想空域相位[Gamma]
wrappedPhaseAllHalfAmplitudeGamma=GetWrapPhase(fringeListAllHalfAmplitudeGamma,moveNumAll);

%% --计算数步相移条纹的空域相位[Gamma]
wrappedPhaseFractionalHalfAmplitudeGamma=GetWrapPhase(fringeListFractionalHalfAmplitudeGamma,moveNumPart);

%% --计算数步相移条纹[Gamma]的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertHalfAmplitudeGamma=HilbertPerRow(fringeListFractionalHalfAmplitudeGamma,moveNumPart);
wrappedPhaseFractionalHilbertHalfAmplitudeGamma=GetWrapPhaseWithHilbert(fringeListFractionalHilbertHalfAmplitudeGamma,moveNumPart);

%% --显示数步相移条纹信号[Gamma]及其Hilbert变换、相位误差图表
% 组合显示原始信号、调制后的组合函数
figure('name',sprintf('1/%d Step of Original Fringe (Half Amplitude) [Gamma]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalHalfAmplitude{1},     plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
plot(fringeListFractionalHalfAmplitudeGamma{1},plotLineType,'Color','m','LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringes (Half Amplitude) [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
legend('Fringe Signal','Fringe Signal [Gamma]','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示理想折叠相位
figure('name','Wrapped Phase (Half Amplitude) [Gamma/Normalization: S-A]','NumberTitle','off');
plot(wrappedPhaseAllHalfAmplitudeGamma,plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Wrapped Phase (Unit0 Amplitude) [Gamma/Normalization: S-A]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma without Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma]','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHalfAmplitudeGamma-wrappedPhaseAllHalfAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertHalfAmplitudeGamma-wrappedPhaseAllHalfAmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
% 平均相位误差
wrappedErrorMean=(wrappedErrorSpace+wrappedErrorHT)/2;
plot(wrappedErrorMean,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Half Amplitude) [Gamma ($\\gamma$=%1.2f)]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% 对信号[Gamma]归一化[Normalization]处理(S/B)及归一化后的理想空域相位、Hilbert域相位、平均相位
% 理想信号[Gamma]
ABAllAmplitudeModulatedGamma=GetABNormal(fringeListAllHalfAmplitudeGamma,wrappedPhaseAllHalfAmplitudeGamma);
% 显示背景光强A
figure('name','Background Illumination A (Unit0 Amplitude) [Gamma/Normalization: S-A]','NumberTitle','off');
plot(ABAllAmplitudeModulatedGamma{1},plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Background Illumination A (Unit0 Amplitude) [Gamma/Normalization: S-A]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
% 显示调制度B
figure('name','Amplitude B (Unit0 Amplitude) [Gamma/Normalization: S/B]','NumberTitle','off');
plot(ABAllAmplitudeModulatedGamma{2},plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Amplitude B (Unit0 Amplitude) [Gamma/Normalization: S/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
fringeListAllHalfAmplitudeGammaNorm=cell(size(fringeListAllHalfAmplitudeGamma));
for k=1:moveNumAll
    fringeListAllHalfAmplitudeGammaNorm{k}=fringeListAllHalfAmplitudeGamma{k}./ABAllAmplitudeModulatedGamma{2};
end
% fringeListAllHalfAmplitudeGammaNorm=NormalizeFringe(fringeListAllHalfAmplitudeGamma,moveNumAll,ABAllAmplitudeModulatedGamma);
wrappedPhaseAllHalfAmplitudeGammaNorm=GetWrapPhase(fringeListAllHalfAmplitudeGammaNorm,moveNumAll);
% 数步信号[Gamma]
ABFractionalAmplitudeModulatedGamma=GetABNormal(fringeListFractionalHalfAmplitudeGamma,wrappedPhaseFractionalHalfAmplitudeGamma);
fringeListFractionalHalfAmplitudeGammaNorm=cell(size(fringeListFractionalHalfAmplitudeGamma));
for k=1:moveNumPart
    fringeListFractionalHalfAmplitudeGammaNorm{k}=fringeListFractionalHalfAmplitudeGamma{k}./ABFractionalAmplitudeModulatedGamma{2};
end
% fringeListFractionalHalfAmplitudeGammaNorm=NormalizeFringe(fringeListFractionalHalfAmplitudeGamma,moveNumPart,ABFractionalAmplitudeModulatedGamma);
wrappedPhaseFractionalHalfAmplitudeGammaNorm=GetWrapPhase(fringeListFractionalHalfAmplitudeGammaNorm,moveNumPart);
% 数步信号Hilbert[Gamma]
fringeListFractionalHilbertHalfAmplitudeGammaNorm=HilbertPerRow(fringeListFractionalHalfAmplitudeGammaNorm,moveNumPart);
wrappedPhaseFractionalHilbertHalfAmplitudeGammaNorm=GetWrapPhaseWithHilbert(fringeListFractionalHilbertHalfAmplitudeGammaNorm,moveNumPart);

% 显示归一化[Normalization]的信号
figure('name',sprintf('1/%d Step of Fringe (Half Amplitude) [Gamma/Normalization: S/B]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalHalfAmplitudeGammaNorm{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringe (Half Amplitude) [Gamma/Normalization: S/B]',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma with Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma/Normalization: S/B] (Half Amplitude)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpaceNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHalfAmplitudeGammaNorm-wrappedPhaseAllHalfAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHTNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertHalfAmplitudeGammaNorm-wrappedPhaseAllHalfAmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 平均相位误差
wrappedErrorMeanNorm=(wrappedErrorSpaceNorm+wrappedErrorHTNorm)/2;
plot(wrappedErrorMeanNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Half Amplitude) [Gamma ($\\gamma$=%1.2f)/Normalization: S/B]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

end


%% @_@{均值为0的Gamma非线性响应条纹条纹标记,归一化：S0/B}☉☉☉☉☉☉☉☉☉☉☉☉
if unit0NormWithDividingBFlag==1
%% --以非对称连续弧段函数调制所有信号
fringeListAllUnit0Amplitude=cell(moveNumAll,1);
fringeListAllUnit0AmplitudeGamma=cell(moveNumAll,1);
% 期望的Hilbert变换
% expectedHilbertOfFringeListFractionalUnit0Amplitude=cell(moveNumPart,1);
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllUnit0Amplitude{k}=floor(255*0.5*(cos(((0:lengthOfSignal-1)-sf)/period*2*pi)+0)/1);
    fringeListAllUnit0AmplitudeGamma{k}=floor(...
        ( (fringeListAllUnit0Amplitude{k}>0)-(fringeListAllUnit0Amplitude{k}<0) ).*...
        filterGamma(abs(fringeListAllUnit0Amplitude{k})+1));
end
% for k=1:moveNumPart
%     fringeListFractionalUnit0Amplitude{k}=floor(fringeListFractional{k}.*filterUnit0Amplitude);
%     sf=-period*(k-1)/moveNumPart;
% %     expectedHilbertOfFringeListFractionalUnit0Amplitude{k}=255.0/2*(sin(((0:lengthOfSignal-1)-sf)/period*2*pi))/2.*filterUnit0Amplitude;
% end

% 抽取出数步相移条纹图像
fringeListFractionalUnit0Amplitude     =SelectNStepFring(fringeListAllUnit0Amplitude     ,moveNumPart);
fringeListFractionalUnit0AmplitudeGamma=SelectNStepFring(fringeListAllUnit0AmplitudeGamma,moveNumPart);

%% --计算理想空域相位[Gamma]
wrappedPhaseAllUnit0AmplitudeGamma=GetWrapPhase(fringeListAllUnit0AmplitudeGamma,moveNumAll);

%% --计算数步相移条纹的空域相位[Gamma]
wrappedPhaseFractionalUnit0AmplitudeGamma=GetWrapPhase(fringeListFractionalUnit0AmplitudeGamma,moveNumPart);

%% --计算数步相移条纹[Gamma]的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertUnit0AmplitudeGamma=HilbertPerRow(fringeListFractionalUnit0AmplitudeGamma,moveNumPart);
wrappedPhaseFractionalHilbertUnit0AmplitudeGamma=GetWrapPhaseWithHilbert(fringeListFractionalHilbertUnit0AmplitudeGamma,moveNumPart);

%% --显示数步相移条纹信号[Gamma]及其Hilbert变换、相位误差图表
% 组合显示原始信号、调制后的组合函数
figure('name',sprintf('1/%d Step of Original Fringe (Unit0 Amplitude) [Gamma]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalUnit0Amplitude{1},     plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
plot(fringeListFractionalUnit0AmplitudeGamma{1},plotLineType,'Color','m','LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringes (Unit0 Amplitude) [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
legend('Fringe Signal','Fringe Signal [Gamma]','Location','NorthWest');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-32,160]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示理想折叠相位
figure('name','Wrapped Phase (Unit0 Amplitude) [Gamma/Normalization: S/B]','NumberTitle','off');
plot(wrappedPhaseAllUnit0AmplitudeGamma,plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Wrapped Phase (Unit0 Amplitude) [Gamma/Normalization: S/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma without Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma]','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalUnit0AmplitudeGamma-wrappedPhaseAllUnit0AmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertUnit0AmplitudeGamma-wrappedPhaseAllUnit0AmplitudeGamma,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
% 平均相位误差
wrappedErrorMean=(wrappedErrorSpace+wrappedErrorHT)/2;
plot(wrappedErrorMean,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Unit0 Amplitude) [Gamma ($\\gamma$=%1.2f)]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% --对信号[Gamma]归一化[Normalization]处理(S/B)及归一化后的理想空域相位、Hilbert域相位、平均相位
% 理想信号[Gamma]
ABAllAmplitudeModulatedGamma=GetABNormal(fringeListAllUnit0AmplitudeGamma,wrappedPhaseAllUnit0AmplitudeGamma);
% 显示背景光强A
figure('name','Background Illumination A (Unit0 Amplitude) [Gamma/Normalization: S/B]','NumberTitle','off');
plot(ABAllAmplitudeModulatedGamma{1},plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Background Illumination A (Unit0 Amplitude) [Gamma/Normalization: S/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
% 显示调制度B
figure('name','Amplitude B (Unit0 Amplitude) [Gamma/Normalization: S/B]','NumberTitle','off');
plot(ABAllAmplitudeModulatedGamma{2},plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Amplitude B (Unit0 Amplitude) [Gamma/Normalization: S/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

fringeListAllUnit0AmplitudeGammaNorm=cell(size(fringeListAllUnit0AmplitudeGamma));
for k=1:moveNumAll
    fringeListAllUnit0AmplitudeGammaNorm{k}=fringeListAllUnit0AmplitudeGamma{k}./ABAllAmplitudeModulatedGamma{2};
end
% fringeListAllUnit0AmplitudeGammaNorm=NormalizeFringe(fringeListAllUnit0AmplitudeGamma,moveNumAll,ABAllAmplitudeModulatedGamma);
wrappedPhaseAllUnit0AmplitudeGammaNorm=GetWrapPhase(fringeListAllUnit0AmplitudeGammaNorm,moveNumAll);
% 数步信号[Gamma]
ABFractionalAmplitudeModulatedGamma=GetABNormal(fringeListFractionalUnit0AmplitudeGamma,wrappedPhaseFractionalUnit0AmplitudeGamma);
fringeListFractionalUnit0AmplitudeGammaNorm=cell(size(fringeListFractionalUnit0AmplitudeGamma));
for k=1:moveNumPart
    fringeListFractionalUnit0AmplitudeGammaNorm{k}=fringeListFractionalUnit0AmplitudeGamma{k}./ABFractionalAmplitudeModulatedGamma{2};
end
% fringeListFractionalUnit0AmplitudeGammaNorm=NormalizeFringe(fringeListFractionalUnit0AmplitudeGamma,moveNumPart,ABFractionalAmplitudeModulatedGamma);
wrappedPhaseFractionalUnit0AmplitudeGammaNorm=GetWrapPhase(fringeListFractionalUnit0AmplitudeGammaNorm,moveNumPart);
% 数步信号Hilbert[Gamma]
fringeListFractionalHilbertUnit0AmplitudeGammaNorm=HilbertPerRow(fringeListFractionalUnit0AmplitudeGammaNorm,moveNumPart);
wrappedPhaseFractionalHilbertUnit0AmplitudeGammaNorm=GetWrapPhaseWithHilbert(fringeListFractionalHilbertUnit0AmplitudeGammaNorm,moveNumPart);

% 显示归一化[Normalization]的信号
figure('name',sprintf('1/%d Step of Fringe (Unit0 Amplitude) [Gamma/Normalization: S/B]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalUnit0AmplitudeGammaNorm{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringe (Unit0 Amplitude) [Gamma/Normalization: S/B]',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma with Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma/Normalization: S/B] (Unit0 Amplitude)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpaceNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalUnit0AmplitudeGammaNorm-wrappedPhaseAllUnit0AmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHTNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertUnit0AmplitudeGammaNorm-wrappedPhaseAllUnit0AmplitudeGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 平均相位误差
wrappedErrorMeanNorm=(wrappedErrorSpaceNorm+wrappedErrorHTNorm)/2;
plot(wrappedErrorMeanNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Unit0 Amplitude) [Gamma ($\\gamma$=%1.2f)/Normalization: S/B]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

end

%% @_@{对非对称连续调制频率信号[完整周期],归一化: S3/B}********************
% 对非对称连续调制频率信号[完整周期]标记,归一化: S3/B
if asymmetricalArcFrequencyModulateNormWithDividingBFlag==1
disp('asymmetricalArcFrequencyModulateDividingB');
asymmetricalArcFactor=1.25;
lengthEnd=lengthOfSignal;
    
%% --生成24幅非对称连续调制频率全部条纹图像并从中抽取出数步相移条纹图像
moveNumAll=24;
fringeListAllAsymmetricalArcFrequency=cell(24,1);
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllAsymmetricalArcFrequency{k}=(255.0/1*(cos(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal)))/period)+1)/2);
    fringeListAllAsymmetricalArcFrequency{k}=fringeListAllAsymmetricalArcFrequency{k}(1:lengthEnd);
end

%% --计算理想空域相位
wrappedPhaseAllAsymmetricalArcFrequency=GetWrapPhase(fringeListAllAsymmetricalArcFrequency,moveNumAll);

%% --利用理想折叠相位，逐周期选择信号
periodNumber=0;
eachPeriodFringeStartIndex=zeros(1,1);
for p=1:lengthEnd-1
    if wrappedPhaseAllAsymmetricalArcFrequency(p+1)-wrappedPhaseAllAsymmetricalArcFrequency(p)<-pi/2
        periodNumber=periodNumber+1;
        eachPeriodFringeStartIndex(periodNumber)=p+1;
    end
end

% 截取信号
fringeListAllAsymmetricalArcFrequencyGamma=cell(size(fringeListAllAsymmetricalArcFrequency));
for k=1:moveNumAll
    fringeListAllAsymmetricalArcFrequency{k}=fringeListAllAsymmetricalArcFrequency{k}(eachPeriodFringeStartIndex(1):eachPeriodFringeStartIndex(end)-1);
    % Gamma非线性响应
    fringeListAllAsymmetricalArcFrequencyGamma{k}=floor(filterGamma(floor(fringeListAllAsymmetricalArcFrequency{k})+1));
end

%% --重新计算理想空域相位
wrappedPhaseAllAsymmetricalArcFrequency     =GetWrapPhase(fringeListAllAsymmetricalArcFrequency,     moveNumAll);
wrappedPhaseAllAsymmetricalArcFrequencyGamma=GetWrapPhase(fringeListAllAsymmetricalArcFrequencyGamma,moveNumAll);

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalAsymmetricalArcFrequency     =SelectNStepFring(fringeListAllAsymmetricalArcFrequency,     moveNumPart);
fringeListFractionalAsymmetricalArcFrequencyGamma=SelectNStepFring(fringeListAllAsymmetricalArcFrequencyGamma,moveNumPart);
% % 期望的Hilbert变换
% expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency=cell(size(fringeListFractionalAsymmetricalArcFrequency));
% for k=1:moveNumPart
%     sf=-period*(k-1)/moveNumPart;
%     expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=255.0/2*(sin(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal)))/period))/2;
%     expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}(eachPeriodFringeStartIndex(1):eachPeriodFringeStartIndex(end)-1);
% end

% 更改参数值
lengthOfSignal=eachPeriodFringeStartIndex(end)-eachPeriodFringeStartIndex(1);
lengthEnd=lengthOfSignal;
eachPeriodFringeStartIndex=eachPeriodFringeStartIndex-eachPeriodFringeStartIndex(1)+1;
eachPeriodFringeStartIndex(end)=eachPeriodFringeStartIndex(end)-1;

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalAsymmetricalArcFrequency     =GetWrapPhase(fringeListFractionalAsymmetricalArcFrequency,     moveNumPart);
wrappedPhaseFractionalAsymmetricalArcFrequencyGamma=GetWrapPhase(fringeListFractionalAsymmetricalArcFrequencyGamma,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAsymmetricalArcFrequency     =HilbertPerRow(fringeListFractionalAsymmetricalArcFrequency,     moveNumPart);
fringeListFractionalHilbertAsymmetricalArcFrequencyGamma=HilbertPerRow(fringeListFractionalAsymmetricalArcFrequencyGamma,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcFrequency     =GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcFrequency,     moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcFrequencyGamma=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcFrequencyGamma,moveNumPart);

%% --显示非对称连续函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 组合显示原始信号、调制后的组合函数
figure('name',sprintf('1/%d Step of Original Fringe & Frequency-Modulated Signal by Asymmetrical Arc Function [Gamma]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalAsymmetricalArcFrequency{1},     plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
plot(fringeListFractionalAsymmetricalArcFrequencyGamma{1},plotLineType,'Color','m','LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringes [Gamma ($\\gamma$=%1.2f)]',moveNumPart,gamma),'Interpreter','latex');
legend('Modulated Signal','Modulated Signal [Gamma]','Location','NorthWest');
xlim([0,lengthOfSignal-1]);grid on;%ylim([-32,160]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% [Gamma without Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma] (Frequency Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcFrequencyGamma    -wrappedPhaseAllAsymmetricalArcFrequencyGamma,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertAsymmetricalArcFrequencyGamma-wrappedPhaseAllAsymmetricalArcFrequencyGamma,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
% 平均相位误差
wrappedErrorMean=(wrappedErrorSpace+wrappedErrorHT)/2;
plot(wrappedErrorMean,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Asymmetrical Arc Frequency) [Gamma ($\\gamma$=%1.2f)]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% 对信号[Gamma]归一化[Normalization]处理S/B及归一化后的理想空域相位、Hilbert域相位、平均相位
% 理想信号[Gamma]
ABAllFrequencyModulatedGamma=GetABNormal(fringeListAllAsymmetricalArcFrequencyGamma,wrappedPhaseAllAsymmetricalArcFrequencyGamma);
% 显示背景光强A
figure('name','Background Illumination A (Asymmetrical Arc Frequency) [Gamma/Normalization: S/B]','NumberTitle','off');
plot(ABAllFrequencyModulatedGamma{1},plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Background Illumination A (Asymmetrical Arc Frequency) [Gamma/Normalization: S/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
% 显示调制度B
figure('name','Frequency B (Asymmetrical Arc Frequency) [Gamma/Normalization: S/B]','NumberTitle','off');
plot(ABAllFrequencyModulatedGamma{2},plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Frequency B (Asymmetrical Arc Frequency) [Gamma/Normalization: S/B]');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
fringeListAllAsymmetricalArcFrequencyGammaNorm=cell(size(fringeListAllAsymmetricalArcFrequencyGamma));
for k=1:moveNumAll
    fringeListAllAsymmetricalArcFrequencyGammaNorm{k}=fringeListAllAsymmetricalArcFrequencyGamma{k}./ABAllFrequencyModulatedGamma{2};
end
% fringeListAllAsymmetricalArcFrequencyGammaNorm=NormalizeFringe(fringeListAllAsymmetricalArcFrequencyGamma,moveNumAll,ABAllFrequencyModulatedGamma);
wrappedPhaseAllAsymmetricalArcFrequencyGammaNorm=GetWrapPhase(fringeListAllAsymmetricalArcFrequencyGammaNorm,moveNumAll);
% 数步信号[Gamma]
ABFractionalFrequencyModulatedGamma=GetABNormal(fringeListFractionalAsymmetricalArcFrequencyGamma,wrappedPhaseFractionalAsymmetricalArcFrequencyGamma);
fringeListFractionalAsymmetricalArcFrequencyGammaNorm=cell(size(fringeListFractionalAsymmetricalArcFrequencyGamma));
for k=1:moveNumPart
    fringeListFractionalAsymmetricalArcFrequencyGammaNorm{k}=fringeListFractionalAsymmetricalArcFrequencyGamma{k}./ABFractionalFrequencyModulatedGamma{2};
end
% fringeListFractionalAsymmetricalArcFrequencyGammaNorm=NormalizeFringe(fringeListFractionalAsymmetricalArcFrequencyGamma,moveNumPart,ABFractionalFrequencyModulatedGamma);
wrappedPhaseFractionalAsymmetricalArcFrequencyGammaNorm=GetWrapPhase(fringeListFractionalAsymmetricalArcFrequencyGammaNorm,moveNumPart);
% 数步信号Hilbert[Gamma]
fringeListFractionalHilbertAsymmetricalArcFrequencyGammaNorm=HilbertPerRow(fringeListFractionalAsymmetricalArcFrequencyGammaNorm,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcFrequencyGammaNorm=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcFrequencyGammaNorm,moveNumPart);

% 显示归一化[Normalization]的信号
figure('name',sprintf('1/%d Step of Fringe (Asymmetrical Arc Frequency) [Gamma/Normalization: S/B]',moveNumPart),'NumberTitle','off');
plot(fringeListFractionalAsymmetricalArcFrequencyGammaNorm{1},plotLineType,'LineWidth',1.5,'MarkerSize',2);
title(sprintf('1/%d Step of Fringe (Asymmetrical Arc Frequency) [Gamma/Normalization: S/B]',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% [Gamma with Normalization)显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error [Gamma/Normalization: S/B] (Asymmetrical Arc Frequency)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpaceNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcFrequencyGammaNorm-wrappedPhaseAllAsymmetricalArcFrequencyGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);  
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHTNorm=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertAsymmetricalArcFrequencyGammaNorm-wrappedPhaseAllAsymmetricalArcFrequencyGammaNorm,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 平均相位误差
wrappedErrorMeanNorm=(wrappedErrorSpaceNorm+wrappedErrorHTNorm)/2;
plot(wrappedErrorMeanNorm,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title(sprintf('Phase Error (Asymmetrical Arc Frequency) [Gamma ($\\gamma$=%1.2f)/Normalization: S/B]',gamma),'Interpreter','latex');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
legend('Space Phase Error','HT Phase Error','Mean Phase Error','Location','SouthOutside','Orientation','Horizontal');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

end


