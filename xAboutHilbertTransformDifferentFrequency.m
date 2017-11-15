% xAboutHilbertTransformDifferentFrequency.m
% 徐文宇，20171106，周一
% 改变条纹信号的频率，分析整数周期采样与非整数周期采样对Hilbert变换以及Hilbert域相位的影响
% ver：---
close all;clear;

%% {设置基本参数}*********************************************************
% 图像条纹参数
width=1024; height=800; period=64;

%% -频率调制
% 阶跃调制频率标记
stepFrequencyModulateFlag=1;
% 对称连续调制频率标记
symmetricalArcFrequencyModulateFlag=1;
% 非对称连续调制频率标记
asymmetricalArcFrequencyModulateFlag=1;

% 信号范围
numOfPeriods=16;
startOfSignal=1;endOfSignal=width;
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

%% {阶跃式调制频率}******************************************************
% 阶跃调制频率标记
if stepFrequencyModulateFlag==1
%% -生成24幅阶跃式调制频率全部条纹图像并从中抽取出数步相移条纹图像
moveNumAll=24;
fringeListAllStepFrequency=cell(24,1);
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    for i=1:lengthOfSignal
        if i<lengthOfSignal/2
            fringeListAllStepFrequency{k}(:,i)=floor(255.0/2*(cos((i-1-sf)/period*2*pi)+1)/2);
        else
            T1=period;T2=period*2;
            zz=(T1-T2)*(lengthOfSignal/2-1-sf)/T1;
            fringeListAllStepFrequency{k}(:,i)=floor(255.0/2*(cos((i-1-sf-zz)/T2 *2*pi)+1)/2);
        end
    end
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalStepFrequency=SelectNStepFring(fringeListAllStepFrequency,moveNumPart);
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalStepFrequency=cell(size(fringeListFractionalStepFrequency));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    for i=1:lengthOfSignal
        if i<lengthOfSignal/2
            expectedHilbertOfFringeListFractionalStepFrequency{k}(:,i)=255.0/2*(sin((i-1-sf)/period*2*pi))/2;
        else
            T1=period;T2=period*2;
            zz=(T1-T2)*(lengthOfSignal/2-1-sf)/T1;
            expectedHilbertOfFringeListFractionalStepFrequency{k}(:,i)=255.0/2*(sin((i-1-sf-zz)/T2*2*pi))/2;
        end
    end
end

%% -计算理想空域相位
wrappedPhaseAllStepFrequency=GetWrapPhase(fringeListAllStepFrequency,moveNumAll);

%% -计算数步相移条纹的空域相位
wrappedPhaseFractionalStepFrequency=GetWrapPhase(fringeListFractionalStepFrequency,moveNumPart);

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertStepFrequency=HilbertPerRow(fringeListFractionalStepFrequency,moveNumPart);
wrappedPhaseFractionalHilbertStepFrequency=GetWrapPhaseWithHilbert(fringeListFractionalHilbertStepFrequency,moveNumPart);

%% -显示阶跃函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示信号及其Hilbert变换
figure('name','Frequency Modulated by Step Function','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(fringeListFractionalStepFrequency{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(expectedHilbertOfFringeListFractionalStepFrequency{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 实际Hilbert变换结果
    plot(imag(hilbert(fringeListFractionalStepFrequency{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    if k==moveNumPart
        legend('Original Fringe','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);ylim([-64-8 128+8]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示Hilbert变换期望值与实际值的条纹误差曲线
figure('name','Fringe Error between Expected HT and Actual HT','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    plot(imag(hilbert(fringeListFractionalStepFrequency{k}))-expectedHilbertOfFringeListFractionalStepFrequency{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Frequency Modulated by Step Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalStepFrequency    -wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertStepFrequency-wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Frequency Modulated by Step Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/阶跃式Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------stepFrequencyModulate-------------\n');
fprintf('          Mean of Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('  Max positive of Space Phase Error: %+f\n',max(wrappedErrorSpace));
fprintf('  Max negative of Space Phase Error: %+f\n',min(wrappedErrorSpace));
fprintf('          RMSE of Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max(wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min(wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);
end


%% {对称连续调制频率}******************************************************
% 对称连续调制频率标记
if symmetricalArcFrequencyModulateFlag==1
%% -生成24幅对称连续调制频率全部条纹图像并从中抽取出数步相移条纹图像
moveNumAll=24;
fringeListAllSymmetricalArcFrequency=cell(24,1);
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllSymmetricalArcFrequency{k}=floor(255.0/2*(cos(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/lengthOfSignal))/period)+1)/2);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalSymmetricalArcFrequency=SelectNStepFring(fringeListAllSymmetricalArcFrequency,moveNumPart);
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalSymmetricalArcFrequency=cell(size(fringeListFractionalStepFrequency));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalSymmetricalArcFrequency{k}=255.0/2*(sin(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/lengthOfSignal))/period))/2;
end

%% -计算理想空域相位
wrappedPhaseAllSymmetricalArcFrequency=GetWrapPhase(fringeListAllSymmetricalArcFrequency,moveNumAll);

%% -计算数步相移条纹的空域相位
wrappedPhaseFractionalSymmetricalArcFrequency=GetWrapPhase(fringeListFractionalSymmetricalArcFrequency,moveNumPart);

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertSymmetricalArcFrequency=HilbertPerRow(fringeListFractionalSymmetricalArcFrequency,moveNumPart);
wrappedPhaseFractionalHilbertSymmetricalArcFrequency=GetWrapPhaseWithHilbert(fringeListFractionalHilbertSymmetricalArcFrequency,moveNumPart);

%% -显示对称连续函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示对称连续
symmetricalArcFrequencyFunction=(0:lengthOfSignal-1)+100*sin(2*pi*(0:lengthOfSignal-1)/width);
figure('name','Symmetrical Arc Function','NumberTitle','off');
plot(symmetricalArcFrequencyFunction,plotLineType,'LineWidth',1,'MarkerSize',2);hold on;
plot((0:lengthOfSignal-1),'g-.','LineWidth',0.5);
title('Symmetrical Arc Function');
xlim([0,lengthOfSignal-1]);ylim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', xTick);set(gca, 'YTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
figure('name','Frequency Modulated by Symmetrical Arc Function','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(fringeListFractionalSymmetricalArcFrequency{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(expectedHilbertOfFringeListFractionalSymmetricalArcFrequency{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 实际Hilbert变换结果
    plot(imag(hilbert(fringeListFractionalSymmetricalArcFrequency{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    if k==moveNumPart
        legend('Original Fringe','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);ylim([-64-8 128+8]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示Hilbert变换期望值与实际值的条纹误差曲线
figure('name','Fringe Error between Expected HT and Actual HT','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 二者结果理论上本应一致，即差值为0，
    plot(zeros(1,lengthOfSignal),plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    % 但受公式计算取整的影响，下式相减结果会有不到1个像素值的误差。
    % plot(imag(hilbert(fringeListFractionalSymmetricalArcFrequency{k}))-expectedHilbertOfFringeListFractionalSymmetricalArcFrequency{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Frequency Modulated by Symmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalSymmetricalArcFrequency    -wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertSymmetricalArcFrequency-wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Frequency Modulated by Symmetrical Arc Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/对称连续式Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------symmetricalArcFrequencyModulate-------------\n');
fprintf('          Mean of Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('  Max positive of Space Phase Error: %+f\n',max(wrappedErrorSpace));
fprintf('  Max negative of Space Phase Error: %+f\n',min(wrappedErrorSpace));
fprintf('          RMSE of Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max(wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min(wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);

end


%% {非对称连续调制频率}******************************************************
% 非对称连续调制频率标记
if asymmetricalArcFrequencyModulateFlag==1
%% -生成24幅非对称连续调制频率全部条纹图像并从中抽取出数步相移条纹图像
moveNumAll=24;
fringeListAllAsymmetricalArcFrequency=cell(24,1);
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllAsymmetricalArcFrequency{k}=floor(255.0/2*(cos(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(1.25*lengthOfSignal)))/period)+1)/2);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalAsymmetricalArcFrequency=SelectNStepFring(fringeListAllAsymmetricalArcFrequency,moveNumPart);
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency=cell(size(fringeListFractionalStepFrequency));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=255.0/2*(sin(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(1.25*lengthOfSignal)))/period))/2;
end

%% -计算理想空域相位
wrappedPhaseAllAsymmetricalArcFrequency=GetWrapPhase(fringeListAllAsymmetricalArcFrequency,moveNumAll);

%% -计算数步相移条纹的空域相位
wrappedPhaseFractionalAsymmetricalArcFrequency=GetWrapPhase(fringeListFractionalAsymmetricalArcFrequency,moveNumPart);

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAsymmetricalArcFrequency=HilbertPerRow(fringeListFractionalAsymmetricalArcFrequency,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcFrequency=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcFrequency,moveNumPart);

%% -显示非对称连续函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示非对称连续
symmetricalArcFrequencyFunction=(0:lengthOfSignal-1)+100*sin(2*pi*(0:lengthOfSignal-1)/(1.25*width));
figure('name','Asymmetrical Arc Function','NumberTitle','off');
plot(symmetricalArcFrequencyFunction,plotLineType,'LineWidth',1,'MarkerSize',2);hold on;
plot((0:lengthOfSignal-1),'g-.','LineWidth',0.5);
title('Asymmetrical Arc Function');
xlim([0,lengthOfSignal-1]);ylim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', xTick);set(gca, 'YTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
figure('name','Frequency Modulated by Asymmetrical Arc Function','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(fringeListFractionalAsymmetricalArcFrequency{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 实际Hilbert变换结果
    plot(imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    if k==moveNumPart
        legend('Original Fringe','Expected HT','Actual HT','Location','SouthOutside','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);ylim([-128-8 128+8]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示Hilbert变换期望值与实际值的条纹误差曲线
figure('name','Fringe Error between Expected HT and Actual HT','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);
    plot(imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{k}))-expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Frequency Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcFrequency    -wrappedPhaseAllAsymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% Hilbert域相位误差
wrappedErrorHT=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertAsymmetricalArcFrequency-wrappedPhaseAllAsymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHT,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Phase Error (Frequency Modulated by Asymmetrical Arc Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/非对称连续式Hilbert域相位误差的平均值、峰值与均方根
fprintf('------------asymmetricalArcFrequencyModulate-------------\n');
fprintf('          Mean of Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('  Max positive of Space Phase Error: %+f\n',max(wrappedErrorSpace));
fprintf('  Max negative of Space Phase Error: %+f\n',min(wrappedErrorSpace));
fprintf('          RMSE of Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max(wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min(wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);

end





