% xAboutHilbertTransformDifferentFrequency.m
% 徐文宇，20171106，周一
% 改变条纹信号的频率，分析整数周期采样与非整数周期采样对Hilbert变换以及Hilbert域相位的影响
% ver：---
close all;clear;

%% @_@{设置基本参数}*********************************************************
% 图像条纹参数
width=1024; height=800; period=128;

%% --阶跃调制频率标记
stepFrequencyModulateFlag=0;
%% --对称连续调制频率标记
symmetricalArcFrequencyModulateFlag=0;
%% --非对称连续调制频率标记
asymmetricalArcFrequencyModulateFlag=0;
%% --对非对称连续调制频率信号[非完整周期]进行逐相位周期Hilbert变换分析
asymmetricalArcFrequencyModulateEachPeriodFlag=0;
%% --非对称连续调制频率延拓标记
asymmetricalArcFrequencyModulateContinuationFlag=0;
%% --对非对称连续调制频率信号[非完整周期]进行逐相位周期Hilbert变换分析
asymmetricalArcFrequencyModulateEachPeriodFlag=1;

% 信号范围
numOfPeriods=8;
startOfSignal=1;endOfSignal=numOfPeriods*period;
lengthOfSignal=endOfSignal-startOfSignal+1;

% 横坐标序列
horizontalIndex=0:lengthOfSignal-1;

% xTick & xTickLabel
if mod(lengthOfSignal,period)==0
    xTick=zeros(1,numOfPeriods+1);
    xTickLabel=cell(1,numOfPeriods+1);
    for k=0:numOfPeriods
        xTick(k+1)=floor(k*period); xTickLabel{k+1}=num2str(xTick(k+1));
    end
    xTick(end)=lengthOfSignal-1; xTickLabel{end}=num2str(lengthOfSignal-1);
else
    xTick=zeros(1,numOfPeriods+2);
    xTickLabel=cell(1,numOfPeriods+2);
    for k=0:numOfPeriods
        xTick(k+1)=floor(k*period); xTickLabel{k+1}=num2str(xTick(k+1));
    end
    xTick(end)=lengthOfSignal-1; xTickLabel{end}=num2str(lengthOfSignal-1);
end
% yTick & yTickLabel
yTickNum=8;
yTick=zeros(1,yTickNum+1);
yTickLabel=cell(1,yTickNum+1);
yTick(yTickNum/2+1)=0;
yTickLabel{yTickNum/2+1}='0';
for k=1:yTickNum/2
    yTick(yTickNum/2+1+k)=floor( k*256/(yTickNum/2)); yTickLabel{yTickNum/2+1+k}=num2str(yTick(yTickNum/2+1+k)); 
    yTick(yTickNum/2+1-k)=floor(-k*256/(yTickNum/2)); yTickLabel{yTickNum/2+1-k}=num2str(yTick(yTickNum/2+1-k));
end
% xTickPart & xTickLabelPart for Fourier Spectrum
partNum=30;
xTickFourierSpectrum=zeros(1,partNum+1);
xTickLabelFourierSpectrum=cell(1,partNum+1);
for k=0:partNum/2
    xTickFourierSpectrum(partNum/2+k+1)=lengthOfSignal/2+1+k+k; xTickLabelFourierSpectrum{partNum/2+k+1}=num2str( k+k);
    xTickFourierSpectrum(partNum/2-k+1)=lengthOfSignal/2+1-k-k; xTickLabelFourierSpectrum{partNum/2-k+1}=num2str(-k-k);
end
xTickFourierSpectrum(partNum/2+1)=lengthOfSignal/2+1; xTickLabelFourierSpectrum{partNum/2+1}=num2str(0);

% 相位误差显示有效区间
upPhaseErrorBound=7; bottomPhaseErrorBound=-7;

% plot画线类型
plotLineType='';        % '' 实线
plotDottedLineType=':'; % ':'虚线

%% @_@{阶跃式调制频率}******************************************************
% 阶跃调制频率标记
if stepFrequencyModulateFlag==1
%% --生成24幅阶跃式调制频率全部条纹图像并从中抽取出数步相移条纹图像
moveNumAll=24;
fringeListAllStepFrequency=cell(24,1);
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    for i=1:lengthOfSignal
        if i<lengthOfSignal/2
%             fringeListAllStepFrequency{k}(:,i)=floor(255.0/2*(cos((i-1-sf)/period*2*pi)+1)/2);
            fringeListAllStepFrequency{k}(:,i)=floor(255.0/2*(cos((i-1-sf)/period*2*pi)+0)/2);
        else
            T1=period;T2=period*2;
            zz=(T1-T2)*(lengthOfSignal/2-1-sf)/T1;
%             fringeListAllStepFrequency{k}(:,i)=floor(255.0/2*(cos((i-1-sf-zz)/T2 *2*pi)+1)/2);
            fringeListAllStepFrequency{k}(:,i)=floor(255.0/2*(cos((i-1-sf-zz)/T2 *2*pi)+0)/2);
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

%% --计算理想空域相位
wrappedPhaseAllStepFrequency=GetWrapPhase(fringeListAllStepFrequency,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalStepFrequency=GetWrapPhase(fringeListFractionalStepFrequency,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertStepFrequency=HilbertPerRow(fringeListFractionalStepFrequency,moveNumPart);
wrappedPhaseFractionalHilbertStepFrequency=GetWrapPhaseWithHilbert(fringeListFractionalHilbertStepFrequency,moveNumPart);

%% --显示阶跃函数、数步相移条纹信号及其Hilbert变换、相位误差图表
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
fprintf('        Mean of   Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('Max positive of   Space Phase Error: %+f\n',max( wrappedErrorSpace));
fprintf('Max negative of   Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('        RMSE of   Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);
end


%% @_@{对称连续调制频率}******************************************************
% 对称连续调制频率标记
if symmetricalArcFrequencyModulateFlag==1
%% --生成24幅对称连续调制频率全部条纹图像并从中抽取出数步相移条纹图像
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
expectedHilbertOfFringeListFractionalSymmetricalArcFrequency=cell(size(fringeListFractionalSymmetricalArcFrequency));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalSymmetricalArcFrequency{k}=255.0/2*(sin(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/lengthOfSignal))/period))/2;
end

%% --计算理想空域相位
wrappedPhaseAllSymmetricalArcFrequency=GetWrapPhase(fringeListAllSymmetricalArcFrequency,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalSymmetricalArcFrequency=GetWrapPhase(fringeListFractionalSymmetricalArcFrequency,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertSymmetricalArcFrequency=HilbertPerRow(fringeListFractionalSymmetricalArcFrequency,moveNumPart);
wrappedPhaseFractionalHilbertSymmetricalArcFrequency=GetWrapPhaseWithHilbert(fringeListFractionalHilbertSymmetricalArcFrequency,moveNumPart);

%% --显示对称连续函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示对称连续
symmetricalArcFrequencyFunction=(0:lengthOfSignal-1)+100*sin(2*pi*(0:lengthOfSignal-1)/lengthOfSignal);
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
fprintf('        Mean of   Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('Max positive of   Space Phase Error: %+f\n',max( wrappedErrorSpace));
fprintf('Max negative of   Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('        RMSE of   Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);

end


%% @_@{非对称连续调制频率,完整周期}******************************************************
% 非对称连续调制频率标记
if asymmetricalArcFrequencyModulateFlag==1
    asymmetricalArcFactor=1.25;
    lengthEnd=lengthOfSignal;
    
%% --生成24幅非对称连续调制频率全部条纹图像并从中抽取出数步相移条纹图像
moveNumAll=24;
fringeListAllAsymmetricalArcFrequency=cell(24,1);
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllAsymmetricalArcFrequency{k}=(255.0/2*(cos(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal)))/period)+0)/2);
    fringeListAllAsymmetricalArcFrequency{k}=fringeListAllAsymmetricalArcFrequency{k}(1:lengthEnd);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalAsymmetricalArcFrequency=SelectNStepFring(fringeListAllAsymmetricalArcFrequency,moveNumPart);
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency=cell(size(fringeListFractionalAsymmetricalArcFrequency));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=255.0/2*(sin(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal)))/period))/2;
    expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}(1:lengthEnd);
end

%% --计算理想空域相位
wrappedPhaseAllAsymmetricalArcFrequency=GetWrapPhase(fringeListAllAsymmetricalArcFrequency,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalAsymmetricalArcFrequency=GetWrapPhase(fringeListFractionalAsymmetricalArcFrequency,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAsymmetricalArcFrequency=HilbertPerRow(fringeListFractionalAsymmetricalArcFrequency,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcFrequency=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcFrequency,moveNumPart);

%% --显示非对称连续函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示非对称连续
symmetricalArcFrequencyFunction=(0:lengthOfSignal-1)+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal));
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
fprintf('        Mean of   Space Phase Error: %+f\n',mean(wrappedErrorSpace));
fprintf('Max positive of   Space Phase Error: %+f\n',max( wrappedErrorSpace));
fprintf('Max negative of   Space Phase Error: %+f\n',min( wrappedErrorSpace));
fprintf('        RMSE of   Space Phase Error: %+f\n',sqrt(sum((wrappedErrorSpace-mean(wrappedErrorSpace)).^2))/lengthOfSignal);
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(wrappedErrorHT));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max( wrappedErrorHT));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min( wrappedErrorHT));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((wrappedErrorHT-mean(wrappedErrorHT)).^2))/lengthOfSignal);

end


%% @_@{测试对非对称调频进行峰值延拓}***********************************
if asymmetricalArcFrequencyModulateContinuationFlag==1
asymmetricalArcFactor=1.5;
%% --生成24幅非对称连续调制频率全部条纹图像并从中抽取出数步相移条纹图像
moveNumAll=24;
fringeListAllAsymmetricalArcFrequency=cell(24,1);
% 初始相位
initialPhase=0;
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    sf=initialPhase-period*(k-1)/moveNumAll;
    fringeListAllAsymmetricalArcFrequency{k}=floor(255.0/2*(cos(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal)))/period)+0)/2);
    fringeListAllAsymmetricalArcFrequency{k}=fringeListAllAsymmetricalArcFrequency{k};
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalAsymmetricalArcFrequency=SelectNStepFring(fringeListAllAsymmetricalArcFrequency,moveNumPart);
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency=cell(size(fringeListFractionalAsymmetricalArcFrequency));
for k=1:moveNumPart
    sf=initialPhase-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=255.0/2*(sin(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal)))/period))/2;
    expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k};
end

%% --计算理想空域相位
wrappedPhaseAllAsymmetricalArcFrequency=GetWrapPhase(fringeListAllAsymmetricalArcFrequency,moveNumAll);

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalAsymmetricalArcFrequency=GetWrapPhase(fringeListFractionalAsymmetricalArcFrequency,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAsymmetricalArcFrequency=HilbertPerRow(fringeListFractionalAsymmetricalArcFrequency,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcFrequency=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcFrequency,moveNumPart);

%% --显示非对称连续函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示非对称连续
symmetricalArcFrequencyFunction=(0:lengthOfSignal-1)+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal));
figure('name','Asymmetrical Arc Function','NumberTitle','off');
plot(horizontalIndex,symmetricalArcFrequencyFunction,plotLineType,'LineWidth',1,'MarkerSize',2);hold on;
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
    plot(horizontalIndex,fringeListFractionalAsymmetricalArcFrequency{k},                       plotLineType,'LineWidth',1.0);hold on;
    % 期望Hilbert变换结果
    plot(horizontalIndex,expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 实际Hilbert变换结果
    plot(horizontalIndex,imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{k})),        plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Original Fringe and its Hilbert Transform',k,moveNumPart));
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',[]);
    if k==moveNumPart
        legend({'Original Fringe','Expected HT','Actual HT'},'Location','NorthEast','Orientation','Horizontal','FontSize',9);
        set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    end
    
    xlim([0,lengthOfSignal-1]);ylim([-64-32 64+32]);grid on;
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

figure('name','Continuation Analysis','NumberTitle','off');
subplot(moveNumPart+1,1,1);hold on;
plot(horizontalIndex,wrappedPhaseFractionalAsymmetricalArcFrequency);
title('Wrapped Phase');
xlim([0,lengthOfSignal]);
ylim([-pi-pi/8,pi+pi/8]);grid on;
set(gca, 'XTick',xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick',-pi:pi/2:pi);
a=gca;
a.TickLabelInterpreter = 'latex';
S = sym(-pi:pi/2:pi); % S = sym(round(vpa(S/pi*2))*pi/2);
a.YTickLabel = strcat('$',arrayfun(@latex, S, 'UniformOutput', false),'$');

%% 搜索折叠相位中边界两端的波峰位置(对应于条纹信号的波谷)，对条纹信号进行整数周期补齐
% 延拓后的折叠相位
wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation=wrappedPhaseFractionalAsymmetricalArcFrequency;
%% 结束端(右端)->>>>->>>>
% 搜索折叠相位中最靠近右边界的临近结束端第一个波峰位置
wrappedPhaseEndPeakIndex=lengthOfSignal;
for p=lengthOfSignal-1:-1:1
    if wrappedPhaseFractionalAsymmetricalArcFrequency(p+1)-wrappedPhaseFractionalAsymmetricalArcFrequency(p)<-pi/2
        wrappedPhaseEndPeakIndex=p;
        break;
    end
end
fprintf('1st Index of Trough at the  End  of Wrapped Phase: %d\n',wrappedPhaseEndPeakIndex);
hold on;plot(wrappedPhaseEndPeakIndex-1,wrappedPhaseFractionalAsymmetricalArcFrequency(p),'ro');
text(wrappedPhaseEndPeakIndex+5,wrappedPhaseFractionalAsymmetricalArcFrequency(p),num2str(wrappedPhaseEndPeakIndex));

% 结束端相位值
wrappedPhaseEnd=wrappedPhaseFractionalAsymmetricalArcFrequency(end);
% 在临近结束端波峰位置左侧与结束端对应的位置索引
continuationEndSinceFromIndex=0;
% 每步相移结束端需要延拓补齐的长度
continuationEndLength=zeros(moveNumPart,1);
% 折叠相位中由结束端映射至前一个周期同等值，二者之间的距离长度(必须>=0)
continuationEndMapDistance=0;
% 从临近结束端波峰的左一个位置处开始向内往左数，第一个小于等于结束端相位值的位置值
for p=wrappedPhaseEndPeakIndex-1:-1:1
    if wrappedPhaseFractionalAsymmetricalArcFrequency(p)<=wrappedPhaseEnd
        continuationEndSinceFromIndex=p+1;
        continuationEndMapDistance=wrappedPhaseEndPeakIndex-1-continuationEndSinceFromIndex+1;
        
        % 若结束端已经是相位波峰，可能会出现相位值比临近结束端波峰左一个位置处的相位值要大，即修正位置与长度
        if continuationEndMapDistance<0
            continuationEndSinceFromIndex=p-continuationEndMapDistance;
            continuationEndMapDistance=0;
        end
        break
    end
end

% 数据太少，未找到位置，直接返回
if continuationEndSinceFromIndex==0
    disp('##Error: continuationEndSinceFromIndex=0');
    return
end

% 标记水平测量线,以便搜索上一个周期折叠相位中的映射位置
plot([lengthOfSignal-1,continuationEndSinceFromIndex-1],[wrappedPhaseEnd,wrappedPhaseEnd],plotDottedLineType,'Color','g','LineWidth',1.5);
% 在折叠相位中，标记首次延拓的复制节段，注意不包含相位的峰值点
plot((continuationEndSinceFromIndex:wrappedPhaseEndPeakIndex-1)-1, ...
    wrappedPhaseFractionalAsymmetricalArcFrequency(continuationEndSinceFromIndex:wrappedPhaseEndPeakIndex-1),...
    plotDottedLineType,'Color','m','LineWidth',2);
% (第一次延拓)从折叠相位结束端向外侧进行首次复制延拓，复制数据为：映射位置至临近结束端波峰波峰(不包含波峰点)
wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation=[wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation,...
    wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(continuationEndSinceFromIndex:wrappedPhaseEndPeakIndex-1)];
plot(lengthOfSignal-1:lengthOfSignal+continuationEndMapDistance-1,...
    wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(lengthOfSignal:lengthOfSignal+continuationEndMapDistance),...
    plotDottedLineType,'Color','m','LineWidth',2);
% (第二次延拓)参考前一次延拓补齐周期进行满幅周期复制延拓,从临近结束端波峰波峰(包含波峰点)至结束端(不包含波峰点)
wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation=[wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation,...
    wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(wrappedPhaseEndPeakIndex:lengthOfSignal+continuationEndMapDistance)];
plot(lengthOfSignal+continuationEndMapDistance-1:lengthOfSignal+continuationEndMapDistance+lengthOfSignal+continuationEndMapDistance-wrappedPhaseEndPeakIndex,...
             wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(lengthOfSignal+continuationEndMapDistance:end),...
             plotDottedLineType,'Color','b','LineWidth',2);

xlim([0,lengthOfSignal+continuationEndMapDistance+lengthOfSignal+continuationEndMapDistance-wrappedPhaseEndPeakIndex]);
fprintf('       Mapped Index at the  End  of Wrapped Phase: %d\n',continuationEndSinceFromIndex);
fprintf('Continuation Length at the  End  of Wrapped Phase: %d\n',continuationEndMapDistance);

% 参照相位补齐小段后，再根据折叠相位分成数段的阈值，根据该阈值位置处对应的索引进行第二次延拓补齐
fringeListContinuationFractionalAsymmetricalArcFrequency=cell(size(fringeListFractionalAsymmetricalArcFrequency));
for k=1:moveNumPart
    fringeListContinuationFractionalAsymmetricalArcFrequency{k}=fringeListFractionalAsymmetricalArcFrequency{k};
end
phaseSegmentThreshold=pi:-2*pi/moveNumPart:-pi+2*pi/moveNumPart;
for k=1:moveNumPart
    subplot(moveNumPart+1,1,k+1);hold on;
    % 原始条纹信号
    plot(horizontalIndex,fringeListContinuationFractionalAsymmetricalArcFrequency{k},plotLineType,'LineWidth',1.0);
    % 标记首次延拓的复制节段：映射位置至临近结束端波峰波峰(不包含波峰点)
    plot((continuationEndSinceFromIndex:wrappedPhaseEndPeakIndex-1)-1,...
        fringeListContinuationFractionalAsymmetricalArcFrequency{k}(continuationEndSinceFromIndex:wrappedPhaseEndPeakIndex-1),...
        plotDottedLineType,'Color','m','LineWidth',2);
    % (第一次延拓)向外侧进行首次复制延拓，复制数据为：相位映射位置至相位临近结束端波峰波峰(不包含波峰点)
    % 此次延拓仅有第一步相移信号能延拓至波谷位置
    fringeListContinuationFractionalAsymmetricalArcFrequency{k}=[fringeListContinuationFractionalAsymmetricalArcFrequency{k},...
        fringeListContinuationFractionalAsymmetricalArcFrequency{k}(continuationEndSinceFromIndex:wrappedPhaseEndPeakIndex-1)];
    plot(lengthOfSignal-1:lengthOfSignal+continuationEndMapDistance-1,...
        fringeListContinuationFractionalAsymmetricalArcFrequency{k}(lengthOfSignal:lengthOfSignal+continuationEndMapDistance),...
        plotDottedLineType,'Color','m','LineWidth',2);
    
    % (第二次延拓)将折叠相位等位分成moveNumPart子段，根据相移步数有区别地进行延拓补齐，延拓后的结束端位于信号波谷处
    % 根据折叠相位分成数段的阈值位置处对应的索引进行第二次延拓补齐
    if k==1 % 第一步相移参考前一次延拓补齐周期进行满幅复制延拓
        fringeListContinuationFractionalAsymmetricalArcFrequency{k}=[fringeListContinuationFractionalAsymmetricalArcFrequency{k},...
            fringeListContinuationFractionalAsymmetricalArcFrequency{k}(wrappedPhaseEndPeakIndex:lengthOfSignal+continuationEndMapDistance)];
        plot(lengthOfSignal+continuationEndMapDistance-1:lengthOfSignal+continuationEndMapDistance+lengthOfSignal+continuationEndMapDistance-wrappedPhaseEndPeakIndex,...
            fringeListContinuationFractionalAsymmetricalArcFrequency{k}(lengthOfSignal+continuationEndMapDistance:end),...
            plotDottedLineType,'Color','b','LineWidth',2);
    else % 后续相移进行部分节段复制延拓
        for p=wrappedPhaseEndPeakIndex+1:lengthOfSignal+continuationEndMapDistance
            if wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(p)>=phaseSegmentThreshold(k)
                fringeListContinuationFractionalAsymmetricalArcFrequency{k}=[fringeListContinuationFractionalAsymmetricalArcFrequency{k},...
                    fringeListContinuationFractionalAsymmetricalArcFrequency{k}(wrappedPhaseEndPeakIndex:p)];
                plot(lengthOfSignal+continuationEndMapDistance-1:lengthOfSignal+continuationEndMapDistance+p-wrappedPhaseEndPeakIndex,...
                    fringeListContinuationFractionalAsymmetricalArcFrequency{k}(lengthOfSignal+continuationEndMapDistance:end),...
                    plotDottedLineType,'Color','b','LineWidth',2);
                break;
            end
        end
    end
    
    % 每步相移结束端需要延拓补齐的长度
    continuationEndLength(k)=length(fringeListContinuationFractionalAsymmetricalArcFrequency{k})-lengthOfSignal;
    
    title(sprintf('%d/%d Step of Original Fringe',k,moveNumPart));
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',[]);
    if k==moveNumPart
        set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    end
    
    xlim([0,length(fringeListContinuationFractionalAsymmetricalArcFrequency{1})]);ylim([-64-32 64+32]);grid on;
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end
% 结束端(右端)-<<<<-<<<<

%% 起始端(左端)-<<<<-<<<<
% 搜索折叠相位中最靠近左边界的临近结束端第一个波峰位置，
wrappedPhaseStartPeakIndex=1;
for p=1:1:lengthOfSignal-1
    if wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(p+1)-wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(p)<-pi/2
        wrappedPhaseStartPeakIndex=p;
        break;
    end
end
fprintf('1st Index of Trough at the Start of Wrapped Phase: %d\n',wrappedPhaseStartPeakIndex);
subplot(moveNumPart+1,1,1);
plot(wrappedPhaseStartPeakIndex-1,wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(p),'ro');
text(wrappedPhaseStartPeakIndex+5,wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(p),num2str(wrappedPhaseStartPeakIndex));

% 起始端相位值
wrappedPhaseStart=wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(1);
% 在临近起始端波峰位置右侧与起始端对应的位置索引
continuationStartSinceFromIndex=0;
% 每步相移结束端需要延拓补齐的长度
continuationStartLength=zeros(moveNumPart,1);
% 折叠相位中由起始端映射至前一个周期同等值，二者之间的距离长度(必须>=0)
continuationStartMapDistance=0;
% 从临近起始端波峰的右一个位置处开始向内往右数，第一个小于等于起始端相位值的位置值
for p=wrappedPhaseStartPeakIndex+1:1:lengthOfSignal-1
    if wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(p)>=wrappedPhaseStart
        continuationStartSinceFromIndex=p-1;
        continuationStartMapDistance=continuationStartSinceFromIndex-1-wrappedPhaseStartPeakIndex+1;
        
        % 若起始端已经是相位波峰，可能会出现相位值比临近起始端波峰右一个位置处的相位值要大，即修正位置与长度
        if continuationStartMapDistance<0
            continuationStartSinceFromIndex=p-continuationStartMapDistance;
            continuationStartMapDistance=0;
        end      
        break
    end
end

% 标记水平测量线,以便搜索上一个周期折叠相位中的映射位置
subplot(moveNumPart+1,1,1);
plot([0,continuationStartSinceFromIndex-1],[wrappedPhaseStart,wrappedPhaseStart],plotDottedLineType,'Color','g','LineWidth',1.5);
% 在折叠相位中，标记首次延拓的复制节段，注意不包含相位的峰值点
plot((wrappedPhaseStartPeakIndex+1:continuationStartSinceFromIndex)-1, ...
    wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(wrappedPhaseStartPeakIndex+1:continuationStartSinceFromIndex),...
    plotDottedLineType,'Color','m','LineWidth',2);
% (第一次延拓)从折叠相位起始端向外侧进行首次复制延拓，复制数据为：临近起始端波峰波峰(不包含波峰点)至映射位置
wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation=[wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(wrappedPhaseStartPeakIndex+1:continuationStartSinceFromIndex),...
    wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation];
plot(-(continuationStartSinceFromIndex-wrappedPhaseStartPeakIndex):1:0,...
    wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(1:continuationStartSinceFromIndex-wrappedPhaseStartPeakIndex+1),...
    plotDottedLineType,'Color','m','LineWidth',2);
% (第二次延拓)参考前一次延拓补齐周期进行满幅周期复制延拓,起始端至从临近起始端波峰(包含波峰点)
wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation=[wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(1:continuationStartSinceFromIndex),...
    wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation];
plot(-continuationStartSinceFromIndex*2+wrappedPhaseStartPeakIndex:1:-(continuationStartSinceFromIndex-wrappedPhaseStartPeakIndex),...
    wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(1:continuationStartSinceFromIndex+1),...
    plotDottedLineType,'Color','b','LineWidth',2);

xlim([-continuationStartSinceFromIndex*2+wrappedPhaseStartPeakIndex,lengthOfSignal+continuationEndMapDistance+lengthOfSignal+continuationEndMapDistance-wrappedPhaseEndPeakIndex]);
fprintf('       Mapped Index at the Start of Wrapped Phase: %d\n',continuationStartSinceFromIndex);
fprintf('Continuation Length at the Start of Wrapped Phase: %d\n',continuationStartMapDistance);

% 折叠相位分成数段的阈值，根据该阈值位置处对应的索引进行第二次延拓补齐
phaseSegmentThreshold=pi:-2*pi/moveNumPart:-pi+2*pi/moveNumPart;
for k=1:moveNumPart
    subplot(moveNumPart+1,1,k+1);hold on;
    % 标记首次延拓的复制节段
    plot((wrappedPhaseStartPeakIndex+1:continuationStartSinceFromIndex)-1, ...
        fringeListContinuationFractionalAsymmetricalArcFrequency{k}(wrappedPhaseStartPeakIndex+1:continuationStartSinceFromIndex),...
        plotDottedLineType,'Color','m','LineWidth',2);
    % (第一次延拓)从信号起始端向外侧进行首次复制延拓，复制数据为：临近起始端波峰波峰(不包含波峰点)至映射位置
    fringeListContinuationFractionalAsymmetricalArcFrequency{k}=[fringeListContinuationFractionalAsymmetricalArcFrequency{k}(wrappedPhaseStartPeakIndex+1:continuationStartSinceFromIndex),...
        fringeListContinuationFractionalAsymmetricalArcFrequency{k}];
    plot(-(continuationStartSinceFromIndex-wrappedPhaseStartPeakIndex):1:0,...
        fringeListContinuationFractionalAsymmetricalArcFrequency{k}(1:continuationStartSinceFromIndex-wrappedPhaseStartPeakIndex+1),...
        plotDottedLineType,'Color','m','LineWidth',2);
    
    % (第二次延拓)将折叠相位等位分成moveNumPart子段，根据相移步数有区别地进行延拓补齐，延拓后的起始端位于信号波谷处
    % 根据折叠相位分成数段的阈值位置处对应的索引进行第二次延拓补齐
    if k==1 % 第一步相移参考前一次延拓补齐周期进行满幅复制延拓
%         fringeListFractionalAsymmetricalArcFrequency{k}=[fringeListFractionalAsymmetricalArcFrequency{k}(1:continuationStartSinceFromIndex),...
%             fringeListFractionalAsymmetricalArcFrequency{k}];
%         plot(-continuationStartSinceFromIndex*2+wrappedPhaseStartPeakIndex:1:-(continuationStartSinceFromIndex-wrappedPhaseStartPeakIndex),...
%             fringeListFractionalAsymmetricalArcFrequency{k}(1:continuationStartSinceFromIndex+1),...
%             plotDottedLineType,'Color','b','LineWidth',2);
    else % 后续相移进行部分节段复制延拓
        for p=2*continuationStartSinceFromIndex:-1:1
            if wrappedPhaseFractionalAsymmetricalArcFrequencyContinuation(p)<=phaseSegmentThreshold(k)
                fringeListContinuationFractionalAsymmetricalArcFrequency{k}=[fringeListContinuationFractionalAsymmetricalArcFrequency{k}(p:2*continuationStartSinceFromIndex),...
                    fringeListContinuationFractionalAsymmetricalArcFrequency{k}];
                plot(-(continuationStartSinceFromIndex-wrappedPhaseStartPeakIndex)-(2*continuationStartSinceFromIndex-p)-1:-(continuationStartSinceFromIndex-wrappedPhaseStartPeakIndex),...
                    fringeListContinuationFractionalAsymmetricalArcFrequency{k}(1:2*continuationStartSinceFromIndex+1-p+1),...
                    plotDottedLineType,'Color','b','LineWidth',2);
                break;
            end
        end
    end
    
    % 每步相移结束端需要延拓补齐的长度
    continuationStartLength(k)=length(fringeListContinuationFractionalAsymmetricalArcFrequency{k})-lengthOfSignal-continuationEndLength(k);

    title(sprintf('%d/%d Step of Original Fringe',k,moveNumPart));
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',[]);
    if k==moveNumPart
        set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    end
    
    xlim([-continuationStartSinceFromIndex*2+wrappedPhaseStartPeakIndex,lengthOfSignal+continuationEndMapDistance+lengthOfSignal+continuationEndMapDistance-wrappedPhaseEndPeakIndex]);
    ylim([-64-32 64+32]);grid on;
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end
% 起始端(左端)-<<<<-<<<<

%% 显示延拓后的条纹信号及其Hilbert变换
lengthOfContinuousSignal=zeros(moveNumPart,1);length(fringeListContinuationFractionalAsymmetricalArcFrequency{k});
htContinuousSignal=cell(size(fringeListContinuationFractionalAsymmetricalArcFrequency));
figure('name','Continuous Fringes & its Hilbert Transfrom (Frequency Modulated by Asymmetrical Arc Function)','NumberTitle','off');
for k=1:moveNumPart
    lengthOfContinuousSignal(k)=length(fringeListContinuationFractionalAsymmetricalArcFrequency{k});

    subplot(moveNumPart,1,k);
    % 条纹信号
    plot(0:lengthOfSignal-1,fringeListContinuationFractionalAsymmetricalArcFrequency{k}(continuationStartLength(k)+1:continuationStartLength(k)+lengthOfSignal),...
        plotLineType,'Color','b','LineWidth',1.0);hold on;
    % 延拓后的实际Hilbert变换结果
    htContinuousSignal{k}=imag(hilbert(fringeListContinuationFractionalAsymmetricalArcFrequency{k}));
    plot(0:lengthOfSignal-1,htContinuousSignal{k}(continuationStartLength(k)+1:continuationStartLength(k)+lengthOfSignal),...
        plotLineType,'Color','m','LineWidth',1.0);hold on;

    % 期望Hilbert变换结果
    plot(horizontalIndex,expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
    % 原始实际Hilbert变换结果
    plot(horizontalIndex,imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{k})),        plotLineType,'Color','g','LineWidth',1.0);

    title(sprintf('%d/%d Step of Continuous Fringes & its Hilbert Transfrom',k,moveNumPart));
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',[]);
    if k==moveNumPart
        set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
        legend('Original Fringe','Actual HT(Continuation)','Expected HT','Actual HT','Location','NorthEast','Orientation','Horizontal');
    end
    
    % 使用虚线绘制界外的信号及其Hilbert变换
    plot(-continuationStartLength(k):0,fringeListContinuationFractionalAsymmetricalArcFrequency{k}(1:continuationStartLength(k)+1),...
        plotDottedLineType,'Color','b','LineWidth',1.5);hold on;
    plot(lengthOfSignal-1:lengthOfContinuousSignal(k)-continuationStartLength(k)-1,fringeListContinuationFractionalAsymmetricalArcFrequency{k}(continuationStartLength(k)+lengthOfSignal:end),...
        plotDottedLineType,'Color','b','LineWidth',1.5);hold on;
    plot(-continuationStartLength(k):0,htContinuousSignal{k}(1:continuationStartLength(k)+1),...
        plotDottedLineType,'Color','m','LineWidth',1.5);hold on;
    plot(lengthOfSignal-1:lengthOfContinuousSignal(k)-continuationStartLength(k)-1,htContinuousSignal{k}(continuationStartLength(k)+lengthOfSignal:end),...
        plotDottedLineType,'Color','m','LineWidth',1.5);hold on;

    xlim([-continuationStartSinceFromIndex*2+wrappedPhaseStartPeakIndex,lengthOfSignal+continuationEndMapDistance+lengthOfSignal+continuationEndMapDistance-wrappedPhaseEndPeakIndex]);
    ylim([-64-32 64+32]);grid on;
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
    
    % 截取Hilbert变换区间结果
    htContinuousSignal{k}=htContinuousSignal{k}(continuationStartLength(k)+1:continuationStartLength(k)+lengthOfSignal);

end

% 显示延拓前后Hilbert变换期望值与实际值的条纹误差曲线
figure('name','Fringe Error between Expected HT and Actual HT (Including Continuation)','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);hold on;
    plot(horizontalIndex,imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{k}))-expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
%     htContinuousSignal=imag(hilbert(fringeListContinuationFractionalAsymmetricalArcFrequency{k}));
    plot(horizontalIndex,htContinuousSignal{k}-expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotLineType,'Color','m','LineWidth',1.0);
    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT',k,moveNumPart));
    if k==moveNumPart
        set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
        legend('Original Fringe Error','Continuous Fringe Error','Location','NorthEast','Orientation','Horizontal');
    end
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

%% --计算数步延拓后相移条纹的Hilbert域相位
% fringeListFractionalHilbertStepFrequency=HilbertPerRow(fringeListFractionalStepFrequency,moveNumPart);
for k=1:moveNumPart
    htContinuousSignal{k}=-htContinuousSignal{k};
end
wrappedPhaseFractionalHilbertContinuous=GetWrapPhaseWithHilbert(htContinuousSignal,moveNumPart);

% 显示空域相位误差、延拓Hilbert域相位误差
figure('name','Continuous Phase Error (Frequency Modulated by Step Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcFrequency-wrappedPhaseAllAsymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 延拓Hilbert域相位误差
wrappedErrorHTContinuous=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertContinuous-wrappedPhaseAllAsymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTContinuous,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Continuous Phase Error (Frequency Modulated by Step Function)');
legend('Space Phase Error','Continuous HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);


%% 对延拓信号进行两边翻转，进行Hilbert变换，再截取Hilbert变换区间结果
fringeListContinuationFlipFractionalAsymmetricalArcFrequency=cell(size(fringeListContinuationFractionalAsymmetricalArcFrequency));
htContinuousFlipSignal=cell(size(fringeListContinuationFractionalAsymmetricalArcFrequency));
figure('name','Continuous Fringes & its Hilbert Transfrom (Frequency Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 两侧均翻转
% for k=1:moveNumPart
%     flip0=fliplr(fringeListContinuationFractionalAsymmetricalArcFrequency{k});
%     fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}=[...
%         flip0(1)*2-flip0(2),...
%         flip0,...
%         fringeListContinuationFractionalAsymmetricalArcFrequency{k}(2:end),...
%         fringeListContinuationFractionalAsymmetricalArcFrequency{k}(end)*2-fringeListContinuationFractionalAsymmetricalArcFrequency{k}(end-1),...
%         flip0(1:end-1)];
%     
%     subplot(moveNumPart,1,k);hold on;
%     plot(0:lengthOfContinuousSignal(k),fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}(1:lengthOfContinuousSignal(k)+1),...
%         plotDottedLineType,'Color','b','LineWidth',1.0);hold on;
%     plot(lengthOfContinuousSignal(k):lengthOfContinuousSignal(k)*2,fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}(lengthOfContinuousSignal(k)+1:lengthOfContinuousSignal(k)*2+1),...
%         plotLineType,'Color','b','LineWidth',1.0);hold on;
%     plot(lengthOfContinuousSignal(k)*2:lengthOfContinuousSignal(k)*3-1,fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}(lengthOfContinuousSignal(k)*2+1:lengthOfContinuousSignal(k)*3),...
%         plotDottedLineType,'Color','b','LineWidth',1.0);hold on;
%     
%     title(sprintf('%d/%d Step of Continuous & Flip Fringes & its Hilbert Transfrom',k,moveNumPart));
%     %     xlim([-lengthOfContinuousSignal(k) lengthOfContinuousSignal(k)*2];
%     %     set(gca, 'XTick', xTick);set(gca, 'XTickLabel',[]);
%     %     if k==moveNumPart
%     %         set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
%     %     end
%     
%     htContinuousFlipSignal{k}=imag(hilbert(fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}));
%     plot(0:lengthOfContinuousSignal(k),htContinuousFlipSignal{k}(1:lengthOfContinuousSignal(k)+1),...
%         plotDottedLineType,'Color','m','LineWidth',1.0);hold on;
%     plot(lengthOfContinuousSignal(k):lengthOfContinuousSignal(k)*2,htContinuousFlipSignal{k}(lengthOfContinuousSignal(k)+1:lengthOfContinuousSignal(k)*2+1),...
%         plotLineType,'Color','m','LineWidth',1.0);hold on;
%     plot(lengthOfContinuousSignal(k)*2:lengthOfContinuousSignal(k)*3-1,htContinuousFlipSignal{k}(lengthOfContinuousSignal(k)*2+1:lengthOfContinuousSignal(k)*3),...
%         plotDottedLineType,'Color','m','LineWidth',1.0);hold on;
%     
%     % 截取区间
%     htContinuousFlipSignal{k}=htContinuousFlipSignal{k}(lengthOfContinuousSignal(k)+continuationStartLength(k)+1:lengthOfContinuousSignal(k)+continuationStartLength(k)+lengthOfSignal);
% end
% 仅于右侧翻转
for k=1:moveNumPart
    flip0=fliplr(fringeListContinuationFractionalAsymmetricalArcFrequency{k});
    fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}=[...
        fringeListContinuationFractionalAsymmetricalArcFrequency{k},...
        fringeListContinuationFractionalAsymmetricalArcFrequency{k}(end)*2-fringeListContinuationFractionalAsymmetricalArcFrequency{k}(end-1),...
        flip0(1:end-1)];
    
    subplot(moveNumPart,1,k);hold on;
    plot(0:lengthOfContinuousSignal(k),fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}(1:lengthOfContinuousSignal(k)+1),...
        plotLineType,'Color','b','LineWidth',1.0);hold on;
    plot(lengthOfContinuousSignal(k):lengthOfContinuousSignal(k)*2-1,fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}(lengthOfContinuousSignal(k)+1:lengthOfContinuousSignal(k)*2),...
        plotDottedLineType,'Color','b','LineWidth',1.0);hold on;
    
    title(sprintf('%d/%d Step of Continuous & Flip Fringes & its Hilbert Transfrom',k,moveNumPart));
    %     xlim([-lengthOfContinuousSignal(k) lengthOfContinuousSignal(k)*2];
    %     set(gca, 'XTick', xTick);set(gca, 'XTickLabel',[]);
    %     if k==moveNumPart
    %         set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    %     end
    
    htContinuousFlipSignal{k}=imag(hilbert(fringeListContinuationFlipFractionalAsymmetricalArcFrequency{k}));
    plot(0:lengthOfContinuousSignal(k),htContinuousFlipSignal{k}(1:lengthOfContinuousSignal(k)+1),...
        plotLineType,'Color','m','LineWidth',1.0);hold on;
    plot(lengthOfContinuousSignal(k):lengthOfContinuousSignal(k)*2-1,htContinuousFlipSignal{k}(lengthOfContinuousSignal(k)+1:lengthOfContinuousSignal(k)*2),...
        plotDottedLineType,'Color','m','LineWidth',1.0);hold on;
    
    % 截取区间
    htContinuousFlipSignal{k}=htContinuousFlipSignal{k}(continuationStartLength(k)+1:continuationStartLength(k)+lengthOfSignal);
end

% 显示[延拓后]的&[延拓及翻转后]的Hilbert变换对比结果
figure('name','Fringe Error between Expected HT and Actual HT (Including Continuation & Flip)','NumberTitle','off');
for k=1:moveNumPart
    lengthOfContinuousSignal(k)=length(fringeListContinuationFractionalAsymmetricalArcFrequency{k});

    subplot(moveNumPart,1,k);hold on;
    % 条纹信号
    plot(horizontalIndex,fringeListContinuationFractionalAsymmetricalArcFrequency{k}(continuationStartLength(k)+1:continuationStartLength(k)+lengthOfSignal),...
        plotLineType,'Color','b','LineWidth',1.0);hold on;
    % 延拓后的实际Hilbert变换结果
    plot(horizontalIndex,-htContinuousSignal{k},plotLineType,'Color','m','LineWidth',1.0);
    % 延拓及翻转后的实际Hilbert变换结果
    plot(horizontalIndex,htContinuousFlipSignal{k},plotLineType,'Color','r','LineWidth',1.5);

    % 期望Hilbert变换结果
    plot(horizontalIndex,expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);
    % 原始实际Hilbert变换结果
    plot(horizontalIndex,imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{k})),        plotLineType,'Color','g','LineWidth',1.0);

    title(sprintf('%d/%d Step of Fringe Error between Expected HT and Actual HT (Including Continuation & Flip)',k,moveNumPart));
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',[]);
    if k==moveNumPart
        set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
        legend('Original Fringe','Actual HT(Continuation)','Actual HT(Continuation & Flip)','Expected HT','Actual HT','Location','NorthEast','Orientation','Horizontal');
    end
    
    % 使用虚线绘制界外的信号及其Hilbert变换
    plot(-continuationStartLength(k):0,fringeListContinuationFractionalAsymmetricalArcFrequency{k}(1:continuationStartLength(k)+1),...
        plotDottedLineType,'Color','b','LineWidth',1.5);hold on;
    plot(lengthOfSignal-1:lengthOfContinuousSignal(k)-continuationStartLength(k)-1,fringeListContinuationFractionalAsymmetricalArcFrequency{k}(continuationStartLength(k)+lengthOfSignal:end),...
        plotDottedLineType,'Color','b','LineWidth',1.5);hold on;

    xlim([-continuationStartSinceFromIndex*2+wrappedPhaseStartPeakIndex,lengthOfSignal+continuationEndMapDistance+lengthOfSignal+continuationEndMapDistance-wrappedPhaseEndPeakIndex]);
    ylim([-64-32 64+32]);grid on;
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

% 显示[延拓后]的&[延拓及翻转后]的Hilbert变换与期望Hilbert变换结果的信号差
figure('name','Hilbert Phase Error (Continuous & Flip Fringes)','NumberTitle','off');
for k=1:moveNumPart
    subplot(moveNumPart,1,k);hold on;
    % 延拓后的实际Hilbert变换结果对比
    plot(horizontalIndex,-htContinuousSignal{k}-expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotLineType,'Color',[0,0,153]/255,'LineWidth',1.0);
    % 延拓及翻转后的实际Hilbert变换结果对比
    plot(horizontalIndex,htContinuousFlipSignal{k}-expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k},plotLineType,'Color','m','LineWidth',1.5);hold on;
%     % 原始实际Hilbert变换结果对比
%     plot(horizontalIndex,imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{k})),        plotLineType,'Color','g','LineWidth',1.0);

    title(sprintf('%d/%d Step of Hilbert Phase Error (Continuous & Flip Fringes)',k,moveNumPart));
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',[]);
    if k==moveNumPart
        set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
        legend('Actual HT(Continuation)','Actual HT(Continuation & Flip)','Location','NorthEast','Orientation','Horizontal');
    end

    xlim([-continuationStartSinceFromIndex*2+wrappedPhaseStartPeakIndex,lengthOfSignal+continuationEndMapDistance+lengthOfSignal+continuationEndMapDistance-wrappedPhaseEndPeakIndex]);
%     ylim([-64-32 64+32]);
    grid on;
%     set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end


%% --计算数步延拓及翻转后相移条纹的Hilbert域相位
% fringeListFractionalHilbertStepFrequency=HilbertPerRow(fringeListFractionalStepFrequency,moveNumPart);
for k=1:moveNumPart
    htContinuousFlipSignal{k}=-htContinuousFlipSignal{k};
end
wrappedPhaseFractionalHilbertContinuousFlip=GetWrapPhaseWithHilbert(htContinuousFlipSignal,moveNumPart);

% 显示空域相位误差、延拓翻转Hilbert域相位误差
figure('name','Continuous & Flip Phase Error (Frequency Modulated by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
wrappedErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAsymmetricalArcFrequency-wrappedPhaseAllAsymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorSpace,plotLineType,'LineWidth',0.5,'MarkerSize',2);hold on;
% 延拓Hilbert域相位误差
wrappedErrorHTContinuousFlip=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertContinuousFlip-wrappedPhaseAllAsymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound);
plot(wrappedErrorHTContinuousFlip,   plotLineType,'LineWidth',0.5,'MarkerSize',2);
title('Continuous & Flip Phase Error (Frequency Modulated by Step Function)');
legend('Space Phase Error','Continuous & Flip HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

end

%% @_@{对非对称连续调制频率信号[非完整周期]进行逐相位周期Hilbert变换分析}********************
% 非对称连续调制频率[非完整周期]逐相位周期Hilbert变换分析标记
if asymmetricalArcFrequencyModulateEachPeriodFlag==1
asymmetricalArcFactor=1.25;
lengthEnd=lengthOfSignal;
    
%% --生成24幅非对称连续调制频率全部条纹图像并从中抽取出数步相移条纹图像
moveNumAll=24;
fringeListAllAsymmetricalArcFrequency=cell(24,1);
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllAsymmetricalArcFrequency{k}=(255.0/2*(cos(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal)))/period)+0)/2);
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
for k=1:moveNumAll
    fringeListAllAsymmetricalArcFrequency{k}=fringeListAllAsymmetricalArcFrequency{k}(eachPeriodFringeStartIndex(1):eachPeriodFringeStartIndex(end)-1);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListFractionalAsymmetricalArcFrequency=SelectNStepFring(fringeListAllAsymmetricalArcFrequency,moveNumPart);
% 期望的Hilbert变换
expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency=cell(size(fringeListFractionalAsymmetricalArcFrequency));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=255.0/2*(sin(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal)))/period))/2;
    expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}=expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{k}(eachPeriodFringeStartIndex(1):eachPeriodFringeStartIndex(end)-1);
end

% 更改参数值
lengthOfSignal=eachPeriodFringeStartIndex(end)-eachPeriodFringeStartIndex(1);
lengthEnd=lengthOfSignal;
eachPeriodFringeStartIndex=eachPeriodFringeStartIndex-eachPeriodFringeStartIndex(1)+1;
eachPeriodFringeStartIndex(end)=eachPeriodFringeStartIndex(end)-1;

%% --计算数步相移条纹的空域相位
wrappedPhaseFractionalAsymmetricalArcFrequency=GetWrapPhase(fringeListFractionalAsymmetricalArcFrequency,moveNumPart);

%% --计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAsymmetricalArcFrequency=HilbertPerRow(fringeListFractionalAsymmetricalArcFrequency,moveNumPart);
wrappedPhaseFractionalHilbertAsymmetricalArcFrequency=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAsymmetricalArcFrequency,moveNumPart);

%% --显示非对称连续函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示非对称连续
symmetricalArcFrequencyFunction=(0:lengthOfSignal-1)+100*sin(2*pi*(0:lengthOfSignal-1)/(asymmetricalArcFactor*lengthOfSignal));
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

% 显示理想折叠相位
figure('name','Wrapped Phase','NumberTitle','off');
plot(wrappedPhaseAllAsymmetricalArcFrequency,plotLineType,'LineWidth',1.0,'MarkerSize',2);
title('Wrapped Phase');
xlim([0,lengthOfSignal-1]);grid on;% ylim([-1,1]);
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 按虚实线逐周期显示第1步条纹信号、Hilbert变换
figure('name','Each Period of First Step Fringe Modulated by Asymmetrical Arc Function','NumberTitle','off');
% htFlip=zeros(1,length(fringeListFractionalAsymmetricalArcFrequency{1}));
for p=1:periodNumber-1
    if mod(p,2)==1
        % 逐周期条纹信号
        plot(eachPeriodFringeStartIndex(p)-1:eachPeriodFringeStartIndex(p+1)-1,fringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)),...
            plotLineType,'Color',[0,128,0]/255,'LineWidth',1.0);hold on;   
        % 全局Hilbert变换
        if p==1
            plot(0:lengthEnd-1,imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{1})),        plotLineType,'Color','b','LineWidth',1.0);hold on;
        end

        % 期望Hilbert变换
        plot(eachPeriodFringeStartIndex(p)-1:eachPeriodFringeStartIndex(p+1)-1,expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)),...
            plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
        % 实际Hilbert变换结果
        plot(eachPeriodFringeStartIndex(p)-1:eachPeriodFringeStartIndex(p+1)-1-2,imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)-2))),...
            plotLineType,'Color',[0,206,209]/255,'LineWidth',1.0);hold on;
        
        % 镜像
        htFlip=imag(hilbert(...
            [fringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)-1),...
            fliplr(fringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)-2))]));
        plot(eachPeriodFringeStartIndex(p)-1:eachPeriodFringeStartIndex(p+1)-1-2,htFlip(1:eachPeriodFringeStartIndex(p+1)-1-eachPeriodFringeStartIndex(p)),...
            plotLineType,'Color',[255,0,0]/255,'LineWidth',1.0);hold on;
    else
        % 逐周期条纹信号
        plot(eachPeriodFringeStartIndex(p)-1:eachPeriodFringeStartIndex(p+1)-1,fringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)),...
            plotDottedLineType,'Color',[0,128,0]/255,'LineWidth',1.5);hold on;
        % 期望Hilbert变换
        plot(eachPeriodFringeStartIndex(p)-1:eachPeriodFringeStartIndex(p+1)-1,expectedHilbertOfFringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)),...
            plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
        % 实际Hilbert变换结果
        plot(eachPeriodFringeStartIndex(p)-1:eachPeriodFringeStartIndex(p+1)-1-2,imag(hilbert(fringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)-2))),...
            plotLineType,'Color',[0,206,209]/255,'LineWidth',1.0);hold on;
        
        % 镜像
        htFlip=imag(hilbert(...
            [fringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)-1),...
            fliplr(fringeListFractionalAsymmetricalArcFrequency{1}(eachPeriodFringeStartIndex(p):eachPeriodFringeStartIndex(p+1)-2))]));
        plot(eachPeriodFringeStartIndex(p)-1:eachPeriodFringeStartIndex(p+1)-1-2,htFlip(1:eachPeriodFringeStartIndex(p+1)-1-eachPeriodFringeStartIndex(p)),...
            plotLineType,'Color',[255,0,0]/255,'LineWidth',1.0);hold on;
    end
    
    if p==1
        legend('Original Fringe','All Actual HT','All Expected HT','Actual HT (Per Period)','Truncated Flip HT (Per Period)','Location','SouthOutside','Orientation','Horizontal');
    end
end
title('Each Period of First Step Fringe Modulated by Asymmetrical Arc Function');
xlim([0,lengthOfSignal-1]);ylim([-64-16 64+16]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

end



