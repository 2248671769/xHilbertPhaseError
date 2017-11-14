% xFringeContinuation.m
% 徐文宇，延拓模块由杨睦尧提供，20171110
% 受非周期信号频谱泄露的影响，Hilbert域相位误差被放大。
% 根据折叠相位边界两端的周期缺失对条纹信号进行延拓，再计算Hilbert域相位，并与空域相位相加求平均值，实现相位误差补偿
% ver:---
% 操作步骤：
% 1.求取[24步理想空域相位wp24]与[4步空域相位wp4]
% 2.求取[归一化后的4步Hilbert域相位wpnHT4]
% 3.计算[wp4]与[wpnHT4]的[均值相位meanOfwp4AndwpnHT4],并计算该值与[wp24]的相位误差
% 4.对归一化后的条纹信号的边界两端分别进行周期延拓
% 5.计算[周期延拓后的Hilbert域相位wpncHT4]，计算[wp4]与[wpncHT4]的[均值相位meanOfwp4AndwpncHT4]，并计算该值与[wp24]的相位误差
% 6.比较上述两个相位误差
close all;clear;

%% {设置基本参数}*************************************************
% 图像条纹参数
width=1024; height=800; period=64;

%% -信号延拓标记
% 阶跃式与非对称连续弧段并行幅度调制标记
amplitudeModulatedFlag=0;
% 全调制：阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制
modulatedGroupFlag=1;

% 信号范围
numOfPeriods=16;
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
% xTickContinuation & xTickLabelkContinuation
xTickContinuation=zeros(1,numOfPeriods+1);
xTickLabelContinuation=cell(1,numOfPeriods+1);
for xt=0:numOfPeriods
    xTickContinuation(xt+1)=floor(xt*period*2); xTickLabelContinuation{xt+1}=num2str(xTickContinuation(xt+1));
end
xTickContinuation(end)=lengthOfSignal*2-1; xTickLabelContinuation{end}=num2str(xTickContinuation(end));

% 相位误差显示有效区间
upPhaseErrorBound=2; bottomPhaseErrorBound=-2;

% plot画线类型
plotLineType='';        % '' 实线
plotDottedLineType=':'; % ':'虚线

%% {生成各类调制信号}********************************************
% -单位半幅度条纹信号
moveNumAll=24;
fringeListAllUnit=cell(moveNumAll,1);
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllUnit{k} = floor(255*0.5*(cos(((0:lengthOfSignal-1)-sf)/period*2*pi)+1)/2);
end

%% -阶跃式与非对称连续弧段并行幅度调制信号
fringeListAllAmplitudeModulated=cell(moveNumAll,1);
% 阶跃函数
filterStepAmplitude=[ones(1,lengthOfSignal/2) 2*ones(1,lengthOfSignal/2)];
% 非对称连续弧段函数
arcCenter=0.625;
filterAsymmetricalArcAmplitude=1-1.0*((0:lengthOfSignal-1)/lengthOfSignal-arcCenter).^2;
for k=1:moveNumAll
    fringeListAllAmplitudeModulated{k}=fringeListAllUnit{k}.*filterStepAmplitude.*filterAsymmetricalArcAmplitude;
end

%% -阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制信号
fringeListAllModulatedGroup=cell(moveNumAll,1);
for k=1:moveNumAll
    sf=-period*(k-1)/moveNumAll;
    fringeListAllModulatedGroup{k}=255.0/2*(cos(2*pi*((0:lengthOfSignal-1)-sf+100*sin(2*pi*(0:lengthOfSignal-1)/(1.25*lengthOfSignal)))/period)+1)/2 ...
         .*filterStepAmplitude.*filterAsymmetricalArcAmplitude;
end

%% -提取数步相称的条纹信号
moveNumPart=4;
% -单位半幅度条纹信号
fringeListFractionalUnit=SelectNStepFring(fringeListAllUnit,moveNumPart);
% -阶跃式与非对称连续弧段并行幅度调制
fringeListFractionalAmplitudeModulated=SelectNStepFring(fringeListAllAmplitudeModulated,moveNumPart);
% -阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制
fringeListFractionalModulatedGroup=SelectNStepFring(fringeListAllModulatedGroup,moveNumPart);

%% -计算理想空域相位
% 单位条纹
wrappedPhaseAllUnit=GetWrapPhase(fringeListAllUnit,moveNumAll);
ABAllUnit=GetABNormal(fringeListAllUnit,wrappedPhaseAllUnit);
% 阶跃式与非对称连续弧段并行幅度调制
wrappedPhaseAllAmplitudeModulated=GetWrapPhase(fringeListAllAmplitudeModulated,moveNumAll);
ABAllAmplitudeModulated=GetABNormal(fringeListAllAmplitudeModulated,wrappedPhaseAllAmplitudeModulated);
% 阶跃式与非对称连续弧段并行幅度调制
wrappedPhaseAllModulatedGroup=GetWrapPhase(fringeListAllModulatedGroup,moveNumAll);
ABAllModulatedGroup=GetABNormal(fringeListAllModulatedGroup,wrappedPhaseAllModulatedGroup);

%% -计算数步相移(单位半幅度条纹信号)的空域相位
wrappedPhaseFractionalUnit=GetWrapPhase(fringeListFractionalUnit,moveNumPart);

%% -显示条纹信号
for k=1:moveNumPart
    figure('name',sprintf('%d/%d Step of Fringe Signal (Modulated Group)',k,moveNumPart),'NumberTitle','off');
    % -单位半幅度条纹信号
    plot(fringeListFractionalUnit{k},':','LineWidth',1.5);hold on;
    % -阶跃式与非对称连续弧段并行幅度调制信号
    plot(fringeListFractionalAmplitudeModulated{k},'b');
    % -阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制信号
    plot(fringeListFractionalModulatedGroup{k},'LineWidth',2);
    title(sprintf('%d/%d Step of Fringe Signal (Modulated Group)',k,moveNumPart));
    legend('Unit Fringe','Amplitude Modulated','Modulated Group');
    xlim([0,lengthOfSignal-1]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
end

%% *****************************************************************
if amplitudeModulatedFlag==1
%% [Hilbert]-[延拓]-(阶跃式与非对称连续弧段并行幅度调制信号)*
%% -计算数步相移(阶跃式与非对称连续弧段并行幅度调制)的空域相位
wrappedPhaseFractionalAmplitudeModulated=GetWrapPhase(fringeListFractionalAmplitudeModulated,moveNumPart);
ABFractionalAmplitudeModulated=GetABNormal(fringeListFractionalAmplitudeModulated,wrappedPhaseFractionalAmplitudeModulated);

%% -[Hilbert]-计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertAmplitudeModulated=HilbertPerRow(fringeListFractionalAmplitudeModulated,moveNumPart);
wrappedPhaseFractionalHilbertAmplitudeModulated=GetWrapPhaseWithHilbert(fringeListFractionalHilbertAmplitudeModulated,moveNumPart);

% %% -显示Hilbert变换
% for k=1:moveNumPart
%     figure('name',sprintf('%d/%d Hilbert Transform of Fringe Signal (Amplitude Modulated)',k,moveNumPart),'NumberTitle','off');
%     % -阶跃式与非对称连续弧段并行幅度调制信号
%     plot(fringeListFractionalAmplitudeModulated{k});hold on;
%     % -Hilbert变换
%     plot(imag(hilbert(fringeListFractionalAmplitudeModulated{k})),'LineWidth',1.5);
%     title(sprintf('%d/%d Hilbert Transform of Fringe Signal (Amplitude Modulated)',k,moveNumPart));
%     legend('Amplitude Modulated','Hilbert Transform','Location','SouthWest');
%     xlim([0,lengthOfSignal-1]);grid on;
%     set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
%     set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% end

%% -归一化处理(fringe-A)/B及归一化后的空域相位
fringeListAllNormalizationAmplitudeModulated=NormalizeFringe(fringeListAllAmplitudeModulated,moveNumAll,ABAllAmplitudeModulated);
wrappedPhaseAllNormalizationAmplitudeModulated=GetWrapPhase(fringeListAllNormalizationAmplitudeModulated,moveNumAll);

%% -归一化处理(fringe-A)/B及归一化后的空域相位
fringeListFractionalNormalizationAmplitudeModulated=NormalizeFringe(fringeListFractionalAmplitudeModulated,moveNumPart,ABFractionalAmplitudeModulated);
wrappedPhaseFractionalNormalizationAmplitudeModulated=GetWrapPhase(fringeListFractionalNormalizationAmplitudeModulated,moveNumPart);
for k=1:moveNumPart
    figure('name',sprintf('%d/%d Step of Normalization of Fringe Signal (Amplitude Modulated)',k,moveNumPart),'NumberTitle','off');
    % -阶跃式与非对称连续弧段并行幅度调制信号
    subplot(2,1,1);
    plot(fringeListFractionalAmplitudeModulated{k});hold on;
    title(sprintf('%d/%d Step of Fringe Signal (Amplitude Modulated)',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);ylim([0,256]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
    % -归一化处理(fringe-A)/B
    subplot(2,1,2);
    plot(fringeListFractionalNormalizationAmplitudeModulated{k},'m','LineWidth',1.5);hold on;
    sf=-period*(k-1)/moveNumPart;
    plot(cos(((0:lengthOfSignal-1)-sf)/period*2*pi),plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);
    title(sprintf('%d/%d Step of Normalization of Fringe Signal (Amplitude Modulated)',k,moveNumPart));
    legend('Normalizaion','Ideal Cosine','Location','SouthOutside','Orientation','Horizontal');
    xlim([0,lengthOfSignal-1]);ylim([-1.05,1.05]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

% %% -{延拓}
% fringeListFractionalContinuationAmplitudeModulated=cell(size(fringeListFractionalAmplitudeModulated));
% for k=1:moveNumPart
%     sf=-period*(k-1)/moveNumPart;
%     fringeListAllUnit{k} = floor(255*0.5*(cos(((0:lengthOfSignal-1)-sf)/period*2*pi)+1)/2);
%     fringeListFractionalContinuationAmplitudeModulated{k}=[ ...
%         fringeListFractionalNormalizationAmplitudeModulated{k}, ...
%         cos(((lengthOfSignal)-sf)/period*2*pi), ...
%         (mod(k,2)*(mod(k,2)+1)-1)*fliplr(fringeListFractionalNormalizationAmplitudeModulated{k}(1:lengthOfSignal))
%         ];
%     
%     figure('name',sprintf('%d/%d Step of Continuation of Fringe Signal (Amplitude Modulated)',k,moveNumPart),'NumberTitle','off');
%     % -阶跃式与非对称连续弧段并行幅度调制信号
%     plot(fringeListFractionalNormalizationAmplitudeModulated{k});hold on;
%     % -Hilbert变换
%     plot(imag(hilbert(fringeListFractionalNormalizationAmplitudeModulated{k})),plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);
%     % -延拓
%     plot(lengthOfSignal:lengthOfSignal*2, fringeListFractionalContinuationAmplitudeModulated{k}(lengthOfSignal:lengthOfSignal*2),'m','LineWidth',1.0);
%     % -Hilbert变换(延拓)
%     plot(imag(hilbert(fringeListFractionalContinuationAmplitudeModulated{k})),'LineWidth',1.5);
%     title(sprintf('%d/%d Step of Continuation of Fringe Signal (Amplitude Modulated)',k,moveNumPart));
%     legend('Amplitude Modulated','Hilbert Transform','Continuation','Hilbert Transform(Continuation)','Location','NorthEast');
%     xlim([0,lengthOfSignal*2-1]);grid on;
%     set(gca, 'XTick', xTickContinuation);set(gca, 'XTickLabel',xTickLabelContinuation);
% %     set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 
% end

%% -[Hilbert]-计算归一化后的数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertNormalizationAmplitudeModulated=HilbertPerRow(fringeListFractionalNormalizationAmplitudeModulated,moveNumPart);
wrappedPhaseFractionalHilbertNormalizationAmplitudeModulated=GetWrapPhaseWithHilbert(fringeListFractionalHilbertNormalizationAmplitudeModulated,moveNumPart);

%% -归一化后的平均相位=(归一化后的空域相位+归一化的Hilbert域相位)/2
wrappedPhaseFractionalMeanNormalization=(wrappedPhaseFractionalNormalizationAmplitudeModulated+wrappedPhaseFractionalHilbertNormalizationAmplitudeModulated)/2.0;

%% -显示归一化后的相位
figure('name','Phase after Normalization (Amplitude Modulated)','NumberTitle','off');
% 24步相移归一化后的理想空域相位
plot(wrappedPhaseAllNormalizationAmplitudeModulated, plotDottedLineType,'LineWidth',1.5);hold on;
% 数步相移归一化后的空域相位
plot(wrappedPhaseFractionalNormalizationAmplitudeModulated, plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 数步相移归一化后的Hilbert域相位
plot(wrappedPhaseFractionalHilbertNormalizationAmplitudeModulated,       plotLineType,'Color','g','LineWidth',1.0);hold on;
% 数步相移归一化后的平均相位
plot(wrappedPhaseFractionalMeanNormalization,       plotLineType,'Color','m','LineWidth',1.0);
title('Phase Error after Normalization (Amplitude Modulated)');
legend( ...
    sprintf('Space Phase after Normalization (Ideal %d Steps)',moveNumAll), ...
    sprintf('Space Phase after Normalization (%d Steps)',moveNumPart), ...
    sprintf('Hilbert Phase after Normalization (%d Steps Amplitude Modulated)',moveNumPart), ...
    sprintf('Mean Phase after Normalization (%d Steps Amplitude Modulated)',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% -显示归一化后的相位误差
figure('name','Phase Error after Normalization (Amplitude Modulated)','NumberTitle','off');
% 数步相移归一化后的空域相位误差
phaseErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalNormalizationAmplitudeModulated-wrappedPhaseAllNormalizationAmplitudeModulated,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseErrorSpace, plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 数步相移归一化后的Hilbert域相位误差
phaseErrorHilbert=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertNormalizationAmplitudeModulated-wrappedPhaseAllNormalizationAmplitudeModulated,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseErrorHilbert,     plotLineType,'Color','g','LineWidth',1.0);hold on;
% 数步相移归一化后的平均相位误差
phaseErrorMean=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalMeanNormalization-wrappedPhaseAllNormalizationAmplitudeModulated,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseErrorMean,        plotLineType,'Color','m','LineWidth',1.0);
title(sprintf('Phase Error after Normalization (%d Steps Amplitude Modulated)',moveNumPart));
legend( ...
    'Space Phase Error after Normalization', ...
    'Hilbert Phase Error after Normalization', ...
    'Mean  Phase Error after Normalization');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/平均相位误差的平均值、峰值与均方根
% 数步相移归一化后的空域相位误差
fprintf('------------Phase Error after Normalization (%d Steps Amplitude Modulated)-------------\n',moveNumPart);
fprintf('          Mean of Space Phase Error: %+f\n',mean(phaseErrorSpace));
fprintf('  Max positive of Space Phase Error: %+f\n',max(phaseErrorSpace));
fprintf('  Max negative of Space Phase Error: %+f\n',min(phaseErrorSpace));
fprintf('          RMSE of Space Phase Error: %+f\n',sqrt(sum((phaseErrorSpace  -mean(phaseErrorSpace))  .^2))/lengthOfSignal);
% 数步相移归一化后的Hilbert域相位误差
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(phaseErrorHilbert));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max(phaseErrorHilbert));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min(phaseErrorHilbert));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((phaseErrorHilbert-mean(phaseErrorHilbert)).^2))/lengthOfSignal);
% 数步相移归一化后的平均相位误差
fprintf('           Mean of Mean Phase Error: %+f\n',mean(phaseErrorMean));
fprintf('   Max positive of Mean Phase Error: %+f\n',max(phaseErrorMean));
fprintf('   Max negetive of Mean Phase Error: %+f\n',min(phaseErrorMean));
fprintf('           RMSE of Mean Phase Error: %+f\n',sqrt(sum((phaseErrorMean   -mean(phaseErrorMean))   .^2))/lengthOfSignal);
end

%% *****************************************************************
if modulatedGroupFlag==1
%% [Hilbert]-[延拓]-(全调制：阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制)*
%% -计算数步相移(全调制)的空域相位
wrappedPhaseFractionalModulatedGroup=GetWrapPhase(fringeListFractionalModulatedGroup,moveNumPart);
ABFractionalModulatedGroup=GetABNormal(fringeListFractionalModulatedGroup,wrappedPhaseFractionalModulatedGroup);

%% -[Hilbert]-计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertModulatedGroup=HilbertPerRow(fringeListFractionalModulatedGroup,moveNumPart);
wrappedPhaseFractionalHilbertModulatedGroup=GetWrapPhaseWithHilbert(fringeListFractionalHilbertModulatedGroup,moveNumPart);

% %% -显示Hilbert变换
% for k=1:moveNumPart
%     figure('name',sprintf('%d/%d Hilbert Transform of Fringe Signal (Amplitude Modulated)',k,moveNumPart),'NumberTitle','off');
%     % -阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制
%     plot(fringeListFractionalModulatedGroup{k});hold on;
%     % -Hilbert变换
%     plot(imag(hilbert(fringeListFractionalModulatedGroup{k})),'LineWidth',1.5);
%     title(sprintf('%d/%d Hilbert Transform of Fringe Signal (Amplitude Modulated)',k,moveNumPart));
%     legend('Amplitude Modulated','Hilbert Transform','Location','SouthWest');
%     xlim([0,lengthOfSignal-1]);grid on;
%     set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
%     set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% end

%% -归一化处理(fringe-A)/B及归一化后的空域相位
fringeListAllNormalizationModulatedGroup=NormalizeFringe(fringeListAllModulatedGroup,moveNumAll,ABAllModulatedGroup);
wrappedPhaseAllNormalizationModulatedGroup=GetWrapPhase(fringeListAllNormalizationModulatedGroup,moveNumAll);

%% -归一化处理(fringe-A)/B及归一化后的空域相位
fringeListFractionalNormalizationModulatedGroup=NormalizeFringe(fringeListFractionalModulatedGroup,moveNumPart,ABFractionalModulatedGroup);
wrappedPhaseFractionalNormalizationModulatedGroup=GetWrapPhase(fringeListFractionalNormalizationModulatedGroup,moveNumPart);
for k=1:moveNumPart
    figure('name',sprintf('%d/%d Step of Normalization of Fringe Signal (Amplitude Modulated)',k,moveNumPart),'NumberTitle','off');
    % -阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制
    subplot(2,1,1);
    plot(fringeListFractionalModulatedGroup{k});hold on;
    title(sprintf('%d/%d Step of Fringe Signal (Amplitude Modulated)',k,moveNumPart));
    xlim([0,lengthOfSignal-1]);ylim([0,256]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
    set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
    % -归一化处理(fringe-A)/B
    subplot(2,1,2);
    plot(fringeListFractionalNormalizationModulatedGroup{k},'m','LineWidth',1.5);hold on;
    sf=-period*(k-1)/moveNumPart;
    plot(cos(((0:lengthOfSignal-1)-sf)/period*2*pi),plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);
    title(sprintf('%d/%d Step of Normalization of Fringe Signal (Amplitude Modulated)',k,moveNumPart));
    legend('Normalizaion','Ideal Cosine','Location','SouthOutside','Orientation','Horizontal');
    xlim([0,lengthOfSignal-1]);ylim([-1.05,1.05]);grid on;
    set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
end

%% -{延拓}
fringeListFractionalContinuationModulatedGroup=cell(size(fringeListFractionalModulatedGroup));
for k=1:moveNumPart
    sf=-period*(k-1)/moveNumPart;
    fringeListAllUnit{k} = floor(255*0.5*(cos(((0:lengthOfSignal-1)-sf)/period*2*pi)+1)/2);
    fringeListFractionalContinuationModulatedGroup{k}=[ ...
        fringeListFractionalNormalizationModulatedGroup{k}, ...
        cos(2*pi*((lengthOfSignal)-sf+100*sin(2*pi*(lengthOfSignal)/(1.25*lengthOfSignal)))/period), ...
        (mod(k,2)*(mod(k,2)+1)-1)*fliplr(fringeListFractionalNormalizationModulatedGroup{k}(1:lengthOfSignal))
        ];
    
    figure('name',sprintf('%d/%d Step of Continuation of Fringe Signal (Amplitude Modulated)',k,moveNumPart),'NumberTitle','off');
    % -阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制
    plot(fringeListFractionalNormalizationModulatedGroup{k});hold on;
    % -Hilbert变换
    plot(imag(hilbert(fringeListFractionalNormalizationModulatedGroup{k})),plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);
    % -延拓
    plot(lengthOfSignal:lengthOfSignal*2, fringeListFractionalContinuationModulatedGroup{k}(lengthOfSignal:lengthOfSignal*2),'m','LineWidth',1.0);
    % -Hilbert变换(延拓)
    plot(imag(hilbert(fringeListFractionalContinuationModulatedGroup{k})),'LineWidth',1.5);
    title(sprintf('%d/%d Step of Continuation of Fringe Signal (Amplitude Modulated)',k,moveNumPart));
    legend('Amplitude Modulated','Hilbert Transform','Continuation','Hilbert Transform(Continuation)','Location','NorthEast');
    xlim([0,lengthOfSignal*2-1]);grid on;
    set(gca, 'XTick', xTickContinuation);set(gca, 'XTickLabel',xTickLabelContinuation);
%     set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

end

%% -[Hilbert]-计算归一化后的数步相移条纹的Hilbert变换与Hilbert域相位
fringeListFractionalHilbertNormalizationModulatedGroup=HilbertPerRow(fringeListFractionalNormalizationModulatedGroup,moveNumPart);
wrappedPhaseFractionalHilbertNormalizationModulatedGroup=GetWrapPhaseWithHilbert(fringeListFractionalHilbertNormalizationModulatedGroup,moveNumPart);

%% -归一化后的平均相位=(归一化后的空域相位+归一化的Hilbert域相位)/2
wrappedPhaseFractionalMeanNormalization=(wrappedPhaseFractionalNormalizationModulatedGroup+wrappedPhaseFractionalHilbertNormalizationModulatedGroup)/2.0;

%% -显示归一化后的相位
figure('name','Phase after Normalization (Modulated Group)','NumberTitle','off');
% 24步相移归一化后的理想空域相位
plot(wrappedPhaseAllNormalizationModulatedGroup, plotDottedLineType,'LineWidth',1.5);hold on;
% 数步相移归一化后的空域相位
plot(wrappedPhaseFractionalNormalizationModulatedGroup, plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 数步相移归一化后的Hilbert域相位
plot(wrappedPhaseFractionalHilbertNormalizationModulatedGroup,       plotLineType,'Color','g','LineWidth',1.0);hold on;
% 数步相移归一化后的平均相位
plot(wrappedPhaseFractionalMeanNormalization,       plotLineType,'Color','m','LineWidth',1.0);
title('Phase Error after Normalization (Modulated Group)');
legend( ...
    sprintf('Space Phase after Normalization (Ideal %d Steps)',moveNumAll), ...
    sprintf('Space Phase after Normalization (%d Steps)',moveNumPart), ...
    sprintf('Hilbert Phase after Normalization (%d Steps Modulated Group)',moveNumPart), ...
    sprintf('Mean Phase after Normalization (%d Steps Modulated Group)',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% -显示归一化后的相位误差
figure('name','Phase Error after Normalization (Modulated Group)','NumberTitle','off');
% 数步相移归一化后的空域相位误差
phaseErrorSpace=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalNormalizationModulatedGroup-wrappedPhaseAllNormalizationModulatedGroup,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseErrorSpace, plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 数步相移归一化后的Hilbert域相位误差
phaseErrorHilbert=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalHilbertNormalizationModulatedGroup-wrappedPhaseAllNormalizationModulatedGroup,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseErrorHilbert,     plotLineType,'Color','g','LineWidth',1.0);hold on;
% 数步相移归一化后的平均相位误差
phaseErrorMean=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalMeanNormalization-wrappedPhaseAllNormalizationModulatedGroup,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseErrorMean,        plotLineType,'Color','m','LineWidth',1.0);
title(sprintf('Phase Error after Normalization (%d Steps Modulated Group)',moveNumPart));
legend( ...
    'Space Phase Error after Normalization', ...
    'Hilbert Phase Error after Normalization', ...
    'Mean  Phase Error after Normalization');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/平均相位误差的平均值、峰值与均方根
% 数步相移归一化后的空域相位误差
fprintf('------------Phase Error after Normalization (%d Steps Modulated Group)-------------\n',moveNumPart);
fprintf('          Mean of Space Phase Error: %+f\n',mean(phaseErrorSpace));
fprintf('  Max positive of Space Phase Error: %+f\n',max(phaseErrorSpace));
fprintf('  Max negative of Space Phase Error: %+f\n',min(phaseErrorSpace));
fprintf('          RMSE of Space Phase Error: %+f\n',sqrt(sum((phaseErrorSpace  -mean(phaseErrorSpace))  .^2))/lengthOfSignal);
% 数步相移归一化后的Hilbert域相位误差
fprintf('        Mean of Hilbert Phase Error: %+f\n',mean(phaseErrorHilbert));
fprintf('Max positive of Hilbert Phase Error: %+f\n',max(phaseErrorHilbert));
fprintf('Max negetive of Hilbert Phase Error: %+f\n',min(phaseErrorHilbert));
fprintf('        RMSE of Hilbert Phase Error: %+f\n',sqrt(sum((phaseErrorHilbert-mean(phaseErrorHilbert)).^2))/lengthOfSignal);
% 数步相移归一化后的平均相位误差
fprintf('           Mean of Mean Phase Error: %+f\n',mean(phaseErrorMean));
fprintf('   Max positive of Mean Phase Error: %+f\n',max(phaseErrorMean));
fprintf('   Max negetive of Mean Phase Error: %+f\n',min(phaseErrorMean));
fprintf('           RMSE of Mean Phase Error: %+f\n',sqrt(sum((phaseErrorMean   -mean(phaseErrorMean))   .^2))/lengthOfSignal);

end

return

%% {Hilbert}-{延拓}-(阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制)*

%% -计算数步相移(阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制)的空域相位
wrappedPhaseFractionalModulatedGroup=GetWrapPhase(fringeListFractionalModulatedGroup,moveNumPart);
% 
% %% -计算数步相移条纹的Hilbert变换与Hilbert域相位
% fringeListFractionalHilbertModulatedGroup=HilbertPerRow(fringeListFractionalModulatedGroup,moveNumPart);
% wrappedPhaseFractionalHilbertModulatedGroup=GetWrapPhaseWithHilbert(fringeListFractionalHilbertModulatedGroup,moveNumPart);

%% {显示图表}*****************************************************

%% -显示折叠相位
figure('name','Wrapped Phase (Amplitude and Frequency Modulated Group)','NumberTitle','off');
% -理想空域相位
plot(wrappedPhaseAllModulatedGroup,           plotLineType,'LineWidth',1.0);hold on;
% 数步相移(单位半幅度条纹信号)的空域相位
plot(wrappedPhaseFractionalUnit, plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 数步相移(阶跃式与非对称连续弧段并行幅度调制)的空域相位
plot(wrappedPhaseFractionalAmplitudeModulated,plotLineType,'Color','g','LineWidth',1.0);hold on;
% 数步相移(阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制)的空域相位
plot(wrappedPhaseFractionalModulatedGroup,plotLineType,'Color','m','LineWidth',1.0);
title('Phase Error (Amplitude and Frequency Modulated Group)');
legend('Ideal Space Phase', ...
    sprintf('Space Phase (%d Steps Unit)',moveNumPart), ...
    sprintf('Space Phase (%d Steps Amplitude Modulated)',moveNumPart), ...
    sprintf('Space Phase (%d Steps Modulated Group)',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

%% -显示相位误差
figure('name','Phase Error (Amplitude and Frequency Modulated Group)','NumberTitle','off');
% 数步相移(单位半幅度条纹信号)的空域相位误差
phaseError1=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalUnit-wrappedPhaseAllModulatedGroup,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseError1, plotDottedLineType,'Color',[0,0,153]/255,'LineWidth',1.5);hold on;
% 数步相移(阶跃式与非对称连续弧段并行幅度调制)的空域相位误差
phaseError2=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalAmplitudeModulated-wrappedPhaseAllModulatedGroup,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseError2,plotLineType,'Color','g','LineWidth',1.0);hold on;
% 数步相移(阶跃式与非对称连续弧段并行幅度调制及非对称连续弧段频率调制)的空域相位误差
phaseError3=extractValidPhaseErrorWithBounds(wrappedPhaseFractionalModulatedGroup-wrappedPhaseAllModulatedGroup,upPhaseErrorBound,bottomPhaseErrorBound);
plot(phaseError3,plotLineType,'Color','m','LineWidth',1.0);
title('Space Phase Error (Amplitude and Frequency Modulated Group)');
legend( ...
    sprintf('Space Phase Error (%d Steps Unit)',moveNumPart), ...
    sprintf('Space Phase Error (%d Steps Amplitude Modulated)',moveNumPart), ...
    sprintf('Space Phase Error (%d Steps Amplitude Modulated)',moveNumPart));
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

return;

grayList=ReadImages(gray_pre,grayOrder+2,moveNumAll);% grayOrder加上2以多读最后两张黑白图片
fringeListAll=ReadImages(pre,moveNumAll,0);


% ----------------------------------------

Pixelwidth = 48;
Fractional = 24;
ImageStart = 1;
pre = 'E:\XUWENYU\3D Reconstruction\Phase Error Compensation\Code\fringeImages\6\Left__';

FringeList24 = ReadImages(pre,Fractional,ImageStart);

% FringeList24 = CutFringe(FringeList24,1272);
FringeList24  = SubFringeListData(FringeList24,200,200,1000,1000);

Width = size(FringeList24{1},2);

Phase24 = GetWrapPhase(FringeList24,Fractional);

AB = GetABNormal(FringeList24,Phase24);

Fractional = 4;

FringeList6 = SelectNStepFring(FringeList24,Fractional);

Phase6 = GetWrapPhase(FringeList6,Fractional);


% 归一化
FringeListNormalize6 = NormalizeFringe(FringeList6,Fractional,AB);


HTFLNorm6 = HilbertPerRow(FringeListNormalize6,Fractional);
% 用截取后的条纹图求相位
Phase6HT = GetWrapPhaseWithHilbert(HTFLNorm6,Fractional);

% q求平均相位
MeanPhase6 = (Phase6 + Phase6HT)/2;

% %--------------------------------------------------------------- 下方是延拓
% 选择一个准确的初始相位
P = Phase6(400,Width);
%  BackCount = zeros(Fractional,1);
[ContiFLNorm6, BackCount]= ContinFringeListBack(FringeListNormalize6,Fractional,Pixelwidth,P,Phase6);

%  FrontCount = zeros(Fractional,1);
P = Phase6(400,1);
[ContiFLNorm6, FrontCount]= ContinFringeListFront(ContiFLNorm6,Fractional,Pixelwidth,P,Phase6);

% 对延拓后的条纹图进行HT
HTConFLNorm6 = HilbertPerRowWithLength(ContiFLNorm6,Fractional,Width,FrontCount,BackCount);

% 对HT后的条纹图截取到延拓前的宽度
HTConFLNormCut6 = CutFringeList(HTConFLNorm6,Width,FrontCount,BackCount);

% 用截取后的条纹图求相位
Phase6HT = GetWrapPhaseWithHilbert(HTConFLNormCut6,Fractional);

% q求平均相位
ConMeanPhase6 = (Phase6 + Phase6HT)/2;


start = 100;
over = 101;
number = 1;
for i = start : over
    figure('Name',sprintf('Line: %d',i));
%     myplot(extractValidPhaseErrorWithBounds(MeanPhase6-Phase24,     upBoundOfPhaseError,bottomBoundOfPhaseError),i,number);
    plot(extractValidPhaseErrorWithBounds(MeanPhase6(i,:)-Phase24(i,:),     upBoundOfPhaseError,bottomBoundOfPhaseError),'--.');
    hold on;
%     myplot(extractValidPhaseErrorWithBounds(ConMeanPhase6-Phase24,  upBoundOfPhaseError,bottomBoundOfPhaseError),i,number);
    plot(extractValidPhaseErrorWithBounds(ConMeanPhase6(i,:)-Phase24(i,:),  upBoundOfPhaseError,bottomBoundOfPhaseError),'LineWidth',1.2);
    grid on;
    xlim([-1,1001]);
    title(sprintf('%d-Steps of Mean Phase Error on Hilbert Domain',Fractional));
    legend('Non-Integate Periods','Integate Periods after Continuation')
    number=number+1;
end

% Spherical Gradient Illumination

%---------------------------------------------------------------------------------------------------------
% x:提取5个元素
% 从phase中从startIndex 开始取5个元素填充arr
% forward 表示是向前搜索还是向后搜索，默认向前
function [arr] = FillArrWith5Ele(phase,startIndex,arr,forward)
if nargin == 3
    forward = 1;
end
if forward == 1
    arr(1,:)=phase(1,startIndex-4:startIndex);
else
    arr(1,:)=phase(1,startIndex:startIndex+4);
end
end

% x:在折叠相位中逐位置移动
% 从phase中去index处的元素放入arr的第一个位置，其余元素依次后移一个位置
% % forward 表示是向前搜索还是向后搜索，默认向前
function [arr] = MoveOneStep(phase,startIndex,arr,forward)
if nargin == 3
    forward = 1;
end
pha = phase(1,startIndex);
% 如果碰到当前元素大于3，说明这时处于相位图的波峰和波谷交界处，此时应该越过
if pha > 3.0
    FillArrWith5Ele(phase,startIndex,arr,forward);
    return;
end
arr = [pha arr(1,1:4)];
end

% x:搜索最靠近边界的周期峰值位置
% 从相位图中找上一个周期中相位的最大值所在的位置
% Pos 是一个数组，每个元素记录一行中该周期里的最大相位的位置
function [Pos] = FindMaxPhasePos(Phase,startIndex,Forward)
if nargin == 2
    Forward = 1;
end;
height = size(Phase,1);
Pos = zeros(height,1);
Width = size(Phase,2);
if Forward == 1%向前查找
    for h = 1 : height
        for i = startIndex : -1 : 1
            pha0 = Phase(h,i);
            pha1 = Phase(h,i-1);
            if (pha1 - pha0) > 0 && abs(pha1 - pha0) > pi
                Pos(h,1) = i - 1;
                %                    return;
                break;
            end
        end
    end
else%向后查找
    for h = 1 : height
        for i = startIndex : Width
            pha0 = Phase(h,i);
            pha1 = Phase(h,i+1);
            if (pha1 - pha0) < 0 && abs(pha1 - pha0) > pi
                Pos(h,1) = i;
                %                        return;
                break;
            end
        end
    end
end
end

% x:--
% 在相位图的上一个周期中找与给定位置处的相位P对应的相位所在的位置索引
% startIndex是开始查找的位置
% Forward 是向前还是向后查找
function [Pos] = FindCorresPhasePos(Phase,startIndex,Forward)
if nargin == 2
    Forward = 1;
end;
height = size(Phase,1);
Pos = zeros(height,1);
Width = size(Phase,2);
if Forward == 1%向前查找
    for h = 1 : height
        if startIndex(h,1)==Width-2
            Pos(h,1) = Width-2;
            continue;
        end;
        for i = startIndex(h,1) : -1 : 1
            pha = Phase(h,i);
            if (pha - Phase(h,Width)) < 0
                Pos(h,1) = i + 1;
                %                             return;
                break;
            end
        end
    end
else%向后查找
    for h = 1 : height
        if startIndex(h,1)==2%如果起始索引查找索引是2，说明第一个相位处于峰值状态
            Pos(h,1) = 2;
            continue;
        end
        for i = startIndex(h,1) : Width
            pha = Phase(h,i);
            if (pha - Phase(h,1)) > 0
                Pos(h,1) = i - 1;
                %                             return;
                break;
            end
        end
    end
end
end


% ---------------------------------------------------------------------------------------------
% 工具函数
% ---------------------------------------------------------------------------------------------
% x:在折叠相位中找出相位对应的索引位置
function [Pos] = FindPosByPhase(Phase,startIndex,threshold,forward)
%
% 用位置 startIndex 处的相位值P 在 Phase 中的上一个周期中找与P值相等或相近的值所对应的像素坐标
if nargin == 3
    forward = 1;
end;
Pos = 0;
Pha = Phase(1,startIndex);
Num = size(Phase,1);
arr = zeros(1,5);
arr = FillArrWith5Ele(Phase,startIndex,arr,forward);
if forward == 1
    startIndex = startIndex - 5;
else
    startIndex = startIndex + 5;
end;

for i  = 1 : Num
    if abs(mean(arr)-Pha) < threshold
        if forward == 1
            Pos = startIndex + 3;
        else
            Pos = startIndex - 3;
        end;
        return ;
    end;
    display(abs(mean(arr)-Pha));
    arr = MoveOneStep(Phase,startIndex,arr);
    
    if forward == 1
        startIndex = startIndex - 1;
    else
        startIndex = startIndex + 1;
    end;
end;

end

% x:复制延拓，根据折叠相位的漏缺部分来对条纹信号进行延拓
function [res] = BlitData(image,src,length,dst,front)
% image 是输入图像
% src 是开始列
% length 是需要从src 开始包括src列，共复制的列数
% dst 是把length个列的数据拷贝到的目标索引
% 函数的作用是从图像的src列开始截取 length 个列添加到以dst为起始索引的后面
% front ： 1 表示在把数据复制到图像头部，0表示把数据复制到图像尾部
if nargin == 4
    front = 0;
end;
res = image;
height = size(image,1);
width = size(image,2);

if front == 0%把拷贝数据添加到图像尾部
    for h = 1 : height
        if src(h,1) > width || length(h,1) <= 0
            %                 display('超出图像宽度');
            continue;
        end;
        res(h,dst:dst+length(h,1)-1) = image(h, src(h,1):src(h,1) + length(h,1) - 1);
    end;
else%把拷贝数据添加到图像头部
    for h = 1 : height
        if src(h,1) > width || length(h,1) <= 0
            %                         display('超出图像宽度');
            continue;
        end;
        res(h,1:length(h,1)+width) = [image(h, src(h,1):src(h,1) + length(h,1) - 1) image(h,:)];
    end;
end;
end

% 
function [images] = CutFringeList(FringeList,Width,FrontCount,BackCount)
%     对 FringList 截取到指定的宽度 width
Num = size(FringeList,1);
RowNum = size(FringeList{1,1},1);
images = cell(Num,1);
for i = 1 : Num
    for r = 1 : RowNum
        front = FrontCount(r,1);
        if front<0
            front = 0;
        end;
        images{i}(r,:) = FringeList{i}(r,front+1:front + Width);
    end;
end;
end

% function [Pos] = MaxIntencityPos(image,startIndex)
%     Num = size(image,2);
%     Inten = image(1,startIndex);
%     for i = startIndex : Num
%         if image(1,startIndex)
%     end;
% end

function [ContiFLNorm, BackCount] = ContinFringeListBack(FringeList,Fractional,Pixelwidth,P,Phase)
% 对条纹图的元胞数组进行延拓，向后方延拓
% Pixelwidth 是条纹周期的像素宽度
% P 是条纹数组中第一幅条纹图的最后一列的相位值
%     BackCount = zeros(Fractional,1);
ContiFLNorm = FringeList;
Width = size(FringeList{1},2);
%     src  = Width - Pixelwidth;
%     threshold = 0.15;
%     src  = FindPosByPhase(Phase,Width,threshold);
pos = FindMaxPhasePos(Phase,Width);%当前周期最大相位所在位置
src = FindCorresPhasePos(Phase,pos);
src = src + 1;
length = pos - src + 1;
BackCount = length;
%     P  = P + pi;
for index = 1 : Fractional
    %             Ph =P + (index-1) * 2*pi/Fractional;
    %             Ph = mod(Ph,2*pi);
    %             length = floor(Pixelwidth -  Ph * Pixelwidth / (2 * pi));
    %             if length<0
    %                 length = 0;
    %             end;
    ContiFLNorm{index} = BlitData(FringeList{index},src,length,Width+1,0);
    %             BackCount(index,1)=length;
end;
end

function [ContiFLNorm,FrontCount] = ContinFringeListFront(FringeList,Fractional,Pixelwidth,P,Phase)
% 对条纹图的元胞数组进行延拓，向前方延拓
% Pixelwidth 是条纹周期的像素宽度
% P 是条纹数组中第一幅条纹图的最后一列的相位值
%     FrontCount = zeros(Fractional,1);
ContiFLNorm = FringeList;
%     P = P + pi;

src =1 + FindMaxPhasePos(Phase,1,0);%上一个周期最小相位所在位置
pos = FindCorresPhasePos(Phase,src,0);
pos = pos -1;
length = pos - src + 1;
FrontCount = length;
for index = 1 : Fractional
    %             Ph =P + (index-1) * 2*pi/Fractional;
    %             Ph = mod(Ph,2*pi);
    %             length = floor(Ph * Pixelwidth / (2 * pi));
    %             if length<0
    %                 length = 0;
    %             end;
    %             src  = Pixelwidth - length;
    ContiFLNorm{index} = BlitData(FringeList{index},src,length,1,1);
    %             FrontCount(index,1) = length;
end;
end


function [FL6Cut] = CutFringe(FringeList,Width)
Num = size(FringeList,1);
FL6Cut = cell(Num,1);
for i = 1 : Num
    FL6Cut{i} = FringeList{i}(:,1:Width);
end;
end


function [FringeListHT] = HilbertPerRowWithLength(FringeList,Fractional,Width,FrontLen,BackLen)

FringeListHT = cell(Fractional,1);
RowNum = size(FringeList{1,1},1);
for k = 1 : Fractional
    for r=1:RowNum
        front = FrontLen(r,1);
        if front < 0
            front =0;
        end;
        rr = FringeList{k}(r,1:Width + BackLen(r,1) + front);
        temp = hilbert(rr);
        FringeListHT{k}(r,1:Width + BackLen(r,1) + front) = -imag(temp);
    end
end
end