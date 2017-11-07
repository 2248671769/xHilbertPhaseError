% xAboutHilbertTransformDifferentAmplitude.m
% 徐文宇，20171102，周四
% 改变条纹信号的调制度，分析整数周期采样与非整数周期采样对Hilbert变换以及Hilbert域相位的影响
% ver：---
close all;clear;

%% 设置基本参数***************************************************
% 图像条纹参数
width=1280; height=800; period=128;

%% -幅度调制
% 阶跃调制幅度标记
stepFunctionModulateFlag=1;
% 对称连续弧段调制幅度标记
symmetricalArcModulateFlag=1;
% 非对称连续弧段调制幅度标记
asymmetricalArcModulateFlag=0;

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
yTick=zeros(1,16+1);
yTickLabel=cell(1,16+1);
yTickLabel{9}='0';
for xt=1:8
    yTick(9+xt)=floor(xt*32);  yTickLabel{9+xt}=num2str(yTick(9+xt)); 
    yTick(9-xt)=floor(-xt*32); yTickLabel{9-xt}=num2str(yTick(9-xt));
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
plotLineType=''; % ''实线，':.'点虚线

%% {单一幅度}******************************************************
%% -读取24幅全部条纹图像并从中抽取出数步相移条纹图像
grayOrder =ceil(log(double(width/period))/log(2.0))+1;
pre='E:\XUWENYU\3D Reconstruction\Phase Error Compensation\Code\fringeImages\41SameAmplitude\I_';
gray_pre='E:\XUWENYU\3D Reconstruction\Phase Error Compensation\Code\grayImages\28\I_';

% 读取全部条纹图像
moveNumAll=24;
grayList=ReadImages(gray_pre,grayOrder+2,moveNumAll);% grayOrder加上2以多读最后两张黑白图片
fringeListAll=ReadImages(pre,moveNumAll,0);
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    fringeListAll{k}=0.5*fringeListAll{k}(400,startOfSignal:endOfSignal);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListMoveNum=SelectNStepFring(fringeListAll,moveNumPart);

%% -计算理想空域相位
wrappedPhaseAll=GetWrapPhase(fringeListAll,moveNumAll);
AB=GetABNormal(fringeListAll,wrappedPhaseAll);

%% -计算数步相移条纹的空域相位
wrappedPhaseMoveNum=GetWrapPhase(fringeListMoveNum,moveNumPart);

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListMoveNumHilbert=HilbertPerRow(fringeListMoveNum,moveNumPart);
wrappedPhaseMoveNumHilbert=GetWrapPhaseWithHilbert(fringeListMoveNumHilbert,moveNumPart);

%% -显示图表
% 显示信号及其Hilbert变换
figure('name','Original fringe','NumberTitle','off');
plot(fringeListMoveNum{2});hold on;
plot(imag(hilbert(fringeListMoveNum{1})),plotLineType,'MarkerSize',8);hold on;
title('Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);


%% {左半右全阶跃调制幅度}****************************************
if stepFunctionModulateFlag==1
%% -调整所有信号幅度为左全右半阶跃状态
% 阶跃函数
filterStepAmplitude=[ones(1,lengthOfSignal/2) 2*ones(1,lengthOfSignal/2)];
fringeListAllStepAmplitude=cell(size(fringeListAll));
fringeListMoveNumStepAmplitude=cell(size(fringeListMoveNum));
for k=1:moveNumAll
    fringeListAllStepAmplitude{k}=floor((fringeListAll{k}-255.0/4).*filterStepAmplitude+255.0/2);
end
for k=1:moveNumPart
    fringeListMoveNumStepAmplitude{k}=floor((fringeListMoveNum{k}-255.0/4).*filterStepAmplitude+255.0/2);
end

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListMoveNumHilbertStepAmplitude=HilbertPerRow(fringeListMoveNumStepAmplitude,moveNumPart);
wrappedPhaseMoveNumHilbertStepAmplitude=GetWrapPhaseWithHilbert(fringeListMoveNumHilbertStepAmplitude,moveNumPart);

%% -显示阶跃函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示阶跃函数及其Hilbert变换
htStepAmplitudeFilter=imag(hilbert(filterStepAmplitude));
figure('name','Step function and its Hilbert Transform','NumberTitle','off');
plot(filterStepAmplitude,  plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(htStepAmplitudeFilter,plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
title('Step function and its Hilbert Transform');
legend('Step function','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
% 显示阶跃函数的Fourier频谱
figure('name','Fourier Spectrum of Step function','NumberTitle','off');
plot(abs(fftshift(fft(filterStepAmplitude))), '','LineWidth',0.5,'MarkerSize',4);hold on;
title('Fourier Spectrum of Step function');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);

% 显示信号及其Hilbert变换
% 第1/4步相移条纹
figure('name','1.Half Amplitude of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepAmplitude{1},                       plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepAmplitude{1})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('1/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第2/4步相移条纹
figure('name','2.Half Amplitude of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepAmplitude{2},                       plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepAmplitude{2})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('2/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第3/4步相移条纹
figure('name','3.Half Amplitude of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepAmplitude{3},                       plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepAmplitude{3})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('3/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第4/4步相移条纹
figure('name','4.Half Amplitude of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepAmplitude{4},                       plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepAmplitude{4})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('4/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Changed by Step Function)','NumberTitle','off');
% 空域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                    -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'Color','m','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',4);hold on;
title('Phase Error (Amplitude Changed by Step Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/阶跃式Hilbert域相位误差的平均值与最大值
fprintf('------------stepFunctionModulate-------------\n');
fprintf('        Mean of Space Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('Max positive of Space Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max negative of Space Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('           Mean of HT Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('   Max positive of HT Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('   Max negetive of HT Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));

end

%% {对称连续弧段调制幅度}****************************************
if symmetricalArcModulateFlag==1
%% -以对称连续弧段函数调制所有信号
% 对称连续弧段函数
arcCenter=0.5;
symmetricalArc=1-2.5*((0:lengthOfSignal-1)/lengthOfSignal-arcCenter).^2;
fringeListAllSymmetricalArcAmplitude=cell(size(fringeListAll));
fringeListMoveNumSymmetricalArcAmplitude=cell(size(fringeListMoveNum));
for k=1:moveNumAll
    fringeListAllSymmetricalArcAmplitude{k}=floor((fringeListAll{k}-255.0/4).*symmetricalArc+255.0/2);
end
for k=1:moveNumPart
    fringeListMoveNumSymmetricalArcAmplitude{k}=floor(fringeListMoveNum{k}.*symmetricalArc);
end

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListMoveNumHilbertSymmetricalArcAmplitude=HilbertPerRow(fringeListMoveNumSymmetricalArcAmplitude,moveNumPart);
wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude=GetWrapPhaseWithHilbert(fringeListMoveNumHilbertSymmetricalArcAmplitude,moveNumPart);

%% -显示对称连续弧段函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 组合显示原始信号、对称连续弧段函数、调制后的组合函数
figure('name','Symmetrical Arc Function','NumberTitle','off');
plot(symmetricalArc, plotLineType,'LineWidth',0.5,'MarkerSize',4);
title('Symmetrical Arc Function');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 组合显示原始信号、调制后的组合函数
figure('name','1/4 Step of Original fringe & Modulated signal by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNum{1},                       plotLineType,'LineWidth',1.5,'MarkerSize',4);hold on;
plot(fringeListMoveNumSymmetricalArcAmplitude{1},plotLineType,'LineWidth',1.5,'MarkerSize',4);
title('1/4 Step of Original fringe & Modulated signal by Symmetrical Arc Function');
legend('Fringe signal','Modulated signal','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 分别显示原始信号、对称连续弧段函数、调制后的组合函数的Fourier变换
% for k=1:moveNumPart 
fftFringeListMoveNum=fftshift(fft(fringeListMoveNum{1}));
fftFringeListMoveNumSymmetricalArcAmplitude=fftshift(fft(fringeListMoveNumSymmetricalArcAmplitude{1}));
fftSymmetricalArcAmplitudeFilter=fftshift(fft(symmetricalArc));
figure('name','(real part) Fourier Transform of Fringe signal, Modulated signal, Symmetrical Arc Function','NumberTitle','off');
subplot(3,1,1)
plot(abs(fftshift(fft(fringeListMoveNum{1}))),                       plotLineType,'LineWidth',1,'MarkerSize',4);
title('Fourier Spectrum of Fringe signal');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
subplot(3,1,2)
plot(abs(fftshift(fft(fringeListMoveNumSymmetricalArcAmplitude{1}))),plotLineType,'LineWidth',1,'MarkerSize',4);
title('Fourier Spectrum of Modulated signal');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
subplot(3,1,3)
plot(abs(fftshift(fft(symmetricalArc))),                             plotLineType,'LineWidth',1,'MarkerSize',4);
title('Fourier Spectrum of Symmetrical Arc Function');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);

% 显示对称连续弧段函数及其Hilbert变换
htSymmetricalArcAmplitudeFilter=imag(hilbert(symmetricalArc));
figure('name','Symmetrical Arc Function and its Hilbert Transform','NumberTitle','off');
plot(symmetricalArc,                 plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(htSymmetricalArcAmplitudeFilter,plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
title('Symmetrical Arc Function and its Hilbert Transform');
legend('Symmetrical Arc Function ','HT','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
% 第1/4步相移条纹
figure('name','1.Amplitude Changed by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcAmplitude{1},             plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcAmplitude{1})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('1/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第2/4步相移条纹
figure('name','2.Amplitude Changed by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcAmplitude{2},             plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcAmplitude{2})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('2/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第3/4步相移条纹
figure('name','3.Amplitude Changed by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcAmplitude{3},             plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcAmplitude{3})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('3/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第4/4步相移条纹
figure('name','4.Amplitude Changed by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcAmplitude{4},             plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcAmplitude{4})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('4/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Changed by Symmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                              -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'Color','m','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1,'MarkerSize',4);hold on;
title('Phase Error (Amplitude Changed by Symmetrical Arc Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/阶跃式Hilbert域相位误差的平均值与最大值
fprintf('------------symmetricalArcModulate-------------\n');
fprintf('        Mean of Space Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('Max positive of Space Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max negative of Space Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('           Mean of HT Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('   Max positive of HT Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('   Max negetive of HT Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
end

%% {非对称连续弧段调制幅度}**************************************
if asymmetricalArcModulateFlag==1
%% -以非对称连续弧段函数调制所有信号
% 非对称连续弧段函数
arcCenter=0.625;
asymmetricalArc=1-1.0*((0:lengthOfSignal-1)/lengthOfSignal-arcCenter).^2;
fringeListAllAsymmetricalArcAmplitude=cell(size(fringeListAll));
fringeListMoveNumAsymmetricalArcAmplitude=cell(size(fringeListMoveNum));
for k=1:moveNumAll
    fringeListAllAsymmetricalArcAmplitude{k}=floor((fringeListAll{k}-255.0/4).*asymmetricalArc+255.0/2);
end
for k=1:moveNumPart
    fringeListMoveNumAsymmetricalArcAmplitude{k}=floor(fringeListMoveNum{k}.*asymmetricalArc);
end

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListMoveNumHilbertAsymmetricalArcAmplitude=HilbertPerRow(fringeListMoveNumAsymmetricalArcAmplitude,moveNumPart);
wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude=GetWrapPhaseWithHilbert(fringeListMoveNumHilbertAsymmetricalArcAmplitude,moveNumPart);

%% -显示非对称连续弧段函数、数步相移条纹信号及其Hilbert变换、相位误差图表
% 组合显示原始信号、对称非连续弧段函数、调制后的组合函数
figure('name','Asymmetrical Arc Function','NumberTitle','off');
plot(asymmetricalArc, plotLineType,'LineWidth',0.5,'MarkerSize',4);
title('Asymmetrical Arc Function');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 组合显示原始信号、调制后的组合函数
figure('name','1/4 Step of Original fringe & Modulated signal by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNum{1},                        plotLineType,'LineWidth',1.5,'MarkerSize',4);hold on;
plot(fringeListMoveNumAsymmetricalArcAmplitude{1},plotLineType,'LineWidth',1.5,'MarkerSize',4);
title('1/4 Step of Original fringe & Modulated signal by Asymmetrical Arc Function');
legend('Fringe signal','Modulated signal','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
 
% 分别显示原始信号、非对称连续弧段函数、调制后的组合函数的Fourier变换
% for k=1:moveNumPart 
fftFringeListMoveNum=fftshift(fft(fringeListMoveNum{1}));
fftFringeListMoveNumAsymmetricalArcAmplitude=fftshift(fft(fringeListMoveNumAsymmetricalArcAmplitude{1}));
fftAsymmetricalArcAmplitudeFilter=fftshift(fft(asymmetricalArc));
% (real part) Fourier Transform
figure('name','(real part) Fourier Transform of Fringe signal, Asymmetrical Arc Function, Modulated signal','NumberTitle','off');
subplot(2,1,1)
plot(real(fftFringeListMoveNum),                         plotLineType,'LineWidth',1,'MarkerSize',4);hold on;
plot(real(fftFringeListMoveNumAsymmetricalArcAmplitude), plotLineType,'LineWidth',1,'MarkerSize',4);
title('(real part) Fourier Transform of Fringe signal & Modulated signal');
legend('Fringe signal','Modulated signal','Location','NorthEast');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
subplot(2,1,2);
plot(real(fftAsymmetricalArcAmplitudeFilter),            plotLineType,'LineWidth',1,'MarkerSize',4);
title('(real part) Fourier Transform of Asymmetrical Arc Function');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);

% (imaginary part) Fourier Transform
figure('name','(imaginary part) Fourier Transform of Fringe signal, Asymmetrical Arc Function, Modulated signal','NumberTitle','off');
subplot(2,1,1);
plot(imag(fftFringeListMoveNum),                         plotLineType,'LineWidth',1,'MarkerSize',4);hold on;
plot(imag(fftFringeListMoveNumAsymmetricalArcAmplitude), plotLineType,'LineWidth',1,'MarkerSize',4);
title('(imaginary part) Fourier Transform of Fringe signal & Modulated signal');
legend('Fringe signal','Modulated signal','Location','SouthEast');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);
subplot(2,1,2);
plot(imag(fftAsymmetricalArcAmplitudeFilter),            plotLineType,'LineWidth',1,'MarkerSize',4);
title('(imaginary part) Fourier Transform of Asymmetrical Arc Function');
xlim([min(xTickFourierSpectrum),max(xTickFourierSpectrum)]);ylim([-1,1]);grid on;
set(gca, 'XTick', xTickFourierSpectrum);set(gca, 'XTickLabel',xTickLabelFourierSpectrum);

% 显示非对称连续弧段函数及其Hilbert变换
htAsymmetricalArcAmplitudeFilter=imag(hilbert(asymmetricalArc));
figure('name','Asymmetrical Arc Function and its Hilbert Transform','NumberTitle','off');
plot(asymmetricalArc,                 plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(htAsymmetricalArcAmplitudeFilter,plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
title('Asymmetrical Arc Function and its Hilbert Transform');
legend('Asymmetrical Arc Function ','HT','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
% 第1/4步相移条纹
figure('name','1.Amplitude Changed by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumAsymmetricalArcAmplitude{1},            plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumAsymmetricalArcAmplitude{1})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('1/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第2/4步相移条纹
figure('name','2.Amplitude Changed by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumAsymmetricalArcAmplitude{2},            plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumAsymmetricalArcAmplitude{2})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('2/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第3/4步相移条纹
figure('name','3.Amplitude Changed by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumAsymmetricalArcAmplitude{3},            plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumAsymmetricalArcAmplitude{3})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('3/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第4/4步相移条纹
figure('name','4.Amplitude Changed by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumAsymmetricalArcAmplitude{4},            plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumAsymmetricalArcAmplitude{4})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('4/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Changed by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                               -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
if symmetricalArcModulateFlag==1
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'Color','g','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1,'MarkerSize',4);hold on;   
end
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'Color','m','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1,'MarkerSize',4);hold on;
title('Phase Error (Amplitude Changed by Asymmetrical Arc Function)');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
if symmetricalArcModulateFlag==1
    legend('Space Phase Error','HT Phase Error (Symmetrical)','HT Phase Error (Asymmetrical)','Location','NorthEast');
else
    legend('Space Phase Error','HT Phase Error','Location','NorthEast');
end
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域/阶跃式Hilbert域相位误差的平均值与最大值
fprintf('------------asymmetricalArcModulate-------------\n');
fprintf('        Mean of Space Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('Max positive of Space Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max negative of Space Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('           Mean of HT Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('   Max positive of HT Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('   Max negetive of HT Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
end




