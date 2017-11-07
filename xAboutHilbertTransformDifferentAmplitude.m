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
asymmetricalArcModulateFlag=1;

% 信号范围
numOfPeriods=8;
startOfSignal=1;
endOfSignal=startOfSignal+period*numOfPeriods-1;
lengthOfSignal=endOfSignal-startOfSignal+1;
% xTick & xTickLabel
yTick=zeros(1,numOfPeriods+1);
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
% 相位误差显示有效区间
upPhaseErrorBound=2; bottomPhaseErrorBound=-2;

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
plot(imag(hilbert(fringeListMoveNum{1})),':.','MarkerSize',8);hold on;
title('Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
grid on;


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
plot(filterStepAmplitude,  ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(htStepAmplitudeFilter,':.','LineWidth',0.5,'MarkerSize',4);hold on;
title('Step function and its Hilbert Transform');
legend('Step function','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
% 第1/4步相移条纹
figure('name','1.Half Amplitude of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepAmplitude{1},':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepAmplitude{1})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('1/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第2/4步相移条纹
figure('name','2.Half Amplitude of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepAmplitude{2},':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepAmplitude{2})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1,'MarkerSize',6);hold on;
title('2/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第3/4步相移条纹
figure('name','3.Half Amplitude of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepAmplitude{3},':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepAmplitude{3})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('3/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第4/4步相移条纹
figure('name','4.Half Amplitude of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepAmplitude{4},':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepAmplitude{4})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('4/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','SouthWest');
xlim([0,lengthOfSignal-1]);ylim([-160 256]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Changed by Step Function)','NumberTitle','off');
% 空域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                    -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    ':.','LineWidth',0.5,'MarkerSize',4);hold on;
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    ':.','Color','m','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',4);hold on;
title('Phase Error (Amplitude Changed by Step Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/HIlbert域/阶跃式Hilbert域相位误差的平均值与最大值
fprintf('Mean of Space Phase Error: %s\n',num2str(mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                   -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max positive of Space Phase Error: %s\n',num2str(max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                 -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Max negative of Space Phase Error: %s\n',num2str(min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                 -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Mean of HT Phase Error: %s\n',num2str(mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepAmplitude  -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max positive of HT Phase Error: %s\n',num2str(max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Max negetive of HT Phase Error: %s\n',num2str(min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));

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
plot(symmetricalArc, ':.','LineWidth',0.5,'MarkerSize',4);
title('Symmetrical Arc Function');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 组合显示原始信号、调制后的组合函数
figure('name','1/4 Step of Original fringe & Demodulated signal by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNum{1},                       ':.','LineWidth',1.5,'MarkerSize',4);hold on;
plot(fringeListMoveNumSymmetricalArcAmplitude{1},':.','LineWidth',1.5,'MarkerSize',4);
title('1/4 Step of Original fringe & Demodulated signal by Symmetrical Arc Function');
legend('Fringe signal','Demodulated signal','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% xTickPart & xTickLabelPart
partNum=30;
xTickPart=zeros(1,partNum+1);
xTickLabelPart=cell(1,partNum+1);
for xt=0:partNum/2
    xTickPart(partNum/2+xt+1)=lengthOfSignal/2+1+xt; xTickLabelPart{partNum/2+xt+1}=num2str( xt);
    xTickPart(partNum/2-xt+1)=lengthOfSignal/2+1-xt; xTickLabelPart{partNum/2-xt+1}=num2str(-xt);
end
xTickPart(partNum/2+1)=lengthOfSignal/2+1; xTickLabelPart{partNum/2+1}=num2str(0);

% 分别显示原始信号、对称连续弧段函数、调制后的组合函数的Fourier变换
% for k=1:moveNumPart 
fftFringeListMoveNum=fftshift(fft(fringeListMoveNum{1}));
fftFringeListMoveNumSymmetricalArcAmplitude=fftshift(fft(fringeListMoveNumSymmetricalArcAmplitude{1}));
fftSymmetricalArcAmplitudeFilter=fftshift(fft(symmetricalArc));
% (real part) Fourier Transform
figure('name','(real part) Fourier Transform of Fringe signal, Symmetrical Arc Function, Demodulated signal','NumberTitle','off');
subplot(2,1,1)
plot(real(fftFringeListMoveNum),                        ':.','LineWidth',1,'MarkerSize',4);hold on;
plot(real(fftFringeListMoveNumSymmetricalArcAmplitude), ':.','LineWidth',1,'MarkerSize',4);
title('(real part) Fourier Transform of Fringe signal & Demodulated signal');
legend('Fringe signal','Demodulated signal','Location','NorthEast');
xlim([min(xTickPart),max(xTickPart)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickPart);set(gca, 'XTickLabel',xTickLabelPart);
subplot(2,1,2);
plot(real(fftSymmetricalArcAmplitudeFilter),            ':.','LineWidth',1,'MarkerSize',4);
title('(real part) Fourier Transform of Symmetrical Arc Function');
xlim([min(xTickPart),max(xTickPart)]);grid on;
set(gca, 'XTick', xTickPart);set(gca, 'XTickLabel',xTickLabelPart);

% (imaginary part) Fourier Transform
figure('name','(imaginary part) Fourier Transform of Fringe signal, Symmetrical Arc Function, Demodulated signal','NumberTitle','off');
subplot(2,1,1);
plot(imag(fftFringeListMoveNum),                        ':.','LineWidth',1,'MarkerSize',4);hold on;
plot(imag(fftFringeListMoveNumSymmetricalArcAmplitude), ':.','LineWidth',1,'MarkerSize',4);
title('(imaginary part) Fourier Transform of Fringe signal & Demodulated signal');
legend('Fringe signal','Demodulated signal','Location','SouthEast');
xlim([min(xTickPart),max(xTickPart)]);grid on;
set(gca, 'XTick', xTickPart);set(gca, 'XTickLabel',xTickLabelPart);
subplot(2,1,2);
plot(imag(fftSymmetricalArcAmplitudeFilter),            ':.','LineWidth',1,'MarkerSize',4);
title('(imaginary part) Fourier Transform of Symmetrical Arc Function');
xlim([min(xTickPart),max(xTickPart)]);ylim([-1,1]);grid on;
set(gca, 'XTick', xTickPart);set(gca, 'XTickLabel',xTickLabelPart);

% 显示对称连续弧段函数及其Hilbert变换
htSymmetricalArcAmplitudeFilter=imag(hilbert(symmetricalArc));
figure('name','Symmetrical Arc Function and its Hilbert Transform','NumberTitle','off');
plot(symmetricalArc,                 ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(htSymmetricalArcAmplitudeFilter,':.','LineWidth',0.5,'MarkerSize',4);hold on;
title('Symmetrical Arc Function and its Hilbert Transform');
legend('Symmetrical Arc Function ','HT','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
% 第1/4步相移条纹
figure('name','1.Amplitude Changed by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcAmplitude{1},             ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcAmplitude{1})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('1/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第2/4步相移条纹
figure('name','2.Amplitude Changed by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcAmplitude{2},             ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcAmplitude{2})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('2/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第3/4步相移条纹
figure('name','3.Amplitude Changed by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcAmplitude{3},             ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcAmplitude{3})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('3/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第4/4步相移条纹
figure('name','4.Amplitude Changed by Symmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcAmplitude{4},             ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcAmplitude{4})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('4/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Changed by Symmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                              -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    ':.','LineWidth',0.5,'MarkerSize',4);hold on;
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    ':.','Color','m','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1,'MarkerSize',4);hold on;
title('Phase Error (Amplitude Changed by Symmetrical Arc Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/HIlbert域/阶跃式Hilbert域相位误差的平均值与最大值
fprintf('Mean of Space Phase Error: %s\n',num2str(mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                   -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max positive of Space Phase Error: %s\n',num2str(max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                 -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Max negative of Space Phase Error: %s\n',num2str(min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                 -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Mean of HT Phase Error: %s\n',num2str(mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude  -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max positive of HT Phase Error: %s\n',num2str(max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Max negetive of HT Phase Error: %s\n',num2str(min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
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
plot(asymmetricalArc, ':.','LineWidth',0.5,'MarkerSize',4);
title('Asymmetrical Arc Function');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 组合显示原始信号、调制后的组合函数
figure('name','1/4 Step of Original fringe & Demodulated signal by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNum{1},                        ':.','LineWidth',1.5,'MarkerSize',4);hold on;
plot(fringeListMoveNumAsymmetricalArcAmplitude{1},':.','LineWidth',1.5,'MarkerSize',4);
title('1/4 Step of Original fringe & Demodulated signal by Asymmetrical Arc Function');
legend('Fringe signal','Demodulated signal','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-32,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% xTickPart & xTickLabelPart
partNum=30;
xTickPart=zeros(1,partNum+1);
xTickLabelPart=cell(1,partNum+1);
for xt=0:partNum/2
    xTickPart(partNum/2+xt+1)=lengthOfSignal/2+1+xt; xTickLabelPart{partNum/2+xt+1}=num2str( xt);
    xTickPart(partNum/2-xt+1)=lengthOfSignal/2+1-xt; xTickLabelPart{partNum/2-xt+1}=num2str(-xt);
end
xTickPart(partNum/2+1)=lengthOfSignal/2+1; xTickLabelPart{partNum/2+1}=num2str(0);

% 分别显示原始信号、非对称连续弧段函数、调制后的组合函数的Fourier变换
% for k=1:moveNumPart 
fftFringeListMoveNum=fftshift(fft(fringeListMoveNum{1}));
fftFringeListMoveNumAsymmetricalArcAmplitude=fftshift(fft(fringeListMoveNumAsymmetricalArcAmplitude{1}));
fftAsymmetricalArcAmplitudeFilter=fftshift(fft(asymmetricalArc));
% (real part) Fourier Transform
figure('name','(real part) Fourier Transform of Fringe signal, Asymmetrical Arc Function, Demodulated signal','NumberTitle','off');
subplot(2,1,1)
plot(real(fftFringeListMoveNum),                         ':.','LineWidth',1,'MarkerSize',4);hold on;
plot(real(fftFringeListMoveNumAsymmetricalArcAmplitude), ':.','LineWidth',1,'MarkerSize',4);
title('(real part) Fourier Transform of Fringe signal & Demodulated signal');
legend('Fringe signal','Demodulated signal','Location','NorthEast');
xlim([min(xTickPart),max(xTickPart)]);ylim([-2*10000,8*10000]);grid on;
set(gca, 'XTick', xTickPart);set(gca, 'XTickLabel',xTickLabelPart);
subplot(2,1,2);
plot(real(fftAsymmetricalArcAmplitudeFilter),            ':.','LineWidth',1,'MarkerSize',4);
title('(real part) Fourier Transform of Asymmetrical Arc Function');
xlim([min(xTickPart),max(xTickPart)]);grid on;
set(gca, 'XTick', xTickPart);set(gca, 'XTickLabel',xTickLabelPart);

% (imaginary part) Fourier Transform
figure('name','(imaginary part) Fourier Transform of Fringe signal, Asymmetrical Arc Function, Demodulated signal','NumberTitle','off');
subplot(2,1,1);
plot(imag(fftFringeListMoveNum),                         ':.','LineWidth',1,'MarkerSize',4);hold on;
plot(imag(fftFringeListMoveNumAsymmetricalArcAmplitude), ':.','LineWidth',1,'MarkerSize',4);
title('(imaginary part) Fourier Transform of Fringe signal & Demodulated signal');
legend('Fringe signal','Demodulated signal','Location','SouthEast');
xlim([min(xTickPart),max(xTickPart)]);grid on;
set(gca, 'XTick', xTickPart);set(gca, 'XTickLabel',xTickLabelPart);
subplot(2,1,2);
plot(imag(fftAsymmetricalArcAmplitudeFilter),            ':.','LineWidth',1,'MarkerSize',4);
title('(imaginary part) Fourier Transform of Asymmetrical Arc Function');
xlim([min(xTickPart),max(xTickPart)]);ylim([-1,1]);grid on;
set(gca, 'XTick', xTickPart);set(gca, 'XTickLabel',xTickLabelPart);

% 显示非对称连续弧段函数及其Hilbert变换
htAsymmetricalArcAmplitudeFilter=imag(hilbert(asymmetricalArc));
figure('name','Asymmetrical Arc Function and its Hilbert Transform','NumberTitle','off');
plot(asymmetricalArc,                 ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(htAsymmetricalArcAmplitudeFilter,':.','LineWidth',0.5,'MarkerSize',4);hold on;
title('Asymmetrical Arc Function and its Hilbert Transform');
legend('Asymmetrical Arc Function ','HT','Location','SouthEast');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 显示信号及其Hilbert变换
% 第1/4步相移条纹
figure('name','1.Amplitude Changed by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumAsymmetricalArcAmplitude{1},            ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumAsymmetricalArcAmplitude{1})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('1/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第2/4步相移条纹
figure('name','2.Amplitude Changed by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumAsymmetricalArcAmplitude{2},            ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumAsymmetricalArcAmplitude{2})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('2/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第3/4步相移条纹
figure('name','3.Amplitude Changed by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumAsymmetricalArcAmplitude{3},            ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumAsymmetricalArcAmplitude{3})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('3/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第4/4步相移条纹
figure('name','4.Amplitude Changed by Asymmetrical Arc Function','NumberTitle','off');
plot(fringeListMoveNumAsymmetricalArcAmplitude{4},            ':.','LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumAsymmetricalArcAmplitude{4})),...
    ':.','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('4/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthWest');
xlim([0,lengthOfSignal-1]);ylim([-96,160]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Amplitude Changed by Asymmetrical Arc Function)','NumberTitle','off');
% 空域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                               -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    ':.','LineWidth',0.5,'MarkerSize',4);hold on;
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
if symmetricalArcModulateFlag==1
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcAmplitude -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    ':.','Color','g','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1,'MarkerSize',4);hold on;   
end
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound),...
    ':.','Color','m','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1,'MarkerSize',4);hold on;
title('Phase Error (Amplitude Changed by Asymmetrical Arc Function)');
% 如果对称连续弧段调制幅度标记为1，亦显示其Hilbert域相位误差
if symmetricalArcModulateFlag==1
    legend('Space Phase Error','HT Phase Error (Symmetrical)','HT Phase Error (Asymmetrical)','Location','NorthEast');
else
    legend('Space Phase Error','HT Phase Error','Location','NorthEast');
end
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/HIlbert域/阶跃式Hilbert域相位误差的平均值与最大值
fprintf('Mean of Space Phase Error: %s\n',num2str(mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                   -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max positive of Space Phase Error: %s\n',num2str(max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                 -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Max negative of Space Phase Error: %s\n',num2str(min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNum                 -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Mean of HT Phase Error: %s\n',num2str(mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude  -wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max positive of HT Phase Error: %s\n',num2str(max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
fprintf('Max negetive of HT Phase Error: %s\n',num2str(min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertAsymmetricalArcAmplitude-wrappedPhaseAll,upPhaseErrorBound,bottomPhaseErrorBound)))));
end




