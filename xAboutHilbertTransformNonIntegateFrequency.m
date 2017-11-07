% xAboutHilbertTransformNonIntegateFrequency.m
% 徐文宇，20171106，周一
% 分析非整数周期采样对Hilbert变换以及Hilbert域相位的影响
% ver：---
close all;clear;

%% 设置基本参数***************************************************
% 图像条纹参数
width=1024; height=800; period=64;

%% -频率调制
% 整数周期标记
integatePeriodsFlag=1;
% 非整数周期标记
nonIntegatePeriodsFlag=1;

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
yTick=zeros(1,16+1);
yTickLabel=cell(1,16+1);
yTickLabel{9}='0';
for xt=1:8
    yTick(9+xt)=floor(xt*32);  yTickLabel{9+xt}=num2str(yTick(9+xt)); 
    yTick(9-xt)=floor(-xt*32); yTickLabel{9-xt}=num2str(yTick(9-xt));
end
% 相位误差显示有效区间
upPhaseErrorBound=2; bottomPhaseErrorBound=-2;

% plot画线类型
plotLineType=''; % ''实线，':'点虚线

%% {整数周期频率}******************************************************
% 整数周期标记
if integatePeriodsFlag==1
grayOrderStepFrequency =ceil(log(double(width/period))/log(2.0))+1;
preStepFrequency='E:\XUWENYU\3D Reconstruction\Phase Error Compensation\Code\fringeImages\39SameFrequency(T64)\I_';
gray_preStepFrequency='E:\XUWENYU\3D Reconstruction\Phase Error Compensation\Code\grayImages\28\I_';

% 读取全部条纹图像
moveNumAll=24;
grayListStepFrequency=ReadImages(gray_preStepFrequency,grayOrderStepFrequency+2,moveNumAll);% grayOrder加上2以多读最后两张黑白图片
fringeListAllStepFrequency=ReadImages(preStepFrequency,moveNumAll,0);
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    fringeListAllStepFrequency{k}=0.5*fringeListAllStepFrequency{k}(400,startOfSignal:endOfSignal);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListMoveNumStepFrequency=SelectNStepFring(fringeListAllStepFrequency,moveNumPart);

%% -计算理想空域相位
wrappedPhaseAllStepFrequency=GetWrapPhase(fringeListAllStepFrequency,moveNumAll);
AB=GetABNormal(fringeListAllStepFrequency,wrappedPhaseAllStepFrequency);

%% -计算数步相移条纹的空域相位
wrappedPhaseMoveNumStepFrequency=GetWrapPhase(fringeListMoveNumStepFrequency,moveNumPart);

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListMoveNumHilbertStepFrequency=HilbertPerRow(fringeListMoveNumStepFrequency,moveNumPart);
wrappedPhaseMoveNumHilbertStepFrequency=GetWrapPhaseWithHilbert(fringeListMoveNumHilbertStepFrequency,moveNumPart);

%% -显示整数周期数步相移条纹信号及其Hilbert变换、相位误差图表
% 显示信号及其Hilbert变换
% 第1/4步相移条纹
figure('name','1.Integate Periods of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepFrequency{1},                       plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepFrequency{1})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('Integate Periods of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96 192]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第2/4步相移条纹
figure('name','2.Integate Periods of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepFrequency{2},                       plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepFrequency{2})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('2/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96 192]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第3/4步相移条纹
figure('name','3.Integate Periods of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepFrequency{3},                       plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepFrequency{3})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('3/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96 192]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第4/4步相移条纹
figure('name','4.Integate Periods of Original fringe','NumberTitle','off');
plot(fringeListMoveNumStepFrequency{4},                       plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumStepFrequency{4})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('4/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96 192]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Integate Periods)','NumberTitle','off');
% 空域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumStepFrequency       -wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepFrequency-wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'Color','g','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',4);hold on;
title('Phase Error (Frequency Changed by Step Function)');
legend('Space Phase Error','HT Phase Error','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域相位误差的平均值与最大值
fprintf('------------integatePeriods-------------\n');
fprintf('        Mean of Space Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumStepFrequency-wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('Max positive of Space Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumStepFrequency-wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max negative of Space Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumStepFrequency-wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('           Mean of HT Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepFrequency-wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('   Max positive of HT Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepFrequency-wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('   Max negetive of HT Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepFrequency-wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound))));

end


%% {非整数周期}******************************************************
% 非整数周期标记
if nonIntegatePeriodsFlag==1
grayOrderSymmetricalArcFrequency =ceil(log(double(width/period))/log(2.0))+1;
preSymmetricalArcFrequency='E:\XUWENYU\3D Reconstruction\Phase Error Compensation\Code\fringeImages\39SameFrequency(T64)\I_';
gray_preSymmetricalArcFrequency='E:\XUWENYU\3D Reconstruction\Phase Error Compensation\Code\grayImages\28\I_';

% 读取全部条纹图像
moveNumAll=24;
grayListSymmetricalArcFrequency=ReadImages(gray_preSymmetricalArcFrequency,grayOrderSymmetricalArcFrequency+2,moveNumAll);% grayOrder加上2以多读最后两张黑白图片
fringeListAllSymmetricalArcFrequency=ReadImages(preSymmetricalArcFrequency,moveNumAll,0);
% 提取指定行列的信号，并直接设置为半幅度
for k=1:moveNumAll
    fringeListAllSymmetricalArcFrequency{k}=0.5*fringeListAllSymmetricalArcFrequency{k}(400,startOfSignal:endOfSignal-period/2);
end

% 抽取出数步相移条纹图像
moveNumPart=4;
fringeListMoveNumSymmetricalArcFrequency=SelectNStepFring(fringeListAllSymmetricalArcFrequency,moveNumPart);

%% -计算理想空域相位
wrappedPhaseAllSymmetricalArcFrequency=GetWrapPhase(fringeListAllSymmetricalArcFrequency,moveNumAll);
AB=GetABNormal(fringeListAllSymmetricalArcFrequency,wrappedPhaseAllSymmetricalArcFrequency);

%% -计算数步相移条纹的空域相位
wrappedPhaseMoveNumSymmetricalArcFrequency=GetWrapPhase(fringeListMoveNumSymmetricalArcFrequency,moveNumPart);

%% -计算数步相移条纹的Hilbert变换与Hilbert域相位
fringeListMoveNumHilbertSymmetricalArcFrequency=HilbertPerRow(fringeListMoveNumSymmetricalArcFrequency,moveNumPart);
wrappedPhaseMoveNumHilbertSymmetricalArcFrequency=GetWrapPhaseWithHilbert(fringeListMoveNumHilbertSymmetricalArcFrequency,moveNumPart);

%% -显示非整数周期数步相移条纹信号及其Hilbert变换、相位误差图表

% 显示信号及其Hilbert变换
% 第1/4步相移条纹
figure('name','1.Non-Integate Periods of Original fringe','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcFrequency{1},             plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcFrequency{1})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('1/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96 192]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第2/4步相移条纹
figure('name','2.Non-Integate Periods of Original fringe','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcFrequency{2},             plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcFrequency{2})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('2/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96 192]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第3/4步相移条纹
figure('name','3.Non-Integate Periods of Original fringe','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcFrequency{3},             plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcFrequency{3})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('3/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96 192]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);
% 第4/4步相移条纹
figure('name','4.Non-Integate Periods of Original fringe','NumberTitle','off');
plot(fringeListMoveNumSymmetricalArcFrequency{4},             plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
plot(imag(hilbert(fringeListMoveNumSymmetricalArcFrequency{4})),...
    plotLineType,'Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',6);hold on;
title('4/4 Step of Original fringe and its Hilbert Transform');
legend('Original fringe','HT','Location','NorthEast');
xlim([0,lengthOfSignal-1]);ylim([-96 192]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);
set(gca, 'YTick', yTick);set(gca, 'YTickLabel',yTickLabel);

% 显示空域相位误差、Hilbert域相位误差
figure('name','Phase Error (Non-Integate Periods)','NumberTitle','off');
% 空域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumSymmetricalArcFrequency       -wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'LineWidth',0.5,'MarkerSize',4);hold on;
% Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertStepFrequency          -wrappedPhaseAllStepFrequency,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'Color','g','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',4);hold on;
% non integate Hilbert域相位误差
plot(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcFrequency-wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound),...
    plotLineType,'Color','m','MarkerEdgeColor',[0.87,0.49,0],'LineWidth',0.5,'MarkerSize',4);hold on;

title('Phase Error (Non-Integate Periods)');
legend('Space Phase Error','HT Phase Error(Integate Periods)','HT Phase Error(Non-Integate Periods)','Location','SouthWest');
xlim([0,lengthOfSignal-1]);grid on;
set(gca, 'XTick', xTick);set(gca, 'XTickLabel',xTickLabel);

% 在命令行中显示空域/Hilbert域相位误差的平均值与最大值
fprintf('------------nonIntegatePeriods-------------\n');
fprintf('        Mean of Space Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumSymmetricalArcFrequency-wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('Max positive of Space Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumSymmetricalArcFrequency-wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('Max negative of Space Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumSymmetricalArcFrequency-wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('           Mean of HT Phase Error: %f\n',mean(extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcFrequency-wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound)));
fprintf('   Max positive of HT Phase Error: %f\n',max((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcFrequency-wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound))));
fprintf('   Max negetive of HT Phase Error: %f\n',min((extractValidPhaseErrorWithBounds(wrappedPhaseMoveNumHilbertSymmetricalArcFrequency-wrappedPhaseAllSymmetricalArcFrequency,upPhaseErrorBound,bottomPhaseErrorBound))));

end

