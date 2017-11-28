% xDiscreteHilbertFilter.m
% 徐文宇，20171106，周一
% 绘制离散Hilbert算子
% ver：---
close all,clear;

% period
period=4; periodWaveNum=12; N=period*periodWaveNum;

% xTick & xTickLabel
xTick=zeros(1,periodWaveNum+1);
xTickLabel=cell(1,periodWaveNum+1);
unit=period;
for xt=0:periodWaveNum-1
    xTick(xt+1)=floor(xt*unit); xTickLabel{xt+1}=num2str(xTick(xt+1));
end
xTick(end)=N-1; xTickLabel{end}=num2str(N-1);

%% @_@{Discrete Hilbert Tranform Filter (Even)}********************
htFilterEven=zeros(1,N);
sinPi2Even=zeros(1,N);
cotPiNEven=zeros(1,N);
for i=0:N-1
    if i==0
        htFilterEven(i+1)=0;
    else
        htFilterEven(i+1)=2/N*power(sin(pi*i/2),2)*cot(pi*i/N);
    end
    sinPi2Even(i+1)=power(sin(pi*i/2),2);
    cotPiNEven(i+1)=2/N*cot(pi*i/N);      
end

figure('name',sprintf('Discrete Hilbert Tranform Filter (even)'));
% sin
plot(0:N-1,sinPi2Even,'.','MarkerSize',10);hold on;
% cot
plot(0:N-1,cotPiNEven,':.','Color',[0.75,0,0.75],'MarkerEdgeColor',[0.75,0,0.75],'LineWidth',0.5);hold on;
% filter: sin*cot
plot(0:N-1,htFilterEven,':o','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1.414);
title('Discrete Hilbert Tranform Filter ($N$=48, even)','Interpreter','latex');
h=legend('$\sin^2 (\frac{\pi}{2}$)','$\frac{2}{N}\cot(\frac{\pi n}{N}$)','DHT Filter (Even)');
set(h,'Interpreter','latex','Location','southwest','FontSize',14);
xlim([0 N-1]);ylim([-1.2 1.2]);
set(gca, 'XTick', xTick); set(gca, 'XTickLabel',xTickLabel)
grid on;

%% @_@{Discrete Hilbert Tranform Filter (Odd)}*********************
N=N+1;
htFilterOdd =zeros(1,N);
for i=1:N+1
    if i<N && sin(pi*i/N)~=0
        htFilterOdd(i)=1/N*(1-cos(pi*i)/sin(pi*i/N));
    end
end
figure('name',sprintf('Discrete Hilbert Tranform Filter (odd)'));
plot(0:N-1,htFilterOdd,':o','Color',[0,0.8078,0.8196],'MarkerEdgeColor',[0.87,0.49,0],'LineWidth',1.414);
title('Discrete Hilbert Tranform Filter ($N$=49, odd)','Interpreter','latex');
xlim([0 N-1]);ylim([-.5 .5]);
xTick(end)=N-1; xTickLabel{end}=num2str(N-1);
set(gca, 'XTick', xTick); set(gca, 'XTickLabel',xTickLabel)
grid on;
