clear;
clc;

load('Fig4_workspace.mat');



mat_BF_N= cell2mat(BF_N);
mat_BF_Imild= cell2mat(BF_Imild);
mat_BF_Isev= cell2mat(BF_Isev);


mat_Thr_N= cell2mat(Thr_N);
mat_Thr_Imild= cell2mat(Thr_Imild);
mat_Thr_Isev= cell2mat(Thr_Isev);


figure(1);
clf;

thresh_boundaries= min(mat_Thr_N)-.1:2:max(mat_Thr_Isev);

subplot(311);
histogram(mat_Thr_N, thresh_boundaries);

subplot(312);
histogram(mat_Thr_Imild, thresh_boundaries);

subplot(313);
histogram(mat_Thr_Isev, thresh_boundaries);

figure(2);
clf;

%% CF boundary plot
CFcenters= [.5 1 2 4 8];
prctPoint= [5 10 50 75];

sp_ax= nan(6, 1);
for cfVar=1:length(CFcenters)
    curCFcent= CFcenters(cfVar);
    curCFstart= curCFcent/sqrt(2);
    curCFend= curCFcent*sqrt(2);
    
    
    sp_ax(cfVar)= subplot(3,2,cfVar);
    hold on;
    
    title(sprintf('%.1f < BF < %.1f', curCFstart, curCFend));
    
    cdf1edf0= 0;
    if cdf1edf0
        hN= cdfplot(mat_Thr_N( (mat_BF_N>curCFstart) & (mat_BF_N<curCFend)));
        hN.LineWidth= 2;
        hImild= cdfplot(mat_Thr_Imild( (mat_BF_Imild>curCFstart) & (mat_BF_Imild<curCFend) ));
        hImild.LineWidth= 2;
        hIsev= cdfplot(mat_Thr_Isev( (mat_BF_Isev>curCFstart) & (mat_BF_Isev<curCFend) ));
        hIsev.LineWidth= 2;
        
    else
        ecdf(mat_Thr_N( (mat_BF_N>curCFstart) & (mat_BF_N<curCFend)),'alpha',0.01, 'bounds', 'on');
        ecdf(mat_Thr_Imild( (mat_BF_Imild>curCFstart) & (mat_BF_Imild<curCFend) ),'alpha',0.01, 'bounds', 'on');
        ecdf(mat_Thr_Isev( (mat_BF_Isev>curCFstart) & (mat_BF_Isev<curCFend) ),'alpha',0.01, 'bounds', 'on');
    end
    
    title(sprintf('CF center = %.1f kHz', curCFcent));
    
    for prcVar=1:length(prctPoint)
        curPrctVal= prctPoint(prcVar);
        nhVal= prctile(mat_Thr_N((mat_BF_N>curCFstart) & (mat_BF_N<curCFend)), curPrctVal);
        ImildVal= prctile(mat_Thr_Imild( (mat_BF_Imild>curCFstart) & (mat_BF_Imild<curCFend) ), curPrctVal);
        IsevVal= prctile(mat_Thr_Isev( (mat_BF_Isev>curCFstart) & (mat_BF_Isev<curCFend) ), curPrctVal);
        if curPrctVal~=10
            text(60, curPrctVal/100, sprintf('Shift@%.0f|Mild=%.0f,Sev=%.0f dB', curPrctVal, ImildVal-nhVal, IsevVal-nhVal));
        else
            text(60, (curPrctVal+5)/100, sprintf('Shift@%.0f|Mild=%.0f,Sev=%.0f dB', curPrctVal, ImildVal-nhVal, IsevVal-nhVal));
        end
    end
    ylabel('');
    xlabel('');
end

sp_ax(6)= subplot(3,2,6);
hold on;


if cdf1edf0
    hN= cdfplot(mat_Thr_N);
    hImild= cdfplot(mat_Thr_Imild);
    hIsev= cdfplot(mat_Thr_Isev);
    hN.LineWidth= 2.5;
    hImild.LineWidth= 2.5;
    hIsev.LineWidth= 2.5;
else
    ecdf(mat_Thr_N,'alpha',0.01, 'bounds', 'on');
    ecdf(mat_Thr_Imild,'alpha',0.01, 'bounds', 'on');
    ecdf(mat_Thr_Isev,'alpha',0.01, 'bounds', 'on');
end
for prcVar=1:length(prctPoint)
    curPrctVal= prctPoint(prcVar);
    nhVal= prctile(mat_Thr_N, curPrctVal);
    ImildVal= prctile(mat_Thr_Imild, curPrctVal);
    IsevVal= prctile(mat_Thr_Isev, curPrctVal);
    if curPrctVal~=10
        text(60, curPrctVal/100, sprintf('Shift@%.0f|Mild=%.0f,Sev=%.0f dB', curPrctVal, ImildVal-nhVal, IsevVal-nhVal));
    else
        text(60, (curPrctVal+5)/100, sprintf('Shift@%.0f|Mild=%.0f,Sev=%.0f dB', curPrctVal, ImildVal-nhVal, IsevVal-nhVal));
    end
    
end
title('All BFs');

axes(sp_ax(5));
ylabel('CDF', 'FontSize', 12);
xlabel('Thresh (dB)', 'FontSize', 12);

linkaxes(sp_ax, 'x')

Xcorner= .05;
Xwidth= .38;
Xshift= .12;
Ycorner= .065;
Ywidth= .27;
Yshift= .05;

% D
set(sp_ax(5),'Position', [Xcorner Ycorner Xwidth Ywidth])
drawnow

% C
set(sp_ax(3),'Position', [Xcorner Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% A
set(sp_ax(1),'Position', [Xcorner Ycorner+2.15*Ywidth+2*Yshift Xwidth .85*Ywidth])
drawnow

% F
set(sp_ax(6),'Position', [Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow

% E
set(sp_ax(4),'Position', [Xcorner+Xwidth+Xshift Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

% B
set(sp_ax(2),'Position', [Xcorner+Xwidth+Xshift Ycorner+2*Ywidth+2*Yshift Xwidth Ywidth])
drawnow

print -dpng Thresh_Shift_Cat