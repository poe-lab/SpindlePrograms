newFs = 30000;
fStart = 300;

fStop = 2000;
fRes = 1;

% window = newFs * winSize;
% noverlap = floor(winOverlap * newFs);
window = 10;
noverlap = 6;

frequencyRange= fStart:fRes:fStop;
eegData = wfsmean_final_normal(200,:);
[~, F, T, sumP] = spectrogram(eegData,window,noverlap,frequencyRange,newFs);

for i = 201:208
    eegData = wfsmean_final_normal(i,:);
    [~, F, T, P] = spectrogram(eegData,window,noverlap,frequencyRange,newFs);
    sumP = sumP + P;
end
avgP = sumP./9;
avgP_dB = 10*log10(avgP);

figure2 = figure('Color',[1 1 1]);
% Create axes
axes('Parent',figure2,...
'Position',[0.13 0.219298245614035 0.775 0.705701754385965],...
'FontSize',20,'FontName','Arial');
surf(tCentered, F,avgP_AcrossUnits,'edgecolor','none');
colormap jet
axis tight;
view(0,90);
xlabel('Time (Seconds)','FontSize',20,'FontName','Arial');
ylabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
title('Spike-Centered Spectrogram','FontSize',16,'FontName','Arial');
colorbar
c = colorbar;
c.Label.String = 'Power/Frequency (\muV^{2}/Hz)';