function ReviewSpindles_11162017
numSpindles = size(EschenkoSpindle.startIdx,1); %Determine number of detected spindles on the selected LFP/EEG channel
% Determine maximum y-axis values (mV) for the spindle filtered signals plus
% added 10%
maxY_HcSpindle = max(abs(HcSpindleSignal)) + 0.1*max(abs(HcSpindleSignal));
maxY_CtxSpindle = max(abs(smooth(ctxSpindleSignal))) + 0.1*max(abs(smooth(ctxSpindleSignal)));

%% Calculate and plot results based on each detected spindle centered in a 3 sec window:
for i = 1:numSpindles     
    % Determine how many datapoints to add to each end of the spindle to
    % make a 3 sec window:
    spindleDataPts =  EschenkoSpindle.stopIdx(i) - EschenkoSpindle.startIdx(i) + 1;
	addPts = 300 - spindleDataPts;
    if addPts < 0   % If the spindle is the max allowed of 3 sec duration,
        addPts = 0; %do not add any data points to the beginning and end.
    end
    start = EschenkoSpindle.startIdx(i) - floor(addPts/2);
    stop = EschenkoSpindle.stopIdx(i) + ceil(addPts/2);

    cells = size(spikeData,2);
    spikesInSpindle = [];
    
    for c = 1:cells                                                         % For each place cell
        spikes = spikeData{c}/1000000;                                      % Keep all its spikes
        spikes = spikes(spikes >= TimeStamps(start) & spikes <= TimeStamps(stop)); % Keep all its spikes that fall within the target interval
        if ~isempty(spikes)                                                 % If there are any spikes
            spikesInSpindle = [spikesInSpindle; [spikes c*ones(size(spikes,1),1)]]; 
        end
    end

    if ~isempty(spikesInSpindle)
        lcSpikes = lcSpikeData/1000000;
        lcSpikes = lcSpikes(lcSpikes >= TimeStamps(start) & lcSpikes <= TimeStamps(stop));
        lcSpikesInSpindle = [];
        if ~isempty(lcSpikes)
            lcSpikesInSpindle = [lcSpikesInSpindle; [lcSpikes (cells+1)*ones(size(lcSpikes,1),1)]];
        end
        targetIdx = rippleFiltTS >= TimeStamps(start) & rippleFiltTS <= TimeStamps(stop);
        
        ax1 = subplot(3,1,1);
        plot(ax1, HcTimeStamps(start:stop), HcSpindleSignal(start:stop), 'b')
        hold on
        plot(ax1, HcTimeStamps(EschenkoSpindle.startIdx(i):EschenkoSpindle.stopIdx(i)),...
            HcSpindleSignal(EschenkoSpindle.startIdx(i):EschenkoSpindle.stopIdx(i)), 'g')
        plot(ax1, TimeStamps(start:stop), HcSignal(start:stop), 'Color', [0.5 0.5 0.5])
        plot(ax1, rippleFiltTS(targetIdx), rippleSignal(targetIdx), 'r')
        hold off
        switch EschenkoSpindle.scoring(i)
            case 2
                stage = 'SWS';
            case 3
                stage = 'REM';
            case 4
                stage = 'QW';
            case 6
                stage = 'TR ';
        end
        title(['Spindle ' num2str(i) ' of ' num2str(numSpindles) '          Stage: ' stage])
        xlim([TimeStamps(start) TimeStamps(stop)])
        ylim([-maxY_HcSpindle maxY_HcSpindle])
        ylabel(ax1, 'HC LFP')
        
        ax2 = subplot(3,1,2);
        map=lines(7);
        colorUnit = map(mod(spikesInSpindle(1,2),7)+1,:);
        plot(ax2, [spikesInSpindle(1,1) spikesInSpindle(1,1)], [spikesInSpindle(1,2)-.5 spikesInSpindle(1,2)+.5], 'Color', colorUnit)
        hold on
        for m = 2:size(spikesInSpindle,1)
            colorUnit = map(mod(spikesInSpindle(m,2),7)+1,:);             %de2bi(mod(spikesInSpindle(m,2),7),3);
            plot(ax2, [spikesInSpindle(m,1) spikesInSpindle(m,1)], [spikesInSpindle(m,2)-.5 spikesInSpindle(m,2)+.5], 'Color', colorUnit)
        end
%         scatter(ax2,spikesInSpindle(:,1),spikesInSpindle(:,2), 'b', 'd')
%         hold on
        if ~isempty(lcSpikesInSpindle)
            scatter(ax2,lcSpikesInSpindle(:,1),lcSpikesInSpindle(:,2), 'r', 's')
        end
        hold off
        ylim([0 cells+1])
        xlim([TimeStamps(start) TimeStamps(stop)])
        ylabel(ax2, 'Place Cells')
%         ax2.TickDirMode = 'manual';
%         ax2.TickDir = 'out';
        
        ax3 = subplot(3,1,3);
        plot(ax3, TimeStamps(start:stop), ctxSpindleSignal(start:stop), 'b')
        hold on
        plot(ax3, TimeStamps(EschenkoSpindle.startIdx(i):EschenkoSpindle.stopIdx(i)),...
            ctxSpindleSignal(EschenkoSpindle.startIdx(i):EschenkoSpindle.stopIdx(i)), 'g')
        plot(ax3, TimeStamps(start:stop), ctxSignal(start:stop), 'Color', [0.5 0.5 0.5])
        hold off
        xlim([TimeStamps(start) TimeStamps(stop)])
        ylim([-maxY_CtxSpindle maxY_CtxSpindle]);
        ylabel(ax3, 'Ctx EEG')
        xlabel(ax3, 'Time(seconds)');
        pause
    end

end