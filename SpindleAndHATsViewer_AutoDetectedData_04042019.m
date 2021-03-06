function SpindleAndHATsViewer_AutoDetectedData_04042019
%% Select .mat file generated by the auto-detection program:
[fileName, filePath] = uigetfile({'*.mat',...
        'Pick data file.'},'Select .MAT File');
dataFile = fullfile(filePath, fileName);
load(dataFile, '-mat')

%% Determine whether or not data includes startTime and endTime variables from spike and spindle analyses:
% If the variables do not exist, create from the time stamp vector:
% if exist('startTime', 'var') == 0 
%     startTime = TimeStamps(1) * 10^6;
%     endTime = TimeStamps(end)  * 10^6;
% end
% clear TimeStamps

%% Set length of viewing window in seconds:
viewTimeWindow = [];
while isempty(viewTimeWindow)
    prompt={'Viewing Window Setting:'};
    dlgTitle='Set length of window (sec)';
    lineNo=1;
    defaultans = {'10'};
    options.Resize = 'on';
    answer = inputdlg(prompt,dlgTitle,lineNo, defaultans, options);
    viewTimeWindow = str2double(answer{1,1});
    clear answer prompt dlgTitle lineNo
end

%% Calculate & plot results for each spindle centered in a defined time window:
numSpindles = size(spindle.startIdx,1); %Determine number of detected spindles on the selected LFP/EEG channel
% Determine maximum y-axis values (mV) for the spindle filtered signals plus
% added 10%
maxY_Signal = max(spindle.maxP2pAmp) + 0.1*max(spindle.maxP2pAmp);
maxY_RMS = max(spindle.maxRMS) + 0.1*max(spindle.maxRMS);

eventNumber = 1;
runProgram = 1;
while runProgram == 1
% for i = 1:numSpindles
    addTime = viewTimeWindow - (spindle.timeBnds(eventNumber,2) - spindle.timeBnds(eventNumber,1));
    if addTime < 0   % If the spindle duration >= viewTimeWindow,
        addTime = 0; %do not add any data points to the beginning and end.
    end
    % Define the time range of the viewing window for the spindle:
    windowStartTime = spindle.timeBnds(eventNumber,1) - addTime/2;
    windowStopTime = spindle.timeBnds(eventNumber,2) + addTime/2;
    
    ax1 = subplot(2,1,1);  % creates a figure with axes
    
    % Find indices of time points within the viewing window:
    windowTargetIdx = info.TimeStamps >= windowStartTime & info.TimeStamps <= windowStopTime;
    
    % Isolate windowed data:
    windowTS = info.TimeStamps(windowTargetIdx);
    windowSignal = info.Signal(windowTargetIdx);
    windowFiltSignal = spindle.Signal(windowTargetIdx);
    windowRMS = spindle.RMS(windowTargetIdx);
    clear windowTargetIdx
    
    % Plot the sigma filtered signal for the viewing widow:
    plot(ax1, windowTS, windowFiltSignal, 'b')
    hold on % keeps plot open for additional signals
    
    % Find all of the events in the time window:
    eventsInWindow = find(spindle.timeBnds(:,1) >= windowTS(1) & spindle.timeBnds(:,2) <= windowTS(end));
    for i =1:length(eventsInWindow)
        % Find indices of time points for the event within the viewing window:
        eventTargetIdx = windowTS >= spindle.timeBnds(eventsInWindow(i),1) & windowTS <= spindle.timeBnds(eventsInWindow(i),2);
        % Plot the event in green:
        plot(ax1, windowTS(eventTargetIdx), windowFiltSignal(eventTargetIdx), 'g')
    end
%     % Find indices of time points for the spindle within the viewing window:
%     eventTargetIdx = windowTS >= spindle.timeBnds(eventNumber,1) & windowTS <= spindle.timeBnds(eventNumber,2);
%     
%     % Plot the spindle in green:
%     plot(ax1, windowTS(eventTargetIdx), windowFiltSignal(eventTargetIdx), 'g')
    
    % Plot the LFP/EEG signal for the viewing widow in gray:
    plot(ax1, windowTS, windowSignal, 'Color', [0.5 0.5 0.5])

    hold off
    switch spindle.state(eventNumber)
        case 2
            stage = 'SWS';
        case 3
            stage = 'REM';
        case 4
            stage = 'QW ';
        case 5
            stage = 'UH ';
        case 6
            stage = 'TR ';
    end
    title(['Spindle ' num2str(eventNumber) ' of ' num2str(numSpindles) '          Stage: ' stage])
    xlim([windowTS(1) windowTS(end)])
    ylim([-maxY_Signal maxY_Signal])
    ylabel(ax1, 'LFP')
    
    % RMS plot:
    ax2 = subplot(2,1,2);
    bar(ax2, windowTS, windowRMS, 'FaceColor', [0.5 0.5 0.5])
    hold on
    plot(ax2, [windowTS(1) windowTS(end)], [spindle.upperThreshold spindle.upperThreshold], 'g')
    plot(ax2, [windowTS(1) windowTS(end)], [spindle.lowerThreshold spindle.lowerThreshold], 'r')
    hold off
    xlim([windowTS(1) windowTS(end)])
    ylim([0 maxY_RMS]);
    ylabel(ax2, 'RMS of Filtered Signal')
    xlabel(ax2, 'Time(seconds)');
    
    checker = 1;
    while checker == 1
        fig = gcf;
        waitforbuttonpress
        u = fig.CurrentCharacter;
        switch u
            case '4'
                if eventNumber > 1
                    eventNumber = eventNumber - 1;
                end
                checker = 0;
            case '6'
                if eventNumber < numSpindles
                    eventNumber = eventNumber + 1;
                end
                checker = 0;
            case 'q'     
                checker = 0;
                runProgram = 0;
            otherwise
        end
    end
end
close(fig)