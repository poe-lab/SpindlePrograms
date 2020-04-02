function SpindleViewer_AutoDetectedData_10152018
%% Select .mat file:
[spindleFilename, spindleFilePath] = uigetfile({'*.mat',...
        'Pick spindle data files.'},'Select .MAT File');
spindleFile = fullfile(spindleFilePath, spindleFilename);
load(spindleFile, '-mat')

%% Determine whether or not data includes startTime and endTime variables from spike and spindle analyses:
% If the variables do not exist, create from the time stamp vector:
if exist('startTime', 'var') == 0 
    startTime = TimeStamps(1) * 10^6;
    endTime = TimeStamps(end)  * 10^6;
end
clear TimeStamps

%% Select CSC file:
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[TimeStamps, SampleFrequencies, Samples, Header] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 1, 4, [startTime endTime] );

% Set the sampling rate:
orginalSampFreq = SampleFrequencies(1);
clear SampleFrequencies

%% Format the LFP/EEG data and interpolate time stamps:
[nsamp,eelen]=size(Samples);

%Reshape the LFP into a single column vector:
newM = nsamp*eelen;
Samples = reshape(Samples, newM, 1);

% Find AD bit volts in Header:
targ= strfind(Header,'-ADBitVolts');
for i=1:length(targ)
    targIdx(i)= isempty(targ{i}); %#ok<AGROW>
end
ADBitVoltsIdx = find(targIdx==0);   
clear targ targIdx
ADBitVolts = str2double(strrep(Header{ADBitVoltsIdx,1}, '-ADBitVolts ', '')); %#ok<FNDSB>

%Convert to microvolts:
signal = Samples * ADBitVolts *10^6;
clear ADBit2uV Samples

%Interpolate time stamps:
interpTimestamps = zeros(eelen*nsamp, 1);
idx = 1;
for i = 1:eelen
  if i < eelen
    t1 = TimeStamps(i);
    t2 = TimeStamps(i+1);
    interval = (t2-t1)/nsamp;
    trange =t1 : interval : t2-interval;
    interpTimestamps(idx:idx+nsamp-1,1) = trange;
  else
    t1 = TimeStamps(i);
    t2 = t1+interval*nsamp;
    trange =(t1 :interval : t2-interval);
    interpTimestamps(idx:idx+nsamp-1,1) = trange;
  end
  idx = idx + nsamp;
end
clear TimeStamps

% Convert from usec to seconds:
TimeStamps = interpTimestamps/1000000;
clear interpTimestamps

%% Calculate downsample target sampling rate:
DS = (1:1:10);
DSampSF = orginalSampFreq./DS;
indSampfactor = find(DSampSF >= 1000);
downsampFactor = DS(indSampfactor(end)); % Determine downsampling factor to get to target sampling frequency
newSampFreq = orginalSampFreq/downsampFactor;  % New downsampled sampling frequency     

%% Design filter for LFP/EEG signal: 
highPassFreq = [];
lowPassFreq = 30;
notchSetting = 1;
[recordHighPass, recordLowPass] = nlxFilterSettings(Header);
clear Header

[highPassFreq, lowPassFreq, sos, g, notchSetting] = filterSettingsCheck(highPassFreq, lowPassFreq, recordHighPass, recordLowPass, orginalSampFreq, newSampFreq, notchSetting);

%% Filter LFP/EEG data:
if ~isempty(g)
    signal = filtfilt(sos, g, signal);
end

%  OPTIONAL 60Hz Notch filter for LFP/EEG signal
if notchSetting ~= 0
    wo = 60/(orginalSampFreq/2);
    [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);
    signal = filtfilt(B_EEG_Notch,A_EEG_Notch, signal);
end

%% Downsample data if needed:
if downsampFactor > 1
    TimeStamps = TimeStamps(1:downsampFactor:end);
    signal = signal(1:downsampFactor:end);
end
clear physInput ADBit2uV 

%% Design filter for sigma band: 
sigmaHighpass = 10;
signmaLowpass = 15;
notchSetting = 0;
[~, ~, sos, g, ~] = filterSettingsCheck(sigmaHighpass, signmaLowpass, highPassFreq, lowPassFreq, newSampFreq, newSampFreq, notchSetting);

%% Filter Signal for Spindle Frequency:
if ~isempty(g)
    spindleSignal = filtfilt(sos, g, signal);
end

%% Set length of viewing window in seconds:
viewTimeWindow = [];
while isempty(viewTimeWindow)
    prompt={'Viewing Window Setting:'};
    dlgTitle='Set length of window (sec)';
    lineNo=1;
    defaultans = {'3'};
    options.Resize = 'on';
    answer = inputdlg(prompt,dlgTitle,lineNo, defaultans, options);
    viewTimeWindow = str2double(answer{1,1});
    clear answer prompt dlgTitle lineNo
end

%% Calculate & plot results for each spindle centered in a defined time window:
numSpindles = size(EschenkoSpindle.startIdx,1); %Determine number of detected spindles on the selected LFP/EEG channel
% Determine maximum y-axis values (mV) for the spindle filtered signals plus
% added 10%
maxY_Spindle = max(abs(spindleSignal)) + 0.1*max(abs(spindleSignal));

for i = 1:numSpindles
    eventStartTS = EschenkoSpindle.timestamp(i);
    eventStopTS = EschenkoSpindle.timestamp(i) + EschenkoSpindle.duration(i);
    addTime = viewTimeWindow - (eventStopTS - eventStartTS);
    if addTime < 0   % If the spindle duration >= viewTimeWindow,
        addTime = 0; %do not add any data points to the beginning and end.
    end
    % Define the time range of the viewing window for the spindle:
    windowStartTime = eventStartTS - addTime/2;
    windowStopTime = eventStopTS + addTime/2;
    
    ax1 = subplot(1,1,1);  % creates a figure with axes
    
    % Find indices of time points within the viewing window:
    windowTargetIdx = TimeStamps >= windowStartTime & TimeStamps <= windowStopTime;
    
    % Isolate windowed data:
    windowTS = TimeStamps(windowTargetIdx);
    windowSpindleSignal = spindleSignal(windowTargetIdx);
    windowLFP = signal(windowTargetIdx);
    clear windowTargetIdx
    
    % Plot the sigma filtered signal for the viewing widow:
    plot(ax1, windowTS, windowSpindleSignal, 'b')
    hold on % keeps plot open for additional signals
    
    % Find indices of time points for the spindle within the viewing window:
    eventTargetIdx = windowTS >= eventStartTS & windowTS <= eventStopTS;
    
    % Plot the spindle in green:
    plot(ax1, windowTS(eventTargetIdx), windowSpindleSignal(eventTargetIdx), 'g')
    
    % Plot the LFP/EEG signal for the viewing widow in gray:
    plot(ax1, windowTS, windowLFP, 'Color', [0.5 0.5 0.5])

    hold off
    switch EschenkoSpindle.scoring(i)
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
    title(['Spindle ' num2str(i) ' of ' num2str(numSpindles) '          Stage: ' stage])
    xlim([windowTS(1) windowTS(end)])
    ylim([-maxY_Spindle maxY_Spindle])
    ylabel(ax1, 'LFP')
    pause
end