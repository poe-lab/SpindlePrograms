function spikeSpindlePowerAnalysis_12202017

%% Select Stage Scored File:
working_dir=pwd;
current_dir='C:\';
cd(current_dir);
scoredCheck = 0;
while isequal(scoredCheck, 0)
    [scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},...
        'Select the Sleep Scored File');
    if isequal(scoredFile,0) || isequal(scoredPath,0)
        uiwait(errordlg('You need to select a file. Please try again',...
            'ERROR','modal'));
    else
        cd(working_dir);
        stageScoredFile= fullfile(scoredPath, scoredFile);
        %Load sleep scored file:
        try
            [numData, stringData] = xlsread(stageScoredFile);
            scoredCheck = 1;
        catch %#ok<*CTCH>
            % If file fails to load, it will notify user and prompt to
            % choose another file.
            uiwait(errordlg('Check if the scored file is saved in Microsoft Excel format.',...
             'ERROR','modal'));
         scoredCheck = 0;
        end

    end
end

%% Detect if states are in number or 2-letter format:
if isequal(size(numData,2),3)
    scoredStates = numData(:,2:3);
    clear numData stringData
else
    scoredStates = numData(:,2);
    clear numData
    stringData = stringData(3:end,3);
    [stateNumber] = stateLetter2NumberConverter(stringData);
    scoredStates = [scoredStates stateNumber];
    clear stateNumber stringData
end
epochInSeconds = scoredStates(2,1) - scoredStates(1,1);
startTime = scoredStates(1,1) * 10^6;
endTime = (scoredStates(end,1) + epochInSeconds) * 10^6;

%% Load spindle data from .MAT file:
working_dir=pwd;
current_dir='C:\';
cd(current_dir);

[spindleFile, spindlePath] = uigetfile({'*.mat',...
        'Detected spindles file (*.MAT)'},'Select the spindle data file:');
if isequal(spindleFile,0) || isequal(spindlePath,0)
    uiwait(errordlg('You need to select a file. Please try again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    spindleFileName = fullfile(spindlePath, spindleFile);
end

load(spindleFileName, '-mat')
targetIdx = EschenkoSpindle.scoring==2 | EschenkoSpindle.scoring==6;
spindleStartTS = TimeStamps(EschenkoSpindle.startIdx(targetIdx));
spindleStopTS = TimeStamps(EschenkoSpindle.stopIdx(targetIdx));
% clear EschenkoSpindle spindlePath spindleFileName CSCFilename scoredFile

%% Select CSC file:
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[TimeStamps, SampleFrequencies, Samples, Header] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 1, 4, [startTime endTime] );

%% Set constants:
Fs = SampleFrequencies(1);
clear SampleFrequencies

% Find AD bit volts in Header:
targ= strfind(Header,'-ADBitVolts');
for i=1:length(targ)
    targIdx(i)= isempty(targ{i}); %#ok<AGROW>
end
ADBitVoltsIdx = find(targIdx==0);   
clear targ targIdx
ADBitVolts = str2double(strrep(Header{ADBitVoltsIdx,1}, '-ADBitVolts ', '')); %#ok<FNDSB>
clear Header

[nsamp,eelen]=size(Samples);

%Reshape the LFP into a single column vector:
newM = nsamp*eelen;
Samples = reshape(Samples, newM, 1);
signal = Samples * ADBitVolts *10^6; %Convert to microvolts
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

%% Filter Signal for EEG/LFP:
if Fs > 1050
    [z, p, k] = ellip(7,1,60, 300/(Fs/2),'low');
    [sos, g] = zp2sos(z,p,k);
    signal = filtfilt(sos, g, signal);
end

%% Downsample Signal:
targetFs = 1000;
DsFactor = floor(Fs/targetFs); % Determine downsampling factor to get to target sampling frequency
newFs = Fs/DsFactor;  % New downsampled sampling frequency     
TimeStamps = TimeStamps(1:DsFactor:end);
signal = signal(1:DsFactor:end);

%% Filter Signal for Spindle Frequency:
sigmaHighpass = 10;
signmaLowpass = 15;
[z, p, k] = ellip(7,1,60, [sigmaHighpass signmaLowpass]/(newFs/2),'bandpass');
[sos, g] = zp2sos(z,p,k);
filtSigma = filtfilt(sos, g, signal);
clear signal

%% Downsample Signal:
targetFs = 100;
DsFactor = floor(newFs/targetFs); % Determine downsampling factor to get to target sampling frequency
newFs = newFs/DsFactor;  % New downsampled sampling frequency     
TimeStamps = TimeStamps(1:DsFactor:end);
filtSigma = filtSigma(1:DsFactor:end);

% %% Assign states to each data point as a new vector:
% lengthSignal = length(filtSigma);
% sleepsamp = zeros(lengthSignal,1);
% for i = 1:size(scoredStates, 1)
%     if isequal(i, size(scoredStates, 1))
%         subIntervalIndx = find(TimeStamps >= scoredStates(i,1) & TimeStamps < (scoredStates(i,1) + 10));
%     else
%         subIntervalIndx = find(TimeStamps >= scoredStates(i,1) & TimeStamps < scoredStates(i+1,1));
%     end
%     if ~isempty(subIntervalIndx)
%         sleepsamp(subIntervalIndx) = scoredStates(i,2);
%         clear subIntervalIndx
%     end
% end

%% Calculate the root mean square at each point with a 0.1 second window size
lengthSignal = length(filtSigma);
squaredSignal = filtSigma.^2;
clear filtSigma
windowSize = 0.1 * newFs;
halfWindow = floor(windowSize/2);
rmsStart = halfWindow+1;
rmsEnd = lengthSignal - halfWindow;
rmsSignal = zeros(lengthSignal,1);
for i = rmsStart:rmsEnd
    rmsSignal(i) = sqrt(sum(squaredSignal(i-halfWindow:i+halfWindow))/(2*halfWindow+1));
end
clear squaredSignal

%% Load tetrode data:
working_dir=pwd;
current_dir='C:\';
cd(current_dir);

[spikeFile, spikePath] = uigetfile({'*.ntt',...
        'Sorted Neuralynx Tetrode File (*.NTT)'},'Select a Spike Sorted Data File');
if isequal(spikeFile,0) || isequal(spikePath,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    spikeFileName = fullfile(spikePath, spikeFile);
end

[spikeTimeStamps, spikeCellNumbers] = Nlx2MatSpike(spikeFileName, [1 0 1 0 0], 0, 1, [] );
clear spikeFileName spikePath
nonZerosIndex = find(spikeCellNumbers);
spikeCellNumbers = spikeCellNumbers(nonZerosIndex)';
spikeTimeStamps = spikeTimeStamps(nonZerosIndex)';
clear nonZerosIndex
spikeTimeStamps = spikeTimeStamps/1000000;
numOfCells = max(spikeCellNumbers);

%% Select cell to analyze:
cellCheck = 0;
while isequal(cellCheck, 0)
    prompt = {['Choose a cell to be analyzed (max=' num2str(numOfCells) ')']};
    def = {'0'};
    dlgTitle = 'Select Cell Number';
    lineNo = 1;
    answer = inputdlg(prompt,dlgTitle,lineNo,def);
    cellNum = str2double(char(answer(1,:)));
    if isempty(cellNum) || cellNum > numOfCells
        
    else
        cellCheck = 1;
    end
end
targetIdx = spikeCellNumbers == cellNum;
clear spikeCellNumbers cellCheck prompt def dlgTitle lineNo answer numOfCells
spikeTimeStamps = spikeTimeStamps(targetIdx);
clear targetIdx

% %% Load spindle data from .MAT file:
% working_dir=pwd;
% current_dir='C:\SleepData';
% cd(current_dir);
% 
% [spindleFile, spindlePath] = uigetfile({'*.mat',...
%         'Detected spindles file (*.MAT)'},'Select the spindle data file:');
% if isequal(spindleFile,0) || isequal(spindlePath,0)
%     uiwait(errordlg('You need to select a file. Please try again',...
%         'ERROR','modal'));
%     cd(working_dir);
% else
%     cd(working_dir);
%     spindleFileName = fullfile(spindlePath, spindleFile);
% end
% 
% load(spindleFileName, '-mat')
% targetIdx = EschenkoSpindle.scoring==2 | EschenkoSpindle.scoring==6;
% spindleStartTS = TimeStamps(EschenkoSpindle.startIdx(targetIdx));
% spindleStopTS = TimeStamps(EschenkoSpindle.stopIdx(targetIdx));
% clear EschenkoSpindle spindlePath spindleFileName CSCFilename scoredFile

% %% Call the sleep scored file name
% working_dir=pwd;
% current_dir='C:\SleepData';
% cd(current_dir);
% scoredCheck = 0;
% while isequal(scoredCheck, 0)
%     [scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},...
%         'Select the Sleep Scored File');
%     if isequal(scoredFile,0) || isequal(scoredPath,0)
%         uiwait(errordlg('You need to select a file. Please try again',...
%             'ERROR','modal'));
%     else
%         cd(working_dir);
%         stageScoredFile= fullfile(scoredPath, scoredFile);
%         %Load sleep scored file:
%         try
%             [numData, stringData] = xlsread(stageScoredFile);
%             scoredCheck = 1;
%         catch %#ok<*CTCH>
%             % If file fails to load, it will notify user and prompt to
%             % choose another file.
%             uiwait(errordlg('Check if the scored file is saved in Microsoft Excel format.',...
%              'ERROR','modal'));
%          scoredCheck = 0;
%         end
% 
%     end
% end
% 
% %% Detect if states are in number or 2-letter format:
% if isequal(size(numData,2),3)
%     scoredStates = numData(:,2:3);
%     clear numData stringData
% else
%     scoredStates = numData(:,2);
%     clear numData
%     stringData = stringData(3:end,3);
%     [stateNumber] = stateLetter2NumberConverter(stringData);
%     scoredStates = [scoredStates stateNumber];
%     clear stateNumber stringData
% end

%% Isolate NREM epochs:
scoredStates = scoredStates(scoredStates(:,2)==2 | scoredStates(:,2)==6, :);

%% Find all spikes within NREM sleep:
numSpikes = length(spikeTimeStamps);
spikesInNREM = zeros(numSpikes,1);
for i = 1:size(scoredStates, 1)
    if isequal(i, size(scoredStates, 1))
        subIntervalIndx = find(spikeTimeStamps >= scoredStates(i,1) & spikeTimeStamps < (scoredStates(i,1) + 10));
    else
        subIntervalIndx = find(spikeTimeStamps >= scoredStates(i,1) & spikeTimeStamps < scoredStates(i+1,1));
    end
    if ~isempty(subIntervalIndx)
        spikesInNREM(subIntervalIndx) = 1;
        clear subIntervalIndx
    end
end

numSpikesInNREM = sum(spikesInNREM);
spikesInNREM = spikesInNREM==1;
spikeTimeStamps = spikeTimeStamps(spikesInNREM);

%% Find all spikes within spindles:
numSpikes = length(spikeTimeStamps);
spikesInSpindles = zeros(numSpikes,1);
for i = 1:size(spindleStartTS, 1)
    subIntervalIndx = find(spikeTimeStamps >= spindleStartTS(i) & spikeTimeStamps <= spindleStopTS(i));
    if ~isempty(subIntervalIndx)
        spikesInSpindles(subIntervalIndx) = 1;
        clear subIntervalIndx
    end
end

numSpikesInSpindles = sum(spikesInSpindles);
numSpikesNotInSpindles = numSpikes - numSpikesInSpindles;
spikesNotInSpindles = spikesInSpindles==0;
spikesNotInSpindlesTS = spikeTimeStamps(spikesNotInSpindles);
spikesInSpindles = spikesInSpindles==1;
spikeTimeStamps = spikeTimeStamps(spikesInSpindles);

%% Find closest index of spike time within continuous signal time stamps:
spikeIdx = zeros(numSpikesInSpindles,1);
spikeIdx = nearestpoint(spikeTimeStamps,TimeStamps);
% for i = 1:numSpikesInSpindles
%     spikeIdx(i) = find(TimeStamps <= spikeTimeStamps(i) , 1, 'last');
% end

spikeNoSpindleIdx = zeros(numSpikesNotInSpindles,1);
spikeNoSpindleIdx = nearestpoint(spikesNotInSpindlesTS,TimeStamps);
% for i = 1:numSpikesNotInSpindles
%     spikeNoSpindleIdx(i) = find(TimeStamps <= spikesNotInSpindlesTS(i) , 1, 'last');
% end

%% Find RMS values 3 seconds before and after each spike:
numDataPts = floor(3 * (1/(TimeStamps(2)-TimeStamps(1))));
spikeRMS = zeros(numSpikesInSpindles, numDataPts * 2);
for i = 1:numSpikesInSpindles
    spikeRMS(i,:) = rmsSignal((spikeIdx(i) - numDataPts + 1):(spikeIdx(i) + numDataPts));
end

avgRms = mean(spikeRMS, 1);
stdRms = std(spikeRMS, 1);

spikeNoSpindRMS = zeros(numSpikesNotInSpindles, numDataPts * 2);
for i = 1:numSpikesNotInSpindles
    spikeNoSpindRMS(i,:) = rmsSignal((spikeNoSpindleIdx(i) - numDataPts + 1):(spikeNoSpindleIdx(i) + numDataPts));
end

avgRmsSpikeNoSpind = mean(spikeNoSpindRMS, 1);
stdRmsSpikeNoSpind = std(spikeNoSpindRMS, 1);

% newSpikeRms = zeros(numSpikesInSpindles,61);
% for i=1:61
%     newSpikeRms(:,i) = sum(spikeRMS(:,(10*(i-1)+1):(10*i)),2);
% end
% avgRms, avgRmsSpikeNoSpind, cellNum, CSCFilename, numSpikesInSpindles, numSpikesNotInSpindles, scoredFile, spikeFile, spikeNoSpindRMS, spikeRMS, spindleFile, stdRms, stdRmsSpikeNoSpind
