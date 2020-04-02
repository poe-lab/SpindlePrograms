function SpindleStats_SpikesVsNot_08212017(UnitsSpikes)

%% Load tetrode data:
working_dir=pwd;
current_dir='C:\SleepData';
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

%% Select cell to analyze:
targetIdx = ismember(spikeCellNumbers, UnitsSpikes);
clear spikeCellNumbers
spikeTimeStamps = spikeTimeStamps(targetIdx);
clear targetIdx

%% Load spindle data from .MAT file:
working_dir=pwd;
current_dir='C:\SleepData';
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
clear spindlePath spindleFileName


%% Find all spindles that have spikes:
numSpikes = length(spikeTimeStamps);
%spikesInSpindles = zeros(numSpikes,1);
numSpindles = size(spindleStartTS, 1);
spindlesWithSpikes = zeros(numSpindles, 1);
spikesPerSpindle = spindlesWithSpikes;
for i = 1:numSpindles
    subIntervalIndx = find(spikeTimeStamps >= spindleStartTS(i) & spikeTimeStamps <= spindleStopTS(i));
    if ~isempty(subIntervalIndx)
        spindlesWithSpikes(i) = 1;
        spikesPerSpindle(i) = length(subIntervalIndx);
%         spikesInSpindles(subIntervalIndx) = 1;
        clear subIntervalIndx
    end
end

spindlesWoutSpikes = spindlesWithSpikes == 0;
spindlesWithSpikes = spindlesWithSpikes == 1;

%% Isolate spindle activity for spindles without spikes:
spindlesNoSpikes.duration = EschenkoSpindle.duration(spindlesWoutSpikes);
spindlesNoSpikes.maxP2pAmp = EschenkoSpindle.maxP2pAmp(spindlesWoutSpikes);
spindlesNoSpikes.symmetry = EschenkoSpindle.symmetry(spindlesWoutSpikes);
spindlesNoSpikes.power = EschenkoSpindle.power(spindlesWoutSpikes);
spindlesNoSpikes.peakFreq = EschenkoSpindle.peakFreq(spindlesWoutSpikes);

%% Isolate spindle activity for spindles with spikes:
spindlesYesSpikes.duration = EschenkoSpindle.duration(spindlesWithSpikes);
spindlesYesSpikes.maxP2pAmp = EschenkoSpindle.maxP2pAmp(spindlesWithSpikes);
spindlesYesSpikes.symmetry = EschenkoSpindle.symmetry(spindlesWithSpikes);
spindlesYesSpikes.power = EschenkoSpindle.power(spindlesWithSpikes);
spindlesYesSpikes.peakFreq = EschenkoSpindle.peakFreq(spindlesWithSpikes);

%% Calculate the statistics:
% [yesMean yesStdDev yesNum noMean noStdDev noNum tTestLogic tTestP]
statsResults = zeros(5, 8);
statsResults(1, :) = [mean(spindlesYesSpikes.duration) std(spindlesYesSpikes.duration)
[h, p] = ttest2(spindlesNoSpikes.duration,spindlesYesSpikes.duration);
% CONTINUE HERE!!!!



% numSpikesInSpindles = sum(spikesInSpindles);
% spikesInSpindles = spikesInSpindles==1;
% spikeTimeStamps = spikeTimeStamps(spikesInSpindles);

% %% Find closest index of spike time within continuous signal time stamps:
% spikeIdx = zeros(numSpikesInSpindles,1);
% for i = 1:numSpikesInSpindles
%     spikeIdx(i) = find(TimeStamps <= spikeTimeStamps(i) , 1, 'last');
% end

% %% Find RMS values 3 seconds before and after each spike:
% numDataPts = floor(3 * (1/(TimeStamps(2)-TimeStamps(1))));
% spikeRMS = zeros(numSpikesInSpindles, numDataPts * 2);
% for i = 1:numSpikesInSpindles
%     spikeRMS(i,:) = rmsSignal((spikeIdx(i) - numDataPts + 1):(spikeIdx(i) + numDataPts));
% end
% 
% avgRms = mean(spikeRMS, 1);
% stdRms = std(spikeRMS, 1);

% newSpikeRms = zeros(numSpikesInSpindles,61);
% for i=1:61
%     newSpikeRms(:,i) = sum(spikeRMS(:,(10*(i-1)+1):(10*i)),2);
% end

