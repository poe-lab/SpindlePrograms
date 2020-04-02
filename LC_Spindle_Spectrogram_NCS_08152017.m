function LC_Spindle_Spectrogram_NCS_08152017
% Find spikes within spindles and calculate spectrogram for each spike.

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
clear EschenkoSpindle spindlePath spindleFileName CSCFilename scoredFile

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
spikesInSpindles = spikesInSpindles==1;
spikeTimeStamps = spikeTimeStamps(spikesInSpindles);

%% Select CSC file:
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[~, SampleFrequencies, Samples, Header] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 1, 1, [] );

%% Set constants:
Fs = SampleFrequencies(1);
clear SampleFrequencies

%% Downsample Signal:
% Automated down-sampling to 1kHz if greater or output an error if below 100 Hz:
if Fs >= 2000
    DS = (1:1:10);
    DSampSF = cscSamplingRate./DS;
    indSampfactor = find(DSampSF >= 1000);
    newFs = DSampSF(indSampfactor(end)); % New downsampled sampling frequency
    DsFactor = DS(indSampfactor(end));  % Determine downsampling factor to get to target sampling frequency
    msgbox({['Recording Sampling Rate:  ' num2str(Fs) 'Hz'];...
        ['Down-Sampled Sampling Rate:  ' num2str(newFs) 'Hz'];...
        ['Sampling Factor:  ' num2str(DsFactor) '']});
    clear DS DSampSF indSampfactor
elseif Fs < 100
    msgbox({['Recording Sampling Rate:  ' num2str(Fs) 'Hz'];...
        'This sampling rate is too low to detect spindles.';...
        'The program will now close.'});
    exit
else  % In practice, we should be sampling at at least 1kHz, but older files may have been recorded at slower frequencies.
    newFs= Fs;
    DsFactor = 1;
    msgbox({['Recording Sampling Rate:  ' num2str(Fs) 'Hz'];...
        'This sampling rate is acceptable to detect spindles.'});
end

%% Find AD bit volts in Header:
targ= strfind(Header,'-ADBitVolts');
for i=1:length(targ)
    targIdx(i)= isempty(targ{i}); %#ok<AGROW>
end
ADBitVoltsIdx = find(targIdx==0);   
clear targ targIdx
ADBitVolts = str2double(strrep(Header{ADBitVoltsIdx,1}, '-ADBitVolts ', '')); %#ok<FNDSB>
clear Header ADBitVoltsIdx

%% Reshape the LFP into a single column vector:
[nsamp,eelen]=size(Samples);
newM = nsamp*eelen;
Samples = reshape(Samples, newM, 1);
signal = Samples * ADBitVolts *10^6; %Convert to microvolts
clear ADBit2uV Samples

%% Interpolate time stamps:
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
  clear t1 t2 trange
  idx = idx + nsamp;
end
clear TimeStamps nsamp eelen idx

% Convert from usec to seconds:
TimeStamps = interpTimestamps/1000000;
clear interpTimestamps

%% Filter Signal for EEG/LFP:
[z, p, k] = ellip(7,1,60, 30/(Fs/2),'low');
[sos, g] = zp2sos(z,p,k);
eegSignal = filtfilt(sos, g, signal);
clear signal z p k sos g

%% Downsample signal data  to target sampling rate if necessary:
TimeStamps = TimeStamps(1:DsFactor:end);
eegSignal = eegSignal(1:DsFactor:end);
clear DsFactor

%% Find closest index of spike time within continuous signal time stamps:
spikeIdx = zeros(numSpikesInSpindles,1);
for i = 1:numSpikesInSpindles
    spikeIdx(i) = find(TimeStamps <= spikeTimeStamps(i) , 1, 'last');
end

%% Find data points 3 seconds before and after each spike:
numDataPts = floor(3 * (1/(TimeStamps(2)-TimeStamps(1))));
spikeSignal = zeros(numSpikesInSpindles, numDataPts * 2);
for i = 1:numSpikesInSpindles
    spikeSignal(i,:) = eegSignal((spikeIdx(i) - numDataPts + 1):(spikeIdx(i) + numDataPts));
end


% winSize = 0.5;
% winOverlap = 0.25;
fStart = ;

fStop = 17;
fRes = 0.1;

% window = newFs * winSize;
% noverlap = floor(winOverlap * newFs);
window = 500;
noverlap = 490;

frequencyRange= fStart:fRes:fStop;
eegData = spikeSignal(1,:);
[~, F, T, sumP] = spectrogram(eegData,window,noverlap,frequencyRange,newFs);

for i = 2:numSpikesInSpindles
    eegData = spikeSignal(i,:);
    [~, F, T, P] = spectrogram(eegData,window,noverlap,frequencyRange,newFs);
    sumP = sumP + P;
end
avgP = sumP./numSpikesInSpindles;
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