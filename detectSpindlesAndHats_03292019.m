function detectSpindlesAndHats_03292019

% Detect sleep spindles based on algorithm in Eschenko et al., "Elevated 
% Sleep Spindle Density after Learning or after Retrieval in Rats ", 
% J. Neurosci., December 13, 2006. 26(50):12914–12920.

% 11.16.2017: Changed parameters for spindle detection consistent with
% Kevin's manual detection parameters. -BAG

% 03.29.2019: Created new program to detect spindles or HATS with two-threshold method.
% -BAG
%% Select Stage Scored File:
working_dir=pwd;
current_dir='C:\';
cd(current_dir);
scoredCheck = 0;
while isequal(scoredCheck, 0)
    [info.scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},...
        'Select the Sleep Scored File');
    if isequal(info.scoredFile,0) || isequal(scoredPath,0)
        uiwait(errordlg('You need to select a file. Please try again',...
            'ERROR','modal'));
    else
        cd(working_dir);
        stageScoredFile= fullfile(scoredPath, info.scoredFile);
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

%% Select CSC file:
[info.CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, info.CSCFilename);

%% Load signal data:

[TimeStamps, SampleFrequencies, Samples, Header] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 1, 4, [startTime endTime] );
Fs = SampleFrequencies(1);
clear SampleFrequencies

%% Find AD bit volts in Header:
targ= strfind(Header,'-ADBitVolts');
for i=1:length(targ)
    targIdx(i)= isempty(targ{i}); %#ok<AGROW>
end
ADBitVoltsIdx = find(targIdx==0);   
clear targ targIdx
ADBitVolts = str2double(strrep(Header{ADBitVoltsIdx,1}, '-ADBitVolts ', '')); %#ok<FNDSB>

%% Find recording filter settings in header:
[hiPass, loPass] = nlxFilterSettings(Header);
clear Header

%% Set frequencies to reduce bandwidth for down-sampling:
info.HighPassFreq = 1.0; % in Hz
info.LowPassFreq = 30; % in Hz
signal_Notch_enable = 0;

%% Determine if signal should be down-sampling
targetFs = 10 * info.LowPassFreq;
DsFactor = floor(Fs/targetFs); % Determine downsampling factor to get to target sampling frequency
newFs = Fs/DsFactor;  % New downsampled sampling frequency



%% Determine correct filter settings:
[info.HighPassFreq, info.LowPassFreq, signal_sos, signal_g, signal_Notch_enable] =...
    filterSettingsCheck(info.HighPassFreq, info.LowPassFreq, hiPass, loPass, Fs, newFs, signal_Notch_enable);
clear hiPass loPass
    

% waitbar(0.6,waithandle,' Converting EEG from Neuralynx CSC to Matlab data ...');
% figure(waithandle),pause(0.2)


%% Reshape the signal into a single column vector & convert to microvolts:
[nsamp,eelen]=size(Samples);
newM = nsamp*eelen;
Samples = reshape(Samples, newM, 1);
signal = Samples * ADBitVolts *10^6; %Convert to microvolts
clear ADBit2uV Samples

%% Interpolate time stamps:
[preciseTS, preciseSamples, medianSampRate] = interpolateTS_Ncsfiles_01242019(TimeStamps,signal, [], [], []);
clear TimeStamps signal

%% Filter signal:
% waitbar(0.8,waithandle,'Filtering the signal ...'); 
% figure(waithandle),pause(0.2),
if ~isempty(signal_g)
    filtered_samples = filtfilt(signal_sos, signal_g, preciseSamples);
else
    filtered_samples = preciseSamples;
end
clear preciseSamples

%  OPTIONAL 60Hz Notch filter for EEG signals
if signal_Notch_enable ~= 0
    wo = 60/(Fs/2);
    [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
    filtered_samples = filtfilt(B_EEG_Notch,A_EEG_Notch, filtered_samples);
end


%% Downsample Signal:
info.TimeStamps = preciseTS(1:DsFactor:end);
clear preciseTS
info.Signal = filtered_samples(1:DsFactor:end);
clear filtered_samples

%% Assign states to each data point as a new vector:
lengthSignal = length(info.Signal);
info.StateArray = zeros(lengthSignal,1);
for i = 1:size(scoredStates, 1)
    if isequal(i, size(scoredStates, 1))
        subIntervalIndx = find(info.TimeStamps >= scoredStates(i,1) & info.TimeStamps < (scoredStates(i,1) + 10));
    else
        subIntervalIndx = find(info.TimeStamps >= scoredStates(i,1) & info.TimeStamps < scoredStates(i+1,1));
    end
    if ~isempty(subIntervalIndx)
        info.StateArray(subIntervalIndx) = scoredStates(i,2);
        clear subIntervalIndx
    end
end

%% Save data to MATLAB file:
%Request user to name output file:
prompt = {'Enter the filename you want to save it as: (just the name)'};
def = {'Rat#_Day'};
dlgTitle = 'Save .MAT file';
lineNo = 1;
answer = inputdlg(prompt,dlgTitle,lineNo,def);
filename = char(answer(1,:));
matFileName = strcat(filename,'_SpindlesAndHats.mat');

% Determine if .MAT file already exists:
matFile = fullfile(CSCFilePath, matFileName);
if exist(matFile, 'file') == 2
    save(matFile, 'info', '-append');
else
    save(matFile, 'info');
end

%% SPINDLE DETECTION
%% Filter Signal for Spindle Frequency:
% Define frequency parameters:
spindle.Highpass = 10;
spindle.Lowpass = 15;
spindle.NotchEnable = 0;

[spindle.Highpass, spindle.Lowpass, spindle.sos, spindle.g, spindle.NotchEnable] =...
    filterSettingsCheck(spindle.Highpass, spindle.Lowpass, info.HighPassFreq, info.LowPassFreq, newFs, newFs, spindle.NotchEnable);

% Filter signal:
% waitbar(0.8,waithandle,'Filtering signal for spindle detection ...'); 
% figure(waithandle),pause(0.2),
if ~isempty(spindle.g)
    spindle.Signal = filtfilt(spindle.sos, spindle.g, info.Signal);
else
    spindle.Signal = info.Signal;
end


%% Calculate the RMS of the sigma-filtered signal:
% Window size in seconds used to calculate RMS:
windowSeconds = 0.1; 
spindle.RMS = calculateRMS(spindle.Signal, windowSeconds, newFs);

%% Detect spindles using 2-threshold of RMS signal method:
% Specify threshold multipliers:
spindle.upperCrit = 5;     
spindle.lowerCrit = 2;
% Enter states to be used to calculate thresholds:
spindle.targetStatesThreshold = [2 6]; 
% Enter states to look for spindles:
spindle.targetStatesDetect = [2 3 6];
% Define max time between spindles to combine adjacent spindles in seconds:
spindle.maxCombineTime = 0.1;
% Define range for spindle length in seconds:
spindle.durationRange = [0.3 3.0];

spindle = twoThresholdRmsEventDetect(spindle, newFs, info.StateArray, info.TimeStamps);

%% Save data to MATLAB file:
if exist(matFile, 'file') == 2
    save(matFile, 'spindle', '-append');
else
    save(matFile, 'spindle');
end
clear spindle

%% HAT DETECTION
%% Filter Signal for HAT Frequency:
% Define frequency parameters:
HAT.Highpass = 5;
HAT.Lowpass = 9;
HAT.NotchEnable = 0;

[HAT.Highpass, HAT.Lowpass, HAT.sos, HAT.g, HAT.NotchEnable] =...
    filterSettingsCheck(HAT.Highpass, HAT.Lowpass, info.HighPassFreq, info.LowPassFreq, newFs, newFs, HAT.NotchEnable);

% Filter signal:
% waitbar(0.8,waithandle,'Filtering signal for HAT detection ...'); 
% figure(waithandle),pause(0.2),
if ~isempty(HAT.g)
    HAT.Signal = filtfilt(HAT.sos, HAT.g, info.Signal);
else
    HAT.Signal = info.Signal;
end


%% Calculate the RMS of the sigma-filtered signal:
% Window size in seconds used to calculate RMS:
windowSeconds = 0.1; 
HAT.RMS = calculateRMS(HAT.Signal, windowSeconds, newFs);

%% Detect HATs using 2-threshold of RMS signal method:
% Specify threshold multipliers:
HAT.upperCrit = 5;     
HAT.lowerCrit = 2;
% Enter states to be used to calculate thresholds:
HAT.targetStatesThreshold = [3]; 
% Enter states to look for HATs:
HAT.targetStatesDetect = [2 3 6];
% Define max time between HATs to combine adjacent spindles in seconds:
HAT.maxCombineTime = 0.1;
% Define range for HAT length in seconds:
HAT.durationRange = [0.5 4.0];

HAT = twoThresholdRmsEventDetect(HAT, newFs, info.StateArray, info.TimeStamps);

%% Save data to MATLAB file:
if exist(matFile, 'file') == 2
    save(matFile, 'HAT', '-append');
else
    save(matFile, 'HAT');
end
clear HAT
end