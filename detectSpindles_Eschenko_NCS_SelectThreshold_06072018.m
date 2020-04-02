function detectSpindles_Eschenko_NCS_SelectThreshold_06072018

% Detect sleep spindles based on algorithm in Eschenko et al., "Elevated 
% Sleep Spindle Density after Learning or after Retrieval in Rats ", 
% J. Neurosci., December 13, 2006. 26(50):12914–12920.

% 11.16.2017: Changed parameters for spindle detection consistent with
% Kevin's manual detection parameters. -BAG

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

%% Assign states to each data point as a new vector:
lengthSignal = length(filtSigma);
sleepsamp = zeros(lengthSignal,1);
for i = 1:size(scoredStates, 1)
    if isequal(i, size(scoredStates, 1))
        subIntervalIndx = find(TimeStamps >= scoredStates(i,1) & TimeStamps < (scoredStates(i,1) + 10));
    else
        subIntervalIndx = find(TimeStamps >= scoredStates(i,1) & TimeStamps < scoredStates(i+1,1));
    end
    if ~isempty(subIntervalIndx)
        sleepsamp(subIntervalIndx) = scoredStates(i,2);
        clear subIntervalIndx
    end
end

%% Calculate the root mean square at each point with a 0.1 second window size
squaredSignal = filtSigma.^2;
% clear filtSigma
windowSize = 0.1 * newFs;
halfWindow = floor(windowSize/2);
rmsStart = halfWindow+1;
rmsEnd = lengthSignal - halfWindow;
rmsSignal = zeros(lengthSignal,1);
for i = rmsStart:rmsEnd
    rmsSignal(i) = sqrt(sum(squaredSignal(i-halfWindow:i+halfWindow))/(2*halfWindow+1));
end
clear squaredSignal
% rmsSignal = rmsSignal(rmsStart:rmsEnd);
% TimeStamps = TimeStamps(rmsStart:rmsEnd);
% sleepsamp = sleepsamp(rmsStart:rmsEnd);
% filtSigma = filtSigma(rmsStart:rmsEnd);

%% Determine threshold
info.threshMultiRMS = [];
while isempty(info.threshMultiRMS)
    prompt={'Threshold = multiplier * RMS standard deviation:'};
    dlgTitle='Select multiplier';
    lineNo=1;
    defaultans = {'3'};
    options.Resize = 'on';
    answer = inputdlg(prompt,dlgTitle,lineNo, defaultans, options);
    info.threshMultiRMS = str2double(answer{1,1});
    clear answer prompt dlgTitle lineNo
end

% meanNremRMS = mean(rmsSignal(sleepsamp==2 | sleepsamp==6));
info.stdDevNremRMS = std(rmsSignal(sleepsamp==2 | sleepsamp==6));
threshold = info.threshMultiRMS * info.stdDevNremRMS;
info.threshold = [];
while isempty(info.threshold)
    prompt={'Threshold = multiplier * RMS standard deviation:'};
    dlgTitle='Select multiplier';
    lineNo=1;
    defaultans = {num2str(threshold)};
    options.Resize = 'on';
    answer = inputdlg(prompt,dlgTitle,lineNo, defaultans, options);
    info.threshold = str2double(answer{1,1});
    clear answer prompt dlgTitle lineNo
end

%% Find where RMS signal crosses threshold
aboveThreshold = rmsSignal > threshold & (sleepsamp==2 | sleepsamp==6);
clear rmsSignal

%% Find the beginning and end of each spindle
startSpindle = find(diff(aboveThreshold)==1)+1;
endSpindle = find(diff(aboveThreshold)==-1);
if endSpindle(1) < startSpindle(1)                                          % Corrects for spindles that start at the beginning of scored period
    startSpindle = [0; startSpindle];
end
if startSpindle(end) > endSpindle(end)                                      % Corrects for spindles that end at the end of scored period
    endSpindle = [endSpindle; lengthSignal];
end
% endSpindle = endSpindle(find(endSpindle > startSpindle(1), 1, 'first'):end);
% startSpindle = startSpindle(1:find(startSpindle < endSpindle(end), 1, 'last'));

%% Combine close spindles:
maxdist = 0.3;                                                              % Maximum distance between spindles to combine into same spindle
timeBetweenSpindles = startSpindle(2:end) - endSpindle(1:(end-1));          % Find the distances between all spindles
k = (timeBetweenSpindles <= maxdist);                                                         % Find when the time between spindles is less than the maximum
clear timeBetweenSpindles
k1 = [k; 0];
k2 = [0; k];                                                                % Add a zero for the first ripple that was not included
clear k
endSpindle(k1==1) = endSpindle(k2==1);                                      % Unite them (use the ending of the second ripple as end of the first one...
endSpindle(k2==1) = [];
startSpindle(k2==1) = [];                                                   % ...and delete the times of the second one)
clear k1 k2

%% Determine scored state for each spindle
scoring = zeros(length(endSpindle)); % sleep state where end lies
for i=1:length(endSpindle)
    scoring(i) = sleepsamp(endSpindle(i));
end
clear sleepsamp

% %% Filter for delta:
% deltaHighpass = 0.5;
% deltaLowpass = 4;
% [z, p, k] = ellip(7,1,60, [deltaHighpass deltaLowpass]/(newFs/2),'bandpass');
% [sos, g] = zp2sos(z,p,k);
% filtDelta = filtfilt(sos, g, signal);
% 
% %% Calculate delta RMS:
% squaredSignal = filtSigma.^2;
% % clear filtSigma
% windowSize = 0.1 * newFs; %0.1 second window size
% halfWindow = floor(windowSize/2);
% rmsStart = halfWindow+1;
% rmsEnd = lengthSignal - halfWindow;
% rmsSignal = zeros(lengthSignal,1);
% for i = rmsStart:rmsEnd
%     rmsSignal(i) = sqrt(sum(squaredSignal(i-halfWindow:i+halfWindow))/(2*halfWindow+1));
% end
% clear squaredSignal

%% Exclude spindles based on duration:
info.minDuration = 0.3;
info.maxDuration = 3;
duration = TimeStamps(endSpindle) - TimeStamps(startSpindle);
duration(duration< info.minDuration | info.maxDuration>5) = NaN;
spindle.startIdx = startSpindle(isnan(duration)~=1);
spindle.stopIdx = endSpindle(isnan(duration)~=1);
spindle.scoring = scoring(isnan(duration)~=1);
spindle.duration = duration(isnan(duration)~=1);
clear duration scoring startSpindle endSpindle

%% Calculate spindle activity parameters for detected spindles:
% sigmaHighpass = 11;
% signmaLowpass = 16;
% [z, p, k] = ellip(7,1,60, [sigmaHighpass signmaLowpass]/(newFs/2),'bandpass');
% [sos, g] = zp2sos(z,p,k);
% filtSigma = filtfilt(sos, g, eegSignal);

for i = 1:length(spindle.duration)
    x = filtSigma(spindle.startIdx(i):spindle.stopIdx(i)); % pull out spindle from signal

    % Find the max peak-to-peak amplitude
    [~, absIdx] = findpeaks(abs(x)); % finds all of the peaks and troughs
    xPeaks = x(absIdx); % Gets values of peaks and troughs
    p2pAmp = abs(diff(xPeaks)); % finds all peak-to-peak amplitudes
    [spindle.maxP2pAmp(i), p2pMaxIdx] = max(p2pAmp); 
    spindle.symmetry(i) = absIdx(p2pMaxIdx)/length(x);
    clear absIdx xPeaks p2pAmp p2pMaxIdx

    % Compute the power spectrum of the Hann tapered data:
    dt = 1/newFs; % Define the sampling interval.
    df = 1/spindle.duration(i); % Determine the frequency resolution.
    fNyquist = newFs/2; % Determine the Nyquist frequency.
    faxis = (0:df:fNyquist); % Construct the frequency axis.
    xh = hann(length(x)).*x;
    Sxx = 2*dt^2/spindle.duration(i) * fft(xh).*conj(fft(xh)); % Compute the power spectrum of Hann tapered data.
    Sxx = real(Sxx); %Ignores imaginary component.
    Sxx = Sxx(1:length(faxis));
    clear dt df fNyquist xh

    % Calculate bandpower
    spindle.power(i) = bandpower(Sxx, faxis, [sigmaHighpass signmaLowpass], 'psd');

    % Find peak frequency 
    range = faxis>=9.9 & faxis<=15.1;
    Sxx = Sxx(range);
    faxis = faxis(range);   
    [~, idxMax] = max(Sxx);
    spindle.peakFreq(i) = faxis(idxMax);
    clear faxis Sxx range idxMax
end

%% Calculate start time for each spindle:
spindle.timestamp = TimeStamps(spindle.startIdx);   %Originally = (spindle.startIdx - 1)/newFs;

%% Save data to MATLAB file:
EschenkoSpindle = spindle; %#ok<NASGU>
clear spindle
%Request user to name output file:
prompt = {'Enter the filename you want to save it as: (just the name)'};
def = {'Rat#_Day'};
dlgTitle = 'Save .MAT file';
lineNo = 1;
answer = inputdlg(prompt,dlgTitle,lineNo,def);
filename = char(answer(1,:));
matFileName = strcat(filename,'_Spindles.mat');

% Determine if .MAT file already exists:
matFile = fullfile(CSCFilePath, matFileName);
if exist(matFile, 'file') == 2
    save(matFile, 'EschenkoSpindle', 'TimeStamps', 'info', '-append');
else
    save(matFile, 'EschenkoSpindle', 'TimeStamps', 'info');
end
end