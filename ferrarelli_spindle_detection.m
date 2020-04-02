function ferrarelli_spindle_detection(batchProcess)
% FERRARELLI Detect sleep spindles using the Ferrarelli algorithm.
% Ferrarelli et al. "Reduced Sleep Spindle Activity in Schizophrenia
% Patients", Am J Psychiatry 164, 2007, pp 483-492
%
% Input is recorded EEG for a whole night of sleep (any channel can be
% used), the sampling frequency, a stage file (STA) that has been loaded in
% to MATLAB.
% Output is a binary vector containing ones at spindle samples.
% Syntax: detection = ferrarelli_spindle_detection(C3,fs,stage_file)
%
% Adopted from Ferrarelli by Sabrina Lyngbye Wendt, July 2013

%% Set constants:
Fs = 250; %250 samples/second. Replace with variable that extracts from header if needed    
lower_thresh_ratio = 2;
upper_thresh_ratio = 8;
epochsize = 30;
%% Select folder with EDF files:
working_dir=pwd;
if batchProcess
    % Select folder to save MATLAB output:
    fileSelectedCheck = 0;
    while isequal(fileSelectedCheck,0)
        matDataFolder = uigetdir;
        if isempty(matDataFolder)
            uiwait(errordlg('You need to select a folder. Please try again',...
                'ERROR','modal'));
        else
            fileSelectedCheck = 1;
        end 
    end
    
    % Select folder and get list of EDF files:
    fileType = '*.edf';
    [dataFolder, fileList, numberOfDataFiles] = batchLoadFiles(fileType);
    
    % Select folder containing the text version of the scored files:
    fileType = '*.txt';
    [annotationFolder, ~, ~] = batchLoadFiles(fileType);
    
else
    numberOfDataFiles = 1;
    dataFolder = [];
    fileName = [];
    matDataFolder = [];
    fileSelectedCheck = 0;
    
    % Select a single EDF file:
    while isequal(fileSelectedCheck,0)
        [fileName, dataFolder] = uigetfile('*.edf', 'Select the EDF data file');
        if isempty(fileName) || isempty(dataFolder)
            uiwait(errordlg('You need to select a file. Please try again',...
                'ERROR','modal'));
        else
            fileSelectedCheck = 1;
        end 
    end
    cd(working_dir);
    
    fileSelectedCheck = 0;
    % Select a single annotation file:
    while isequal(fileSelectedCheck,0)
        [annotationFileName, annotationFolder] = uigetfile('*.txt', 'Select the scored text file');
        if isempty(annotationFileName) || isempty(annotationFolder)
            uiwait(errordlg('You need to select a file. Please try again',...
                'ERROR','modal'));
        else
            fileSelectedCheck = 1;
        end 
    end
    
    % Select folder to save MATLAB output:
    fileSelectedCheck = 0;
    while isequal(fileSelectedCheck,0)
        matDataFolder = uigetdir;
        if isempty(matDataFolder)
            uiwait(errordlg('You need to select a folder. Please try again',...
                'ERROR','modal'));
        else
            fileSelectedCheck = 1;
        end 
    end
    cd(working_dir);
    
end

%% BANDPASS_FILTER_FERRARELLI Bandpass filter used by Ferrarelli et al.
% This function creates a 12th order (if the sampling frequency is 100 Hz)
% Chebyshev Type II bandpass filter with passband between 10 and 16 Hz. The
% filter is -3 dB at 10.7 and 15 Hz.
% The input signal is filtered with the created filter and the filtered
% signal is returned as output.
Wp=[11 15]/(Fs/2);
Ws=[10 16]/(Fs/2);
Rp=3;
Rs=40;
[n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bbp, abp]=cheby2(n,Rs,Wn);

%% Filter used for spindle activity analyses for MDD v HC paper
sigmaHighpass = 11;
signmaLowpass = 16;
[z, p, k] = ellip(7,1,60, [sigmaHighpass signmaLowpass]/(Fs/2),'bandpass');
[sos, g] = zp2sos(z,p,k);

%% Batch process spindle detection:
for m = 1:numberOfDataFiles
    if batchProcess
        fileName = strtrim(fileList(m,:)); %Removes any white space at end of file name string.
        matFileName = strrep(fileName, '_edited.edf','.mat'); 
        annotationFileName = strrep(fileName, '_edited.edf', '_edited_annotations.txt');
    else
        matFileName = fileName(1:8);
    end
    
    edfFile = fullfile(dataFolder,fileName);
    [~, ~, signalCell] = blockEdfLoad(edfFile);
    signal = (signalCell{1,5});

    clear signalCell edfFile fileName
  
    %% Load scored stages:    
    annotationFile = fullfile(annotationFolder,annotationFileName);
    delimiter = ',';
    startRow = 2;
    clear annotationFileName

    % Format string for each line of text:
    %   column1: text (%s)
    %	column2: double (%f)
    %   column3: text (%s)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%s%f%s%[^\n\r]';

    % Open the text file.
    fileID = fopen(annotationFile,'r');
    if fileID == -1
        continue
    end
        

    % Read columns of data according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

    % Close the text file.
    fclose(fileID);
    
    %Convert the contents of column with dates to vector of date components:
    dataArray{1} = datevec(dataArray{1}, 'HH:MM:SS.FFF');

    % Allocate imported array to column variable names
    Onset = dataArray{:, 1};
    Annotation = dataArray{:, 3};

    % Clear temporary variables
    clearvars annotationFileName annotationFile delimiter startRow formatSpec fileID dataArray ans
    
    %Convert Onset to seconds from start of recording:
    scoredTimestamp = Onset(:,4)*3600 + Onset(:,5)*60 + Onset(:,6);
    clear Onset
    addStagesToStart = zeros((scoredTimestamp(1)/epochsize),1); %add to unscored gap at beginning of recording
    clear scoredTimestamp
    
    %Convert scored stages to strings:
    C = char(Annotation);
    clear Annotation
    stagesString = C(:, 1);
    clear C
    
    %Convert stages to numbers:
    numStages = length(stagesString);
    sleep = zeros(numStages,1);
    for j = 1:numStages
        switch stagesString(j)
            case '1'
                sleep(j) = 1;
            case '2'
                sleep(j) = 2;
            case '3'
                sleep(j) = 3;
            case '4'
                sleep(j) = 4;
            otherwise
                sleep(j) = 0;
        end
    end
    sleep = [addStagesToStart; sleep]; %#ok<AGROW>
    clear stagesString numStages addStagesToStart

    %% Redefine sleep stage numbers 
    sleep(sleep==1)=-1;sleep(sleep==2)=-2;sleep(sleep==3)=-3;sleep(sleep==4)=-4;
    sleepsamp = reshape(repmat(sleep,1,epochsize*Fs)',1,length(sleep)*epochsize*Fs)';
    
    %% Resize EEG signal and sleep scoring vectors to be same length in case
    % file not scored to end of recording or the last epoch was short:
    if length(sleepsamp) < length(signal)
        signal = signal(1:length(sleepsamp));
    else
        sleepsamp = sleepsamp(1:length(signal));
    end
    
    %% Get the desired nrem samples
    nremsamples = find(sleepsamp<=-2); % Use data from stage S2+S3+S4
    clear sleep

    %% Bandpass filter from 11-15 Hz and rectify filtered signal
    BandFilteredData = filtfilt(bbp, abp, signal);
    RectifiedData = abs(BandFilteredData);
    clear BandFilteredData

    %% Create envelope from the peaks of rectified signal (peaks found using zero-crossing of the derivative)
    datader = diff(RectifiedData); % x(2)-x(1), x(3)-x(2), ... + at increase, - at decrease
    posder = zeros(length(datader),1);
    posder(datader>0) = 1; % index of all points at which the rectified signal is increasing in amplitude
    diffder = diff(posder); % -1 going from increase to decrease, 1 going from decrease to increase, 0 no change
    envelope_samples = find(diffder==-1)+1; % peak index of rectified signal
    Envelope = RectifiedData(envelope_samples); % peak amplitude of rectified signal

    %% Finds peaks of the envelope
    datader = diff(Envelope);
    posder = zeros(length(datader),1);
    posder(datader>0) = 1; % index of all points at which the rectified signal is increasing in amplitude
    diffder = diff(posder);
    envelope_peaks = envelope_samples(find(diffder==-1)+1); % peak index of Envelope signal
    envelope_peaks_amp = RectifiedData(envelope_peaks); % peak amplitude of Envelope signal
    clear datader posder Envelope

    %% Finds troughs of the envelope
    envelope_troughs = envelope_samples(find(diffder==1)+1); % trough index of Envelope signal
    envelope_troughs_amp = RectifiedData(envelope_troughs); % peak trough of Envelope signal
    clear diffder envelope_samples

    %% Determine upper and lower thresholds
    nrem_peaks_index=sleepsamp(envelope_peaks)<=-2; % extract samples that are in NREM stage S2+S3+S4
    [counts, amps] = hist(envelope_peaks_amp(nrem_peaks_index),120); % divide the distribution peaks of the Envelope signal in 120 bins
    [~,maxi] = max(counts); % select the most numerous bin
    ampdist_max = amps(maxi); % peak of the amplitude distribution
    lower_threshold = lower_thresh_ratio*ampdist_max;
    upper_threshold = upper_thresh_ratio*mean(RectifiedData(nremsamples)); %#ok<FNDSB>
    clear nremsamples nrem_peaks_index counts amps maxi ampdist_max

    %% Find where peaks are higher/lower than threshold
    below_troughs = envelope_troughs(envelope_troughs_amp<lower_threshold); % lower threshold corresponding to 4* the power of the most numerous bin
    %above_peaks=envelope_peaks(envelope_peaks_amp>upper_threshold & sleepsamp(envelope_peaks)<=-2); % Use this line insted of next if spindles should only be detected in S2+S3+S4
    above_peaks = envelope_peaks(envelope_peaks_amp>upper_threshold);
    clear envelope_troughs envelope_peaks envelope_peaks_amp envelope_troughs...
        envelope_troughs_amp lower_threshold upper_threshold
    %% For each of peaks above threshold
    spistart = NaN(length(above_peaks),1); % start of spindle (in 100Hz samples)
    spiend = NaN(length(above_peaks),1); % end of spindle (in 100Hz samples)

    nspi=0; % spindle count
    % for all indexes of peaks (peaks of peaks)
    i = 1;
    while i <= length(above_peaks)
        current_peak = above_peaks(i);
        % find troughs before and after current peak
        trough_before = below_troughs(find(below_troughs > 1 & below_troughs < current_peak,1,'last'));
        trough_after  = below_troughs(find(below_troughs < length(RectifiedData) & below_troughs > current_peak,1,'first'));

        if ~isempty(trough_before) && ~isempty(trough_after)  % only count spindle if it has a start and end
            nspi=nspi+1;
            spistart(nspi)=trough_before;
            spiend(nspi)=trough_after;
            % if there are multiple peaks, pick the highest and skip the rest
            potential_peaks = above_peaks(above_peaks > trough_before & above_peaks < trough_after);
    %         [~, maxpki]=max(RectifiedData(potential_peaks));
    %         current_peak=potential_peaks(maxpki);

            i = i+length(potential_peaks); % adjust the index to account for different max
        else
            i = i+1;
        end
    end
    clear RectifiedData below_troughs above_peaks

    %% Determine scored state for each spindle
    scoring = NaN(nspi,1);  % sleep state where end lies
    for j=1:nspi
        scoring(j) = sleepsamp(spiend(j));
    end
    clear sleepsamp

    %Remove all unwanted sleep states:
    scoring(scoring<-2 | scoring>-2) = NaN;
    spistart = spistart(isnan(scoring)~=1);
    spiend = spiend(isnan(scoring)~=1);
    scoring = scoring(isnan(scoring)~=1);

    %% Exclude spindles based on duration
    duration = spiend-spistart;
    duration(duration<0.3*Fs | duration>3*Fs) = NaN;
    spindle.startIdx = spistart(isnan(duration)~=1);
    spindle.stopIdx = spiend(isnan(duration)~=1);
    spindle.scoring = scoring(isnan(duration)~=1);
    spindle.duration = duration(isnan(duration)~=1)/250;
    clear duration scoring spistart spiend

    %% Create the binary output vector
    % detection = zeros(size(signal));
    % spistart = spistart(isnan(spistart)~=1);
    % spiend = spiend(isnan(spiend)~=1);
    % for k = 1:length(spistart)
    %     detection(spistart(k):spiend(k)) = 1;
    % end

    %% Calculate Spindle Characteristics
    filtSignal = filtfilt(sos,g, signal);
    clear signal

    for i = 1:length(spindle.duration)
        x = filtSignal(spindle.startIdx(i):spindle.stopIdx(i)); % pull out spindle from signal

        % Find the max peak-to-peak amplitude
        [~, absIdx] = findpeaks(abs(x)); % finds all of the peaks and troughs
        xPeaks = x(absIdx); % Gets values of peaks and troughs
        p2pAmp = abs(diff(xPeaks)); % finds all peak-to-peak amplitudes
        [spindle.maxP2pAmp(i), p2pMaxIdx] = max(p2pAmp); 
        spindle.symmetry(i) = absIdx(p2pMaxIdx)/length(x);
        clear absIdx xPeaks p2pAmp p2pMaxIdx

        % Compute the power spectrum of the Hann tapered data:
        dt = 1/Fs; % Define the sampling interval.
        df = 1/spindle.duration(i); % Determine the frequency resolution.
        fNyquist = Fs/2; % Determine the Nyquist frequency.
        faxis = (0:df:fNyquist); % Construct the frequency axis.
        xh = hann(length(x)).*x;
        Sxx = 2*dt^2/spindle.duration(i) * fft(xh).*conj(fft(xh)); % Compute the power spectrum of Hann tapered data.
        Sxx = real(Sxx); %Ignores imaginary component.
        Sxx = Sxx(1:length(faxis));
        clear dt df fNyquist xh

        % Calculate bandpower
        spindle.power(i) = bandpower(Sxx, faxis, [11 16], 'psd');

        % Find peak frequency 
        range = faxis>=10 & faxis<=17;
        Sxx = Sxx(range);
        faxis = faxis(range);   
        [~, idxMax] = max(Sxx);
        spindle.peakFreq(i) = faxis(idxMax);
        clear faxis Sxx range idxMax
    end

    %% Calculate start time for each spindle:
    spindle.timestamp = (spindle.startIdx - 1)/Fs;

    %% Save data to MATLAB file:
    ferrarelliSpindle = spindle; %#ok<NASGU>
    clear spindle
    % Determine if .MAT file already exists:
    matFile = fullfile(matDataFolder,matFileName);
    if exist(matFile, 'file') == 2
        save(fullfile(matDataFolder,matFileName), 'ferrarelliSpindle', '-append');
    else
        save(fullfile(matDataFolder,matFileName), 'ferrarelliSpindle');
    end
end
end