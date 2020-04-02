function ReviewOnlySpindles_03262018
%% Select .mat file:
[spindleFilename, spindleFilePath] = uigetfile({'*.mat',...
        'Pick spindle data files.'},'Select .MAT File');
spindleFile = fullfile(spindleFilePath, spindleFilename);
load(spindleFile, '-mat')

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
spindleSignal = filtfilt(sos, g, signal);

%% Calculate and plot results based on each detected spindle centered in a 3 sec window:
numSpindles = size(EschenkoSpindle.startIdx,1); %Determine number of detected spindles on the selected LFP/EEG channel
% Determine maximum y-axis values (mV) for the spindle filtered signals plus
% added 10%
maxY_Spindle = max(abs(spindleSignal)) + 0.1*max(abs(spindleSignal));

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
      
    ax1 = subplot(1,1,1);
    plot(ax1, TimeStamps(start:stop), spindleSignal(start:stop), 'b')
    hold on
    plot(ax1, TimeStamps(EschenkoSpindle.startIdx(i):EschenkoSpindle.stopIdx(i)),...
        spindleSignal(EschenkoSpindle.startIdx(i):EschenkoSpindle.stopIdx(i)), 'g')
    plot(ax1, TimeStamps(start:stop), signal(start:stop), 'Color', [0.5 0.5 0.5])

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
    ylim([-maxY_Spindle maxY_Spindle])
    ylabel(ax1, 'LFP')
    pause
end