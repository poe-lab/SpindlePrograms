numSpindles = size(EschenkoSpindle.startIdx,1);
SpindleBins = zeros(numSpindles, 2);
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
    SpindleBins(i,:) = [TimeStamps(start) TimeStamps(stop)];
end

timeIdx = find((SpindleBins(2:end,1) - SpindleBins(1:end-1,2)) <= 0);
newStart = SpindleBins(:,1);
newStart(timeIdx+1) = [];
newStop = SpindleBins(:,2);
newStop(timeIdx) = [];
newTimeBins = [newStart newStop];

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

numEpochs = length(scoredStates);
boutCategory = scoredStates(1,2);
boutStart = scoredStates(1,1);
boutStop = scoredStates(2,1);
c=1;
for i = 2:numEpochs
  if   isequal(scoredStates(i,2), scoredStates(i-1,2))
  else
      c = c+1;
      boutCategory(c) = scoredStates(i,2);
      boutStart(c) = scoredStates(i,1);
  end
  if isequal(i, numEpochs)
      boutStop(c) = scoredStates(i,1) + 10;
  else
      boutStop(c) = scoredStates(i+1,1);
  end
end

boutIdx = boutCategory==2;
boutCategory = boutCategory(boutIdx);
boutStart = boutStart(boutIdx);
boutStop = boutStop(boutIdx);
numBouts = size(boutCategory,2);
nremBins = [];
for i = 1:numBouts
    subIntervalIndx = find(newTimeBins(:,1) >= boutStart(i) & newTimeBins(:,1) <= boutStop(i));
    if ~isempty(subIntervalIndx)
        for j = 1:length(subIntervalIndx)
            %CONTINUE HERE
        end
        clear subIntervalIndx
    end
    
end
