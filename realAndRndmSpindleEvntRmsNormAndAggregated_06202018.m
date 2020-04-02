% (rms - min rms)/(max rms - min rms)
%% Select all of the files to include in the analyses and define variables:
[fileSet,pathName] = uigetfile('*.mat','Select file for Place Field Analyses','MultiSelect', 'on');
numFiles = size(fileSet,2); % # of files selected

%% Determine the smallest number of shuffles among the data files:
shuffleNum = zeros(numFiles,1);
for i = 1:numFiles
    shuffleNum(i) = load(fullfile(pathName, fileSet{1,i}), 'numShuffles', '-mat');
end
minNumShuffles = min(shuffleNum);
clear i

normRealRMS = [];
normRandRMS = cell(minNumShuffles, 1);
for i = 1:numFiles
    load(fullfile(pathName, fileSet{1,i}), '-mat');
    % Normalize the real RMS data and concatenate data from all animals:
    maxRMS = max(realData.spikeRMS, [], 2);
    minRMS = min(realData.spikeRMS, [], 2);
    normRMS = (realData.spikeRMS - minRMS)./(maxRMS - minRMS);
    normRealRMS = [normRealRMS; normRMS];
    clear minRMS maxRMS normRMS
    
   % Normalize the randomly RMS data per shuffle and concatenate data from all animals:
   for j = 1:minNumShuffles
       maxRMS = max(rndmData.RMSperShuffle{j, 1}, [], 2);
       minRMS = min(rndmData.RMSperShuffle{j, 1}, [], 2);
       normRMS = (rndmData.RMSperShuffle{j, 1} - minRMS)./(maxRMS - minRMS);
       normRandRMS{j,1} = [normRandRMS{j,1}; normRMS];
   end
end