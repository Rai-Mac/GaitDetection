function [WalkRestDPLM2G, WalkBoutsLM2G, VertAxis] = WBDetect_Generic(AcData, Fs, OKPause, LFi)
%Example code of identifying walking independent of oreintation
%Author: Mhairi K MacLean
%Comments: This is the long handed version of the code. It is written in
%this long format to assist with understanding the methods used.
%Should be rewritten for speed.
if nargin <3 %Characteristics of Very very low pass filter. If not given as inpput, define here
    LFi = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 0.25, 'SampleRate', Fs);
end
MaxWalkThresh = 0.4;
LPD = 0.025;

LinearAcc = filtfilt(LFi, AcData); %filter all 3 axes
MinusAcc = AcData - LinearAcc; %Subtract filtered data to adjust acceleration
Magnitude = sqrt(MinusAcc(:,1).^2+MinusAcc(:,2).^2+MinusAcc(:,3).^2); %find the magnitude of activity accross axes


TestThresholdMin = 0.05; %If the median is too high, replace with a constant threshold

MagnitudeThresh = Magnitude>TestThresholdMin; %1 if above walking threshold, 0 if below walking threshold
SMTG = smoothdata(MagnitudeThresh, 'gaussian');

%Identifying starts and ends of walking bouts
MTDiff = diff(MagnitudeThresh); %0 = no change, 1 = start walking, -1 = end walking
WalkEnd = find(MTDiff == -1)+1; %+1 to account for change in index due to the diff function
WalkStart = find(MTDiff == 1)+1; %index 1 of diff function means change occured in original data at dp 2

if WalkEnd(1) < WalkStart(1) %If the first end is before the first start
    WalkStart = [1; WalkStart]; %make the first start be the first dp
end
if WalkEnd(end) < WalkStart(end) %If the last start is after the last end
    WalkEnd = [WalkEnd; length(MTDiff)]; %make the last end by the final dp
end

%To make "orientation independent", look at the low pass filter of all
%axes, and find the one closest to -1G
%First, ttranslate bout starts and stops to logical
WalkRestDPEarly = zeros(size(Magnitude)); %Initialize WalkRest variable as all 0, same length as data
WalkIndiciesEarly = arrayfun(@(s,f)s:f,WalkStart,WalkEnd,'UniformOutput',false); %using array fun to compile indicies of each walking bout
WalkRestDPEarly([WalkIndiciesEarly{:}]) = 1; %Every element in the bouts is set to equal 1.

MedFiltAcc = median(LinearAcc(logical(WalkRestDPEarly),:));
%Is it upside or correct side up? Just in case, do the absolute.
[~,VertAxis]=min(abs(abs(MedFiltAcc)-1));
UpOrDownCorrection = -1 *  sign(MedFiltAcc(VertAxis));
%If x axis average is not within 0.8 and 1.2 g, then exclude 
% for BN = 1:length(WalkStart)
%     AvV(BN,1) = UpOrDownCorrection * mean(AcData(WalkStart(BN):WalkEnd(BN), VertAxis));
% end
% 
%     StraightB = find((AvV > (-1-GR)) & (AvV < (-1+GR)));
%     WalkStartS = WalkStart(StraightB);
%     WalkEndS = WalkEnd(StraightB);
WalkStartS = WalkStart; %NO VERT IS VERY SIMILAR. MIGHT AS WELL USE IT
WalkEndS = WalkEnd;

GapIndex = [WalkEndS(1:end-1), WalkStartS(2:end) ];
SmThr = 0.2; % Smoothed data compared to this threshold. Anything below this considered a pause in walking.
for GapN = 1:length(GapIndex)
    %Looks at the smoothed data in every gap. Any element less than SmThr is
    %marked as 1, then elements are summed. 0 means the smoothed data
    %does NOT reach 0 during the gap.
    GapIndex(GapN,3) = sum(SMTG([GapIndex(GapN,1): GapIndex(GapN,2)]) < SmThr);
end
%OKPause = 20; %Anything greater this is considered broken walking (i.e. there is a gap)
OKGapsG = GapIndex(:,3)>OKPause;
WalkEndNoGaps2G = [WalkEndS(OKGapsG); WalkEndS(end)]; %Extract all ends where the gap did not come below the threshold
WalkStartNoGaps2G = WalkStartS([true; OKGapsG]); %Extract all starts where the preceeding gap did not drop below the threshold

%Add section that rejects any of the data above ThresholdMax
%TestThresholdMax = 0.4;
LBN = [];

SmMag = smoothdata(Magnitude, 'gaussian');
for BN = 1:length(WalkStartNoGaps2G)
    DT=SmMag(WalkStartNoGaps2G(BN):WalkEndNoGaps2G(BN)); %The magnitude of the identified activity bout
    
    DTThresh = DT>MaxWalkThresh; %1 if above the max threshold
    LP = round(length(DT)*LPD);
    if (sum(DTThresh)) < LP %If above the max threshold for less than 5% of the bout
        %Include this walking bout
        LBN = [LBN, BN];
    end
end
WalkStartNoGapsLessMax2G = WalkStartNoGaps2G(LBN);
WalkEndNoGapsLessMax2G = WalkEndNoGaps2G(LBN);


%If walking bout is less than 2 seconds, then remove walking
BoutDurationLM2G = WalkEndNoGapsLessMax2G - WalkStartNoGapsLessMax2G; %Duration of walking bout
OKBoutLength = 2*Fs; %Value of minimum acceptable bout length. (same as gap length, but can be changed)
LongBoutsLM2G = BoutDurationLM2G>OKBoutLength; %1 if bout is long enough, 0 if bout is too short
WalkStartNoGapsLessMaxLongBouts2G = WalkStartNoGapsLessMax2G(LongBoutsLM2G); %Extract all starts where bout is over threshold
WalkEndNoGapsLessMaxLongBouts2G = WalkEndNoGapsLessMax2G(LongBoutsLM2G); %Extract all ends where bout is over threshold
WalkBoutsLM2G = [WalkStartNoGapsLessMaxLongBouts2G WalkEndNoGapsLessMaxLongBouts2G];

%All datapoints, 1=walking, 0=no walking, with small gaps and short bouts

WalkRestDPLM2G = zeros(size(Magnitude)); %Initialize WalkRest variable as all 0, same length as data
WalkIndicies2G = arrayfun(@(s,f)s:f,WalkStartNoGapsLessMaxLongBouts2G,WalkEndNoGapsLessMaxLongBouts2G,'UniformOutput',false); %using array fun to compile indicies of each walking bout
WalkRestDPLM2G([WalkIndicies2G{:}]) = 1; %Every element in the bouts is set to equal 1.

LogicalWalkRestDP2G = logical(WalkRestDPLM2G);
end
%%
