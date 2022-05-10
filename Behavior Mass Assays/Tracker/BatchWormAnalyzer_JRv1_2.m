
function Tracks = BatchWormAnalyzer_JRv1_2(sParam, filename)
%-- This function extracts the major parameters for all tracks in a Tracks file.

global Prefs;
global Tracks;
global Current;

polishRevEvntIncldeFrms = sParam.AnalyzerPolishRevEvntIncldeFrms;
polishRevDirectionChangeThresh = sParam.AnalyzerPolishRevDirectionChangeThresh;
polishRevDirectionChangeUpperThresh = sParam.AnalyzerPolishRevDirectionChangeUpperThresh;
polishRevMaxLength = sParam.AnalyzerPolishRevMaxLength;
polishRevAngularSpeedThresh = sParam.AnalyzerPolishRevAngularSpeedThresh;
LRThresh = sParam.AnalyzerPolishRevLRThresh;
SRThresh = sParam.AnalyzerPolishRevSRThresh;

%-- Initialize Preferences
Prefs = struct();

Prefs.SmoothWinSize = sParam.AnalyzerSmoothWinSize;                  % 3; Size of Window for smoothing track data (in frames)
Prefs.StepSize = sParam.AnalyzerStepSize;                            % 3; Size of step for calculating changes in X and Y coordinates (in frames): 2 for 4M pixels, 3 for 1M pixels
Prefs.SampleRate = sParam.AnalyzerSampleRate;                        % 3; Movie frame rate (frames/sec)
Prefs.TrackFileName = sParam.AnalyzerTrackFileName;                  % ''; Trackfile name

Prefs.TransThresh = sParam.AnalyzerTransThresh;                      % 20; For Makoto's update 09-18-05 - threshold for detecting Omegas
Prefs.RevTransThresh = sParam.AnalyzerRevTransThresh;                % 110; For Makoto's update 09-18-05 - threshold for detecting Reversals
Prefs.RoundThresh = sParam.AnalyzerRoundThresh;                      % 1.47; Cutoff for omega bends
Prefs.EccentricityThresh = sParam.AnalyzerEccentricityThresh;        % 0.75; Cutoff for omega bends

Prefs.FFSpeed = sParam.AnalyzerFFSpeed;                              % 6; Speed of FF (and RW) track playback
Prefs.PixelSize = sParam.AnalyzerPixelSize;                          % 0.0276; 0.0392; new small square setting; old settings 0.0416); %1/(17.64));

Prefs.ROI = ones([2056,2056]);  % defines Arena Size

Prefs.AnalyzerUseDynamicEccentricityThresh = sParam.AnalyzerUseDynamicEccentricityThresh;

% Image calibration - pixels/mm mesured in both vertical and horizontal directions) -
% Microscope in Lowest Magnification setting
% (Leica Mic. cal - 68.5, Navatar Lens Min Zoom cal - 36.5,
% Navatar Lens Exps (N2)401-404 - 49.2)

load (filename);
disp(strcat('... start', 32, mfilename, 32, ', analyzing', 32, filename));
disp(strcat('... (wormanalyzer), size tracks', 32, num2str(size(Tracks,1)), 32, num2str(size(Tracks,2))));

%-- Set up for processing, running in Batch Analysis Mode
Current.BatchAnalysis = 1;

%%%%%%%
%-- Analyze Tracks: basic and old stuff:  size,omegas & angular speed---
NumTracks=length(Tracks);
bi=NaN;
for TN = 1:NumTracks
    
    Tracks(TN).Analyzed = 0;
    Current.Analyzed = 0;
    Current.TempAnalyzed = 0;
    
    %too small?
    if nanmedian(Tracks(TN).FilledArea)<sParam.MinWormArea
        bi(TN)=1;
    end

    OmegaTrackAnalysis_HighMagn_JR1(TN);
    
    % the way angular velocity is calculated is insanely complicated, new version
    % here (deg/frame):
    Tracks(TN).AngSpeed=diff(Tracks(1,TN).Direction);
    
    Tracks(TN).Analyzed = 1;
    
    
    %     end;
end;

Tracks(find(bi))=[];

for currTrack = 1:size(Tracks,2)
    
    if currTrack==7
        currTrack;
    end
    %-- Direction 360 -----------------------------------------------------
    
    %-- the worm analyzer saves the direction of a worm as follows:
    %-- 0 is north or up in a circle.
    %-- to the right are positive angles until +180 degrees
    %-- to the left are negative angles until -180 degrees
    %-- get direction of the worm immediately before the turn and immediately
    %-- after the turn
    %-- change this to 360 degrees and add this to Tracks,
    %-- but remember, 0 is still at north
    try
        
    Tracks(currTrack).Direction360 = Tracks(currTrack).Direction;
    Tracks(currTrack).Direction360(Tracks(currTrack).Direction360<0) = 360 + ...
        Tracks(currTrack).Direction360(Tracks(currTrack).Direction360<0);
    %---------------------------------------------------------------------
    %-- get approximated worm length --------------------------------------
    %----------------------------------------------------------------------
    %-- sort Major axes by size to calculate quartiles
    sortAxes = sort(Tracks(currTrack).MajorAxes);
    %-- get median (quartil 50)
    q50 = median(sortAxes);
    %-- get quartil 75 from median
    q75 = median(sortAxes(find(sortAxes > q50)));
    %-- get mean of upper quartil major axes values to approximate worm length
    approxWormLength = floor(mean(sortAxes(find(sortAxes > q75))));
    Tracks(currTrack).ApproxWormLength = approxWormLength;
    catch ME
    end
    
    if currTrack==17
        1;
    end
    %%%%%%%
    %-- polished reversals ------------------------------------------------
    %----------------------------------------------------------------------
    getRevEvntCndidts = [];
    checkLastEvent = -1;
    
    for i = (1+polishRevEvntIncldeFrms):size(Tracks(currTrack).Direction360,2)-polishRevEvntIncldeFrms
        %-- get events where direction changes drastically indicating a possible reversal
        %-- exclude events, if the direction jump is too large, indicating
        %-- the worm solely moving back and forth between 0 and 360 degrees
        if abs(Tracks(currTrack).Direction360(i-polishRevEvntIncldeFrms) - Tracks(currTrack).Direction360(i)) > polishRevDirectionChangeThresh & abs(Tracks(currTrack).Direction360(i-polishRevEvntIncldeFrms) - Tracks(currTrack).Direction360(i)) < polishRevDirectionChangeUpperThresh
            %-- exclude events, that occur one frame after the last event
            if i-1 > checkLastEvent
                getRevEvntCndidts = [getRevEvntCndidts; Tracks(currTrack).Direction360(i-polishRevEvntIncldeFrms)-Tracks(currTrack).Direction360(i) i];
            end;
            checkLastEvent = i;
        end
    end
    
    finalRevCandidates = [];
    inactiveEvents = [];
    %-- identify all pairs in directional change within specified number
    %-- of frames which are our potential candidates for reversals.
    %=start and end points of reversals
    for currDirChangeIdx = 1:size(getRevEvntCndidts,1)
        clear currNeighbour;
        
        currDirChange = getRevEvntCndidts(currDirChangeIdx,2);
        currIdx = nearestneighbour(currDirChange, getRevEvntCndidts(:,2)', 'Radius', (polishRevMaxLength*Prefs.SampleRate));
        
        %-- keep only neighbours in front of us, we already addressed
        %-- neighbours behind us.
        currIdx = currIdx(find(currIdx > currDirChangeIdx));
        currEvents = getRevEvntCndidts(currIdx, 2);
        
        %-- get all Candidate pairs plus the average of the absolute angular
        %-- speed over the current timeframe
        
        %-- 1) check neighbours within spec number of frames.
        %-- 2) if yes, check mean absolute angular speed above thresh
        %-- 3) if above thresh, take only the very next neighbour and exclude it from
        %-- further analysis as well, because we assume its the ending event of the
        %-- reversal.
        if ~isempty(currEvents)
            for currNeighbour = 1:size(currEvents, 2)
                if isempty(find(inactiveEvents == currDirChange)) & isempty(find(inactiveEvents == currEvents(currNeighbour)))
                    currMeanAngSpeed = nanmean(abs(Tracks(currTrack).AngSpeed(currDirChange+3:currEvents(currNeighbour)-3)));
                       
                    eccentr=1;
                    try
                        eccentr = mean(abs(Tracks(currTrack).Eccentricity(currDirChange-10:currDirChange+10)));
                    catch
                    end
                    if currMeanAngSpeed < polishRevAngularSpeedThresh & eccentr>sParam.AnalyzerEccentricityThresh
                        
                        %-- distinguishing between large and short reversals
                        
                        currRevType = NaN;
                        %normal reversals
                        if  (currEvents(currNeighbour) - currDirChange)/Prefs.SampleRate > SRThresh & (currEvents(currNeighbour) - currDirChange)/Prefs.SampleRate <=LRThresh &...
                                mean (Tracks(1, currTrack).Speed(currDirChange+2:currEvents(currNeighbour)-2))> sParam.RevSpeedThr ...
                                 & max (Tracks(1, currTrack).Speed(currDirChange:currEvents(currNeighbour)))> (sParam.RevSpeedThr)*1.5;
                            currRevType = 1;  %normal reversal
                            
                            inactiveEvents = [inactiveEvents getRevEvntCndidts(find(getRevEvntCndidts(:,2) >=currDirChange & getRevEvntCndidts(:,2) <= currEvents(currNeighbour)),2)'];
                        %long reversals >LRThreshold
                        elseif  (currEvents(currNeighbour) - currDirChange)/Prefs.SampleRate > SRThresh & (currEvents(currNeighbour) - currDirChange)/Prefs.SampleRate >LRThresh &...
                                mean (Tracks(1, currTrack).Speed(currDirChange+2:currEvents(currNeighbour)-2))> sParam.RevSpeedThr ...
                                 & max (Tracks(1, currTrack).Speed(currDirChange:currEvents(currNeighbour)))> (sParam.RevSpeedThr)*1.5
                            currRevType = 2;  %long reversal
                            
                            inactiveEvents = [inactiveEvents getRevEvntCndidts(find(getRevEvntCndidts(:,2) >=currDirChange & getRevEvntCndidts(:,2) <= currEvents(currNeighbour)),2)'];
                            %stops
                        elseif  (currEvents(currNeighbour) - currDirChange)/Prefs.SampleRate > SRThresh & (currEvents(currNeighbour) - currDirChange)/Prefs.SampleRate <=LRThresh &...
                                max (Tracks(1, currTrack).Speed(currDirChange+2:currEvents(currNeighbour)-2))< sParam.RevSpeedThr*1.5 ...
                                & mean (Tracks(1, currTrack).Speed(currDirChange+2:currEvents(currNeighbour)-2)) < sParam.RevSpeedThr
                            
                            currRevType = 3; %slow reversal or brief stop
                            
                            inactiveEvents = [inactiveEvents getRevEvntCndidts(find(getRevEvntCndidts(:,2) >=currDirChange & getRevEvntCndidts(:,2) <= currEvents(currNeighbour)),2)'];
                            
                            %very short reversals
                        elseif  (currEvents(currNeighbour) - currDirChange)/Prefs.SampleRate > SRThresh & (currEvents(currNeighbour) - currDirChange)/Prefs.SampleRate <=LRThresh &...
                                max (Tracks(1, currTrack).Speed(currDirChange+2:currEvents(currNeighbour)-2))> sParam.RevSpeedThr*1.8 ...
                                & mean (Tracks(1, currTrack).Speed(currDirChange+2:currEvents(currNeighbour)-2)) > sParam.RevSpeedThr
                            
                            currRevType = 0;
                            inactiveEvents = [inactiveEvents getRevEvntCndidts(find(getRevEvntCndidts(:,2) >=currDirChange & getRevEvntCndidts(:,2) <= currEvents(currNeighbour)),2)'];
                            
                        end;
                        
                        finalRevCandidates = [finalRevCandidates; currDirChange currEvents(currNeighbour) currMeanAngSpeed currRevType];
                        
                        %check if there was an omega inbetween:
                        try
                            if ~isempty(Tracks(1,currTrack).OmegaTrans) %find omega indices which happened before reversal end
                                oi=find (Tracks(1,currTrack).OmegaTrans(:,1)>currDirChange & Tracks(1, currTrack).OmegaTrans(:,1)<currEvents(currNeighbour));
                                if length(oi)>1
                                    oi=oi(1);
                                end
                                %ignore if omega was >4s before reversal
                                %end
                                if ~isempty(oi) & currEvents(currNeighbour)-Tracks(1, currTrack).OmegaTrans(oi,1)<(4*Prefs.SampleRate)  ...&
                                    currEvents(currNeighbour)-Tracks(1, currTrack).OmegaTrans(oi,1)>10 & Tracks(1, currTrack).Speed(Tracks(1, currTrack).OmegaTrans(oi,1))<0.08;
                                    
                                    if Tracks(1, currTrack).OmegaTrans(oi,1)-currDirChange <LRThresh
                                        currRevType = 1;
                                    end
                                    
                                    finalRevCandidates(end,:)=[];
                                    currMeanAngSpeed = nanmean(abs(Tracks(currTrack).AngSpeed(currDirChange+3:Tracks(1, currTrack).OmegaTrans(oi,1)-3)));
                                    finalRevCandidates = [finalRevCandidates; currDirChange Tracks(1, currTrack).OmegaTrans(oi,1) currMeanAngSpeed currRevType];
                                    inactiveEvents(inactiveEvents>Tracks(1, currTrack).OmegaTrans(oi,1))=[];
                                    inactiveEvents = [inactiveEvents getRevEvntCndidts(find(getRevEvntCndidts(:,2) >=currDirChange & getRevEvntCndidts(:,2) <= Tracks(1, currTrack).OmegaTrans(oi,1)),2)'];
                                    
                                end
                            end
                        catch
                        end
                    else
                        cmas=currMeanAngSpeed;
                    end
                end
            end
        end
    end;
    % end of polished reversals-----
    Tracks(currTrack).polishedReversals = finalRevCandidates;
    %----------------------------------------------------------------------
    
end;

%remove unnecessary fields:
Tracks=rmfield(Tracks,'RingEffect');
Tracks=rmfield(Tracks,'Round');

%-- Save analysis results
[~, savename, ~] = fileparts(filename);
%-- account for file tag - if file tag is not already in the name from
%-- previous steps, insert it before '_tracks'.
if isempty(strfind(savename, sParam.addFileTag))
    savename = strrep(savename, '_tracks', strcat(sParam.addFileTag,'_tracks'));
end;

disp('... save _als.mat');
sParam.FileVersion = mfilename;
sParam.AnalysisDate=datestr(now);

AnalysisParameters = sParam;
%check for filesize:
filesize=vsize(Tracks,'-r','-bs','mb');
if filesize{2}>1500
    disp('tracks file>1 GB');
    save([savename '_JR_als.mat'], 'AnalysisParameters', 'Tracks','-v7.3'); % we introduced the v7-3 switch for saving large files
else
    save([savename '_JR_als.mat'], 'AnalysisParameters', 'Tracks');
end

%-- Return to Non-Batch Mode
Current.BatchAnalysis = 0;