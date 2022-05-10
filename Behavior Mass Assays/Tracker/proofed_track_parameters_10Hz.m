%% for 10 hz!!!--------------------------------------------------------------------------
%-- for the general use of this parameter file:
%-- 1) keywords not in file getKeyList.m will be ignored
%-- 2) this parameter file is actually a normal textfile, but everything 
%--         in here will be treated like matlab code. therefore using comments works.
%--------------------------------------------------------------------------

%-- define if you run the script on the cluster or on a local machine
%-- this affects whether movies to be tracked are temporarily copied from 
%-- the working directory to the local temporary folder to speed up 
%-- processing on the cluster.
%-- useCluster, 1 ... run on 15cluster, 0 ... run locally
useCluster = 0;

%-- working directory containing movies to be tracked or files to be
%-- processed
%-- e.g. workingDirectory = 'D:\Sonntag_Michael\testingGround';
workingDirectory = 'current';
%-- ACTIVATE TRACKER MODULE -----------------------------------------------
%-- use worm tracker, 1 ... yes, 0 ... no
%-- requires *.avi files in the workingDirectory folder
vUseWormTracker = 1;

%-- if you want all movies in the current folder analyzed, use '*.avi'
%-- if you want only a specific movie analyzed, paste the whole name
%-- including file extension below. Then only this file will be found.
%-- e.g. sParam.MovieSearchString = 'Onlythismovie.avi';
%-- if you want only specific movies analysed, which can be distinguished
%-- by filename, you can use windows search patterns. 
%-- e.g. sParam.MovieSearchString = '2013*eat4*.avi'; retrieves all movies
%-- starting with 2013, contain eat4 and end with .avi.
%-- you can always try a search pattern in windows explorer before using it
%-- here.
MovieSearchString = '*jpeg*.avi';

%-- Movie codec BUG WORKAROUND --------------------------------------------
%-- if you use the raw, uncompressed avi recordings, use 
%-- setting MovieTypeUsed = 'AVI'
%-- if you use mp4 or already compressed avi's, use MovieTypeUsed = 'MP4'
MovieTypeUsed = 'AVI';

%--------------------------------------------------------------------------

%-- ACTIVATE ANALYSER MODULE ----------------------------------------------
%-- use worm analyzer, 1 ... yes, 0 ... no
%-- requires *_track.m files in the workingDirectory folder
%-- will be automatically created, if vUseWormTracker = 1
vUseWormAnalyzer = 1;

%-- FIND FILES-------------------------------------------------------------
%-- search pattern for tracks files. Only files included in the search
%-- string will be analyzed. 
%-- default: singleTracksMatSearchString = '*_tracks.mat';
TracksMatFileSearchString = '*_tracks.mat';
%--------------------------------------------------------------------------

%-- search pattern for als files. Only files included in the search
%-- string will be analyzed. 
%-- default: AlsMatFileSearchString = '*_als.mat';
AlsMatFileSearchString = '*_als.mat';

%-- if you want to average data, specify filename for the average output files
%-- e.g. avgSaveFilename = 'N2_shift_21_04_21';
%-- current date will be added automatically to the beginning of the
%-- average filename
avgSaveFilename = 'N2_s_21_10_21_avg';
%--------------------------------------------------------------------------

%-- GENERAL SETTINGS ------------------------------------------------------
%-- decide whether plots should be displayed within the matlab session
%-- ... use 'on' or 'off' accordingly. Use 'off' ALWAYS if you are working
%-- on the cluster.
vDisplayPlots = 'off';

%-- with this you can mark ALL created files with a tag, e.g. if you try
%-- different parameters with the same dataset to distinguish between
%-- the individual files. The tag will be added before the last tags like
%-- '_tracks.mat', '_als.mat' or '_dats.mat'. if you dont need this, use
%-- default: addFileTag = '';
addFileTag = '';

%-- decide whether altogether plots should be saved or not
%-- 1 ... yes, 0 ... no
checkSave = 1;

%-- supported file types for saving figures:
%-- 'fig', 'pdf', 'eps', 'ai'
%-- define, in which file types plots should be saved
saveFileTypes = { 'pdf'};

%-- decide, whether single plots should be saved (single plots are plots of 
%-- e.g. omega turns only or reversals only)
%-- switching useSinglePlots off still creates and saves the multiplot figures
%-- 1 ... yes, 0 ... no
useSinglePlots = 0;
%-- END GENERAL SETTINGS --------------------------------------------------


%-- Tracking and ANALYZING PARAMETERS -------------------------------

%-- frames per s
SampleRate = 10;
%-- tickmarks along the axes
tickmarks = 180;

%-- WORM TRACKER PARAMETERS -----------------------------------------------
%-- Min Length of valid track (in frames)
MinTrackLength = 10*SampleRate;
%-- Min area for object to be a valid worm
MinWormArea = 80; %LR:40; HR 100
%-- Max area for object to be a valid worm
MaxWormArea = 600; %LR 120;% HR 500
%-- Max wormsize change between frames (in pixels)
SizeChangeThreshold = 100;   %% adapted to HR

BackgroundLevel = 0.3;

%-- IMPORTANT SETTING ---------------------
%-- Level for identifying shapes of individual worms when converting to 
%-- binary image, dependent on movie intensity levels. Dynamic level uses
%-- an automated thresholding(recomended).
%-- Read detailed description below, if worms are not recognised properly.
%-- default setting: BWLevel = 0.2;
BWLevel = 0.2;
DynamicLevel = 1;

%-- Display tracking results every 'PlotFrameRate' frames
PlotFrameRate = 100;
%-- binning for worm size histogram
HistoTrackSize = 100; %-- 50 for high magnification

%-- if the following option is set to 1, the background used for analysing
%-- the movie will be created dynamically using function getbackgroundMZ64.
%-- if the following is set to 0, then the current folder has to contain 
%-- a tiff image named 'background.tiff', which will be used instead.
%-- default setting: TrackerUseDynamicBackground = 1;
TrackerUseDynamicBackground = 1;

%-- END WORM TRACKER PARAMETERS -------------------------------------------

%-- WORM ANALYZER JR PARAMETERS ----------------------------------------
%-- default frame rate ... 3
%-- using a different frame rate will affect the analysis of reverse omegas, 
%-- omega and shallow turns a factor (floor(framerate/3)) will be used to 
%-- account for the higher frame rate.
AnalyzerFrameRate = 10;
AnalyzerStepSize = 1;                             % Size of step for calculating changes in X and Y coordinates (in frames): 2 for 4M pixels, 3 for 1M pixels
AnalyzerSampleRate = 10;                           % Movie frame rate (frames/sec)
AnalyzerTrackFileName = '';                       % Trackfile name
AnalyzerPixelSize = 0.0155;                       % mm/pixel...0.025 for stage bheavior ; 0.155 for HR camera setting with maxi
AnalyzerSmoothWinSize = 2;                        %2 for 10 hz; Size of Window for smoothing track data (in seconds?)
%--Omegas
AnalyzerTransThresh = 20;                         % For Makoto's update 09-18-05 - threshold for detecting Omegas
AnalyzerRevTransThresh = 110;                     % For Makoto's update 09-18-05 - threshold for detecting Reversals
AnalyzerRoundThresh = 1.125;                      % Cutoff for omega bends
AnalyzerEccentricityThresh = 0.88;                % Cutoff for omega bends

AnalyzerFFSpeed = 6;                              % Speed of FF (and RW) track playback

%-- use dynamic eccentricity threshold to identify omega turns
%-- if this is set to 0 the value of parameter AnalyzerEccentricityThresh
%-- will be used with every track
%-- if this is set to 1, then the maximal value of either the parameter
%-- AnalyzerEccentricityThresh or a dynamically calculated value will be used
%-- to account for tracks displaying large differences in mean eccentricity
AnalyzerUseDynamicEccentricityThresh = 0;

%-- this parameter defines the threshold for Omega turns ----------
%-- identified turns above this threshold will be counted as omega turns
AnalyzerOmegaTurnAngSpeedThresh = 45;
%------------------------------------------------------------------

%-- parameters to identify shallow turns ------
AnalyzerShallowTurnEccThresh = 0.96;
AnalyzerShallowTurnAngSpeedTresh = 60;
AnalyzerShallowTurnDirecThresh = 150;
%----------------------------------------------

%-- parameters to identify reversals ----------
RevSpeedThr=0.02;
AnalyzerPolishRevEvntIncldeFrms = 1;
AnalyzerPolishRevDirectionChangeThresh = 50;
AnalyzerPolishRevDirectionChangeUpperThresh = 340;
AnalyzerPolishRevMaxLength = 15; % in seconds; 
AnalyzerPolishRevAngularSpeedThresh = 20; % should be low between reversal start and end, set to relatively high value 20 to include reversals within omegas.
% threshold for long reversals: min length in sec
AnalyzerPolishRevLRThresh = 8.6;

% short reversals
AnalyzerPolishRevSRThresh = 0.5;   % changed from previously 1.5

%-- set time window for identifying reverse omegas
%-- maximum number of frames from the end of a reversal to the
%-- beginning of an omega turn to be still counted as a reverse omega
AnalyzerRevOmegaTimeWindow = 11;

%--------------------------------------------------------------------------
%%%removed parameters%%%

