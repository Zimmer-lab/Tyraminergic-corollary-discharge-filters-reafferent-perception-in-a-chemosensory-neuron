%-- remove clear all to use breakpoints!

%Julias version, just tracking and analyzing, most old&uneccessary stuff
%removed

%-- IMPORTANT! if you want to add new parameters in the parameter file,
%-- don't forget to include them in the wormtracker keyList as well!
%-- Otherwise they won't be processed...

function Batch_track_compress_copy_JR_v1(dirParameterFile)
warning off
display(['...using' dirParameterFile])
%-- reset workspace
close all;
clearvars -except dirParameterFile;

%-- variable for error and intensity report at the end of the logfile
logReport = '';

%-- catch errormessages
try
    %-- dock figures inside matlab window in none cluster mode
    %set(0, 'DefaultFigureWindowStyle', 'docked');
    set(0, 'DefaultFigureWindowStyle', 'normal');
    
    %-- suppress annoying warning messages to increase logfile readability
    %-- supress message 'Image is too big to fit on screen; displaying at 33%'
    warning('off', 'images:initSize:adjustingMag');
    %-- supress message 'The '-dill' print device will be removed in a future release. Use Encapsulated PostScript ('-depsc') instead.'
    warning('off', 'MATLAB:print:Illustrator:DeprecatedDevice');
    
    %import Useful_functions_MS_v104.copyfiles;
    %import Useful_functions_MS_v104.all2str;
    
    %-- USER INPUT ------------------------------------------------------------
    %--------------------------------------------------------------------------
    
    keyList = getKeyList_JR1('wormtracker');
    
    %-- retrieve input parameter structure from parameter file
    %-- comment out in compiled version
    sParam = parseParameterFile(dirParameterFile, keyList);
    %-- comment in in compiled version
    %sParam = parseParameterFile(dirParameterFile, keyList);
    
    %--------------------------------------------------------------------------
    %-- USER INPUT END --------------------------------------------------------
    %-- set main figure visibility on or off
    set(0, 'DefaultFigureVisible', sParam.vDisplayPlots);
    
    %-- get local working directory
    %-- needed on cluster to copy movie to cluster node that
    %-- is doing the actual computation to reduce reading time for
    %-- each frame when tracking movies.
    sParam.localWorkingDirectory = pwd;
    sParam.ParameterFile=dirParameterFile;
    
    %-- change to working directory containing movies
    if sParam.workingDirectory=='current'
        display('...current directory')
    else
    cd(sParam.workingDirectory);
    end
    
    %-- print messages to logfile; logfile will be saved in the current working
    %-- directory
    logfile = sprintf('%slogfile.txt',datestr(now,'yyyymmdd_HHMMSS'));
    diary(logfile);
    
    disp('...');
    disp(strcat('... you are running', 32, mfilename));
    
    
    %-- track worms from movies-----------------------------
    
    try
        if sParam.vUseWormTracker == 1
            disp(strcat('... local working directory', 32, sParam.localWorkingDirectory));
            [logReport, Tracks] = BatchWormTrackerHist_JRv1_largeFiles(sParam, logReport);
        else
            disp('tracker module = off')
        end;
    catch err
        disp(getReport(err));
        logReport = strcat(logReport, 10, getReport(err));
    end;
    
    %% -- analyze tracks------------------------------------
    try
        if sParam.vUseWormAnalyzer == 1
            disp('...analyzing tracks..');
            trackFiles = dir(sParam.TracksMatFileSearchString)
            
            for currentTrackFile = 1 : size(trackFiles)
                
                trackFileName = (trackFiles(currentTrackFile).name);
                
                analyzed_name= [trackFileName(1:end-4)  '_JR_als.mat'];
                %---check if the als file already exists:
                if exist(analyzed_name)~=2
                    display(currentTrackFile)
                    Tracks = BatchWormAnalyzer_JRv1_2(sParam, trackFileName);
                    Tracks = [];
                    
                else display ([analyzed_name,'...already exists']);
                end
            end;
        end;
    catch err
        disp(getReport(err));
        logReport = strcat(logReport, 10, getReport(err));
    end;
    set(0, 'DefaultFigureVisible', 'on');
    
catch err
    disp(getReport(err));
    logReport = strcat(logReport, 10, getReport(err));
end;

disp('...');
disp('... Final Logging Report');
disp('...');
disp(logReport);
disp('...');
diary off;

close all;

end