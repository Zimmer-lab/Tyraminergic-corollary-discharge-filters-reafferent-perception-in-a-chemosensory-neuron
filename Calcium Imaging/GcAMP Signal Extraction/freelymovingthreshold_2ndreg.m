%(c) Julia Riedl 2017
%extracts flourescence intensity values of the second largest object in *.stk files of a rois given by the threshold
%supplied by the param input file.
%param format:
% row 1:worm number; row 2 cherry channel threshold; row 3: gcamp
% background; row 4: cherry background
%example: [1 2;4800 3800;0 0;517 517;517 517]
%worms input vector indicates which recordings within a folder to analyze
%Data should be organized: each recording in a subfolder with the
%searchterm in the folder name( e.g. W1, W2,...etc). *stk files dor Gcamp
%must have 'gc' in file name, for cherry 'ch', or change search string accordingly in
%line 62
%output: thiswormdata2: column 1 :GcaMP values, column 2: mCherry values
%        Trackparameters: parameter input file
%        proofpics: code diplays the thresholded gcam chanel image every
%        100th frame. If indicated in the code these images are stored for
%        posthoc proofreading
%param=xlsread('\\storage.imp.ac.at\groups\zimmer\Julie\CaImaging\Analyze CaImaging Julie.xlsx',-1);
function [imgnum] = freelymovingthreshold_2ndreg(param,worms)


tic
date=datestr(now,29);
%save(['Trackparameters_' date],'-v6')
searchterm='W' %input('folder name part?','s');
folders=dir(['*' searchterm '*']);

home=cd;
remove=0;


counter = 0;
foldercounter = 1; %this counts the number of folders generated, so that even if multiple folders are made for one worm
%(separate frames tracked with different thresholds) then this increases by
%1 for each folder - allows the correct threshold and bckgd values to be
%used!!!!

%parameters:
thresh1=param(2,worms);
bckgds=param(3:4,worms);
%save parameters
tic
date=datestr(now,29);
try
    save(['Trackparameters_' date],'param')
catch
    disp('err')
end

for  w_idx= 1:length(worms)
    
    cd(home)
    wormnum=worms(w_idx)
    stkfolder = [cd '\' folders(wormnum).name];
    cd(stkfolder)
    all_stks=dir('*.stk');
    
    
    for stack = 1:length(all_stks)/2
        
        current_stks= dir(['*.stk']);
        
        % find cherry and gcamp stk:
        stknames={current_stks.name};
        
        ch_i=find( cellfun('isempty', strfind( stknames,'_gc')));
        gc_i=find( cellfun('isempty', strfind( stknames,'_ch')));
        ch_stk=current_stks(ch_i(stack)).name
        gc_stk=current_stks(gc_i(stack)).name
        
        numframes = 7900;% size(gcamp,2); %Get the number of images/frames in the recording. SUBTRACT 3 to get the real number...
        foldername = folders(wormnum).name;
        
        fn=dir([foldername '_' num2str(stack)  '_R1_log*']);
        
        %analyze:
        if isempty(dir([foldername '_' num2str(stack)  '_R2_log*']))
            
            
            currentthreshold=thresh1(w_idx);
            thistrackbckgd = bckgds(:,w_idx)';
            thiswormdata2=NaN;
            [thiswormdata2,imgnum] =getregionvalues_jr1_region2(ch_stk,gc_stk,foldername,...
                currentthreshold,thistrackbckgd,numframes,stack,thiswormdata2);
            
            
        elseif length(dir([foldername '_' num2str(stack)  '_R2_log*']))>0 % if R2 log file already exists, load and start analysis at fist NaN
            
            fn=dir([foldername '_' num2str(stack)  '_R2_log*']);
            load (fn(1).name)
            disp('...log file already exists')
            disp(['...loading  ' fn(1).name])
            currentthreshold = thresh1(w_idx);
            thistrackbckgd = bckgds(:,w_idx)';
            [thiswormdata2,imgnum] =getregionvalues_jr1_region2(ch_stk,gc_stk,foldername,...
                currentthreshold,thistrackbckgd,numframes,stack,thiswormdata2);
            
        end
        
        
        counter = counter + 1;
        
    end % end substack loop
    cd (home)
    toc
    display (['W' num2str(wormnum) '...done'])
    foldercounter = foldercounter + 1;
    
end % end folder loop

%% save parameters

save(['Trackparameters_' date],'param')

end