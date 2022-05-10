
%---(c) Julia Riedl 2017----
%extracts flourescence intensity values of *.stk files of a roi given by the threshold
%supplied by the param input file.
%param format:
% row 1:worm number; row 2 cherry channel threshold; row 3: optional 2nd threshold; row 4: gcamp
% background; row 5: cherry background
%example: [1 2;4800 3800;0 0;517 517;517 517]
%worms input vector indicates which recordings within a folder to analyze
%Data should be organized: each recording in a subfolder with the
%searchterm in the folder name( e.g. W1, W2,...etc). *stk files dor Gcamp
%must have 'gc' in file name, for cherry 'ch', or change search string accordingly in
%line 62
%output: thiswormdata: column 1 :GcaMP values, column 2: mCherry values
%        Trackparameters: parameter input file
%        proofpics: code diplays the thresholded gcam chanel image every
%        100th frame. If indicated in the code these images are stored for
%        posthoc proofreading


function [imgnum] = freelymovingthreshold_offset_crct(param,worms)

date=datestr(now,29);
save(['Trackparameters_' date],'param','-v6')
searchterm='W' ;%input('folder name part?','s');
folders=dir([searchterm '*']);

home=cd;

%parameters:
thresh1=param(2,worms);
thresh2=param(3,worms);
bckgds=param(4:5,worms);

for  w_idx= 1:length(worms)
    cd(home)
    wormnum=worms(w_idx);
    stkfolder = [cd '\' folders(wormnum).name];
    cd(stkfolder)
    all_stks=dir('*.stk');
    

    for stack = 1:length(all_stks)/2
        
        current_stks= dir(['*.stk']);
        
        % find cherry and gcamp stk:
        stknames={current_stks.name};
        for i=1:length(stknames)
            stknames{i}=stknames{i}(end-10:end);
        end
        
        ch_i=find( cellfun('isempty', strfind( stknames,'_gc')));
        gc_i=find( cellfun('isempty', strfind( stknames,'_ch')));
        if isempty(ch_i)
        ch_i=find( cellfun('isempty', strfind( stknames,'_gc1')));
        gc_i=find( cellfun('isempty', strfind( stknames,'_ch1')));
        end
        
            
        ch_stk=current_stks(ch_i(stack)).name
        gc_stk=current_stks(gc_i(stack)).name
        
        numframes = 15000;% size(gcamp,2); %Get the number of images/frames in the recording. 
        
        foldername = folders(wormnum).name;
        trackinstruction = 0;
        
        fn=dir([foldername '_' num2str(stack)  '_R1_log*']);
        
        leadzeros = '5'; %assuming if it's over 10,000 it's still below 100,000
       
        %analyze:
        if  length(fn)<1 & length(dir([foldername '_' num2str(stack)  '_R1_log*']))<1
                        
            thistrackthreshold = thresh1(w_idx);
            currentthreshold2=thresh2(w_idx);
            thistrackbckgd = bckgds(:,w_idx);
            thiswormdata=NaN;
            
            [thiswormdata,thiswormdata2,imgnum] =getregionvalues_ofs(ch_stk,gc_stk,foldername,trackinstruction,...
                thistrackthreshold,currentthreshold2,thistrackbckgd,numframes,stack,thiswormdata);            
            
        elseif length(dir([foldername '_' num2str(stack)  '_R1_log*']))>0
            %if log file already exists, load it and start only from first
            %NaN
            load (fn(1).name)
            disp('...log file already exists')
            disp(['...loading  ' fn(1).name])
            thistrackthreshold = thresh1(w_idx);
            currentthreshold2=thresh2(w_idx);
            thistrackbckgd = bckgds(:,w_idx);
            [thiswormdata,thiswormdata2,imgnum] =getregionvalues_ofs(ch_stk,gc_stk,foldername,wormnum,leadzeros,trackinstruction,...
                thistrackthreshold,currentthreshold2,thistrackbckgd,numframes,stack,thiswormdata);
            
        end        
    
end % end substack loop
cd (home)
toc
display (['W' num2str(wormnum) '...done'])

end % end folder loop


end