%(c) Julia Riedl 2017

%param=xlsread('\\storage.imp.ac.at\groups\zimmer\Julie\CaImaging\Analyze CaImaging Julie.xlsx',-1);
function [imgnum] = freelymovingthreshold_jr1_2ndreg(param,worms)

tic
%save(['Trackparameters_' date],'-v6')
searchterm='W' %input('folder name part?','s');
folders=dir(['*' searchterm '*']);

home=cd;
counter = 0;
foldercounter = 1; %this counts the number of folders generated

%parameters:
thresh1=param(2,worms);
thresh2=param(3,worms);
bckgds=param(4:5,worms);
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
    separateframes = 0;
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
        
        % (1)if there is yet no log file
               
        if isempty(dir([foldername '_' num2str(stack)  '_R2_log*']))
            
            
            thistrackthreshold = thresh1(w_idx);
            currentthreshold2=thresh2(w_idx);
            thistrackbckgd = bckgds(:,w_idx)';
            thiswormdata2=NaN;
            [thiswormdata2,imgnum] =getregionvalues_jr1_region2(ch_stk,gc_stk,foldername,...
                currentthreshold2,numframes,stack,thiswormdata2);
            
%---(2) if there is a log file from prevoius analysis-----
        elseif length(dir([foldername '_' num2str(stack)  '_R2_log*']))>0 
            
            fn=dir([foldername '_' num2str(stack)  '_R2_log*']);            
            load (fn(1).name)
            disp('...log file already exists')
            disp(['...loading  ' fn(1).name])
            thistrackthreshold = thresh1(w_idx);
            currentthreshold2=thresh2(w_idx);
            thistrackbckgd = bckgds(:,w_idx)';
            [thiswormdata2,imgnum] =getregionvalues_jr1_region2(ch_stk,gc_stk,foldername,...
    currentthreshold2,numframes,stack,thiswormdata2)
            
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