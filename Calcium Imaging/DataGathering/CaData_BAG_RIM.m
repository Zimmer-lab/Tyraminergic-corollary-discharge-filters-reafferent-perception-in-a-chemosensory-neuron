
%(c) Julia Riedl 2018
%gather all imaging data and corresponding behavior annotation (*log.mat files) into one struct file "ImagingDataAll":

clear
home=cd;
cc=1;

worms=dir('*W*');

for W=1:length(worms)
    
    cd(worms(W).name)
    
    revfiles=dir('*SF*');
    rev_end_files=dir('*rev_end*');
    
    if ~isempty(revfiles)
        rev=1;
    else
        rev=0;
    end
    
    %ratios:
    files=dir('*R1*_log.mat');
    files2=dir('*R2*_log.mat');
    % go through various movies for this worm:
    if ~isempty (files2)
        for movie = 1:length(files)
            
            disp(files(movie).name)
            load (files(movie).name);
            
            disp(files2(movie).name)
            load (files2(movie).name);
            stks=dir('*.stk');
            
            %BAG:
            gcamp_bag=thiswormdata(:,1);
            cherry_bag=thiswormdata(:,2);
            
            %RIM:
            gcamp_rim=thiswormdata2(:,1);
            cherry_rim=thiswormdata2(:,2);
            
            %background:
            if size(thiswormdata2,2)>4
                backgr=thiswormdata2(:,4)./thiswormdata2(:,5);               
            else
                backgr=zeros(length(cherry_rim),1);
            end
                    
            %stimulus:
            if exist('stimulus.mat')==2
                load('stimulus.mat');
            else
                stimulus=NaN;
            end
            
            
            CaData{cc}.cherry_bag=cherry_bag;
            CaData{cc}.gcamp_bag=gcamp_bag;
            CaData{cc}.cherry_rim=cherry_rim;
            CaData{cc}.gcamp_rim=gcamp_rim;
            CaData{cc}.backgr=backgr;
            CaData{cc}.stimulus=stimulus;
            CaData{cc}.TrialLabel=[stks(1).name(1:end-4)];
            CaData{cc}.date=datestr(now);
            cc=cc+1;
        end
    else
        disp('no R2 log file')
    end
    cd ..\
    
    
end % end Worm Loop


%% save:
ImagingData=CaData;
save ImagingDataAll ImagingData







