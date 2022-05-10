%(c) Julia Riedl 2018
%gather all imaging data and corresponding behavior annotation (*log.mat files)
 %into one struct file "ImagingDataAll":

clear
files = dir('*BAG*');
CaData=cell(1,1);

home=cd;
cc=1;
ml=NaN;

worms=dir('*W*');


for W=1:length(worms)
    
    if ~isdir(worms(W).name)
        continue
    end
    
    cd(worms(W).name)
    
    revfiles=dir('*SF*');
    rev_end_files=dir('*rev_end*');
    
    if ~isempty(revfiles)% & ~isempty(rev_end_files)
        rev=1;
    else
        rev=0;
    end
    
    %ratio:
    files=dir('*R1*_log.mat');
    if isempty(files)
        disp([worms(W).name ': no log file'])
        cd ..\
        continue
    end
    % go through various movies for this worm and concatenate if needed:
    
    for movie = 1:length(files)
        
        disp(files(movie).name)
        load (files(movie).name);
        
        gcamp=thiswormdata(:,1);
        cherry=thiswormdata(:,2);
        
        if movie==1
            gcamp_o=thiswormdata(:,1);
            cherry_o=thiswormdata(:,2);
        else
            gcamp_o=[gcamp_o;gcamp];
            cherry_o=[cherry_o;cherry];
        end
        
        ratio=(thiswormdata(:,1)./thiswormdata(:,2));
        %remove out of focus episodes: find very rapid and drastic changes in size of cherry area
        %and remove this data from both cherry and gcamp:
        areaR1=thiswormdata(:,3);
        areaR1=nandetrend(areaR1);
        areaJumps1=find(diff(medfilt1(areaR1))>300 | diff(medfilt1(areaR1))<-300);
        areaJumps1=vertcat(areaJumps1, find(areaR1>2*nanmedian(areaR1) | areaR1<0.2*nanmedian(areaR1)));
        if numel(areaJumps1)>100
%             figure
%             hold on
%             plot(areaR1)
%             scatter(areaJumps1,(areaR1(areaJumps1)))
%             ylabel('area')
%             title(W)
%             pause(0.5)
        end
        gcamp(areaJumps1)=NaN;
        cherry(areaJumps1)=NaN;
        
        
        if length(files)>1 & movie==1
            if isnan(ml)
                ml=input('split movie part length?');
            end
            if length(gcamp)<ml
                ml=length(gcamp)
            end
            gcamp1=gcamp(1:ml);
            cherry1=cherry(1:ml);
            ratio1=ratio(1:ml);
            
        end
        
        if length(files)>1 & movie==2
            
            gcamp=[gcamp1; gcamp];
            cherry= [cherry1; cherry];
            ratio=[ratio1; ratio];
        end
        
        %annotated reversals:
        RevFrames30hz=NaN;
        if rev==1  & movie==1 %get reversals for first part or if only 1 reversal file
            load (revfiles(movie).name)
            disp (revfiles(movie).name)
            
            revs1=RevFrames30hz(:);
            
        elseif rev==1 & movie==2 & length(revfiles)==2
            
            load (revfiles(movie).name)
            disp (revfiles(movie).name)
            %concatenate with first:
            RevAll=[revs1 ; RevFrames30hz(:)+ml];
            
        elseif rev==1 & movie==2 & length(revfiles)==1
            RevAll=revs1;
        else
            revs1=NaN;
            beep
            disp('incorrect reversal files')
        end
        
 
            
            if  ~isempty(RevFrames30hz) & movie==2 &rev==1
                CaData{cc}.RevFrames30hz=RevAll;
            else
                CaData{cc}.RevFrames30hz=revs1;

            end
        
    end
    
    %XY:
    posfiles=dir('*stagepos*');
    
    if ~isempty(posfiles)
        if length(posfiles)==2
            disp(posfiles(1).name);
            XYpos=load(posfiles(1).name);
            if isstruct(XYpos)
                XYpos=XYpos.XY;
            end
            X1=1-XYpos(:,1);
            Y1=XYpos(:,2);
            disp(posfiles(2).name);
            XYpos=load(posfiles(2).name);
            if isstruct(XYpos)
                XYpos=XYpos.XY;
            end
            X=1-XYpos(:,1);
            Y=XYpos(:,2);
            X=[X1;X];
            Y=[Y1;Y];
            CaData{cc}.XY=[X Y];
        elseif length(posfiles)==1
            disp(posfiles(1).name);
            XYpos=load(posfiles(1).name);
            if isstruct(XYpos)
                XYpos=XYpos.XY;
            end
            X1=1-XYpos(:,1);
            Y1=XYpos(:,2);
            CaData{cc}.XY=[X1 Y1];
        end
    end
    
    CaData{cc}.cherry_bag_raw=cherry_o;
    CaData{cc}.gcamp_bag_raw=gcamp_o;
    CaData{cc}.cherry_bag=cherry;
    CaData{cc}.gcamp_bag=gcamp;
    CaData{cc}.ratio=ratio;
    CaData{cc}.TrialLabel=files(movie).name(1:end-7);
    cc=cc+1;
    
    cd ..\
    
end % end Worm Loop




%% save:
ImagingData=CaData;
save ImagingDataBAG ImagingData
save ImagingDataAll ImagingData






