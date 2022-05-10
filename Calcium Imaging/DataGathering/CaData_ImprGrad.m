%(c) Julia Riedl 2018
%gather all Imaging Data from *.stk data (Visiview) into one struct file
%for coverslip experiments with imprinted gradient:

clear
warning off

plots=1; %proof plots?

ImagingData=cell(1,1);
load CO2fit_impr
XB=3800; %maximum X stage position;
%ml=input('split movie part length?');
ml=7500; % hard coded stk part frame length
mX=input('offset X stage value?');
mY=input('offset Y stage value?');

%maxX=22000; %stagepos maximum value
%imprinted gradient:
gradient=vertcat(CO2fit_impr,CO2fit_impr,CO2fit_impr);
gradient1=gradient;
for i=1:12
    gradient=vertcat(gradient,gradient(1:end/1.05,:));
end

home=cd;
cc=1;

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
            disp( ['split movie part length:' num2str(ml)])
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
            disp('incorrect reversal files')
        end
        
        
        
        if  ~isempty(RevFrames30hz) & movie==2 &rev==1
            ImagingData{cc}.RevFrames30hz=RevAll;
        else
            ImagingData{cc}.RevFrames30hz=revs1;
            
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
            %correct for exceeding x values:
            X=X+mX;
            Y=Y+mY;
            if min(X)<0
                X=X-min(X)+1;
            end
            ImagingData{cc}.XY=[X Y];
        elseif length(posfiles)==1
            disp(posfiles(1).name);
            XYpos=load(posfiles(1).name);
            if isstruct(XYpos)
                XYpos=XYpos.XY;
            end
            X1=1-XYpos(:,1);
            Y1=XYpos(:,2);
            %correct for exceeding x values:
            if min(X1)<0
                X1=X1-min(X1)+1;
            end
            ImagingData{cc}.XY=[X1 Y1];
            
        end
    else
        ImagingData{cc}.XY=NaN(length(ratio),2);
        
    end
    
    %----CO2 values:----
    %normalize to 0 X and round:
    bin=10;
    Xd=round((ImagingData{cc}.XY(:,1))/bin);
    Yd=round((ImagingData{cc}.XY(:,2))/bin);
    if min(Xd)<=0
        Xd=Xd-min(Xd)+1;
    end
    if min(Yd)<=0
        Yd=Yd-min(Yd)+1;
    end
    
    if max(Xd)>size(gradient,2)
        disp('Xd beyond gradient size limits!')
        Xd=Xd-(max(Xd)-size(gradient,2))-bin;
    end
    
    %     %fit gradient to stagepos dimensions:
    gf=gradient(:,1:(size(gradient,2)/((XB))):end);
    
    %Map XY on gradient:
    %     X(X>length(gradient))=NaN;
    C_ind=sub2ind(size(gradient),Yd,Xd);
    
    %get CO2 concentration at this position of path:
    nan_idx=isnan(C_ind);
    C_ind=C_ind(nan_idx~=1);
    
    %get CO2 concentration at this position of path and reinsert Nans:
    gradV=gf(round(C_ind));
    gradVal=NaN(1,length(nan_idx));
    gradVal(nan_idx~=1)=gradV;
    
    ImagingData{cc}.cherry_bag_raw=cherry_o;
    ImagingData{cc}.gcamp_bag_raw=gcamp_o;
    ImagingData{cc}.cherry_bag=cherry;
    ImagingData{cc}.gcamp_bag=gcamp;
    ImagingData{cc}.ratio=ratio;
    ImagingData{cc}.CO2=gradVal;
    ImagingData{cc}.TrialLabel=files(movie).name(1:end-7);
    cc=cc+1;
    
    if plots==1
        if W==1
            figure
        end
        scatter(Xd,Yd,10,gradVal)
        hold on
        scatter(Xd(1),Yd(1),'k')
        text(Xd(1),Yd(1),num2str(W))
        colorbar
    end
    
    
    cd ..\
    
end % end Worm Loop


%% save:
disp('saving data')
save ImagingDataAll ImagingData






