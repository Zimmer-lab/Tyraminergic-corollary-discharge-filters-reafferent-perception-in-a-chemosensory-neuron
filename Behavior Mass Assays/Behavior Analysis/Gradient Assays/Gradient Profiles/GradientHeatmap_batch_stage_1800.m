%(c) Julia Riedl 2016
% bins worm's occupancy and reversal frequency, speed and angular speed over area of gradient arena
%and saves as 3D matrix for given amount of timepoints 

clearvars -except Tracks

% ---for any avi movies:---
als_folder=dir('*analysis*');
files= dir('*JR_als.mat');

%%%%things to edit:%%%
dim=[1800 1800];%input('Video dimensions XY?');
frame_number=12000;
hz=10; % video frame rate--> for 1 s bins
timebin2=10; % total time bin in sec


for Movie=1:2% length(files)
    Movie
    display( strcat('...current file...',files(Movie).name))
    
    %find and load als file
    
    files= dir('*JR_als.mat');
%     load (files(Movie).name);
    
    %% ---(1)get all XY positions and behavior events for X frames
    cc=1;
    display ('...gather paths of worms')

    for frame=1:hz:frame_number-hz %  nr of frames in the video
        XYpos1=NaN(1,2);
        XYspeed1=NaN(1,1);
        XYangSpeed1=NaN(1,1);
        c=1; %counter
        for i=1:length(Tracks)
            
            if find (Tracks(i).Frames== frame)  & ~isempty(Tracks(i).Speed)
               
                fi=find (Tracks(i).Frames== frame);

                if fi<= length(abs(Tracks(i).AngSpeed))-hz;
                    XYpos1(c,1:2)=round(nanmean(Tracks(i).Path(fi:fi+hz,:))*10)/10;
                    XYspeed1(c,1)=nanmean(Tracks(i).Speed(fi:fi+hz));
                    XYangSpeed1(c,1)=nanmean(abs(Tracks(i).AngSpeed(fi:fi+hz)));
                    c=c+1;
                end
            end
        end
        % store binned values in cell:
        XYpos{cc}=XYpos1;
        XYspeed{cc}=XYspeed1;
        XYangSpeed{cc}=XYangSpeed1;
        cc=cc+1;
        
        if mod(frame,1+hz*10)==0
            disp(['frame..' num2str(frame)])
        end
    end
    
    %% ---put them into XY bins----
    display ('...bin')
    ap=max(cell2mat(XYpos')); %maximum x dimension
    binsize=floor(round(ap(2))/30);
    c=0;
    for kk=1:timebin2:length(XYpos)-timebin2
        
        c=c+1;
        XYbin=[];
        XYspeedbin=[];
        XYangSpeedbin=[];
        XYbinC=[];
        for k=kk:kk+timebin2-1
            clear XYbin1
            XY=(XYpos{1,k});
            XYs=(XYspeed{1,k});
            XYas=(XYangSpeed{1,k});
            XY=XY/binsize;
            XYbin1=ceil(XY);
            XYbin=vertcat(XYbin,XYbin1);
            XYspeedbin=vertcat(XYspeedbin,XYs);
            XYangSpeedbin=vertcat(XYangSpeedbin,XYas);
        end
        
        
        M=zeros(floor(dim/binsize));
        m=zeros(1,length(XYbinC));
        for a=1:ceil(dim(1)/binsize)
            for b=1:ceil(dim(2)/binsize)
                
                m=find(XYbin(:,1)==a & XYbin(:,2)==b);
                M(b,a)=length(m);
                SpeedM(b,a)=sum(XYspeedbin(m))/length(m);
                try
                    angSpeedM(b,a)=sum(XYangSpeedbin(m))/length(m);
                end
            end
        end
        M_norm=M/sum(sum(M)/100);
        HM_all(:,:,c)=M;
        HMnorm_all(:,:,c)=M_norm;
        SpeedM_all(:,:,c)=SpeedM;
        angSpeedM_all(:,:,c)=angSpeedM;
        
    end
    HM_batch{Movie}=HM_all;
    HMnorm_batch{Movie}=HMnorm_all;
    SpeedM_batch{Movie}=SpeedM_all;
    angSpeedM_batch{Movie}=angSpeedM_all;
    
end

%% put all data into analysis folder:
d= strfind(cd, '\');
nd=(cd);
name=nd(d(end)+1:end);
mkdir(strcat(name,'_profile_analysis'));
cd(strcat(name,'_profile_analysis'));
display('...save')
save HM_batch HM_batch
save HMnorm_batch HMnorm_batch
save SpeedM_batch SpeedM_batch
save angSpeedM_batch angSpeedM_batch






