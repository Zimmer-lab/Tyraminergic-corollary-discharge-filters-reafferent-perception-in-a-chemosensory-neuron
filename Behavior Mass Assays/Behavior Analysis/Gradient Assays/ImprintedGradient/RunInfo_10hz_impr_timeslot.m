%(C) Julia Riedl
% for 10hz high resolution behavioral recordings
% extracts all forward run episodes of all Tracks, downsamples to 3 hz and stores the
%corresponding parameters (X,Y,speed,etc) in a "runinfo" file.

clearvars -except fo folders
hz=10;
max_time=30*60*hz;
min_time=1*60*hz;
center=900;
bin=3; %binning factor
minrl=3; %minimum reversal length in sec

disp(['max time: ' num2str(max_time) 'frames'])
boundary=890;
files=dir('*NI_als.mat'); %search string
if isempty(files)
    files=dir('*als.mat');
end
controlplots=0;
CO2=1800; % x of CO2 max


% cd ..\
for F=1:length(files)
    
    disp(F)
    clearvars -except CO2 bin hz minrl F center boundary exp  home name files bearing dBearing DistanceTrav speed_ctr RunFrames...
        max_time min_time fo folders controlplots
    cc=0;
    
    fname=(files(F,1).name);
    disp(fname)
    load (files(F).name);
    %load([fname(1:end-4), '_anglesV8.mat']);
    
    
    pixToMM=0.0155;
    
    Ctr_dBearing=cell(1);
    Ctr_Bearing=cell(1);
    Ctr_Speed=cell(1);
    Ctr_Dist=cell(1);
    dSensCtr={NaN};
    dSensHead={NaN};
    HeadRunsX={NaN};
    HeadRunsY={NaN};
    SensoryRunsHead={NaN};
    head_speed={NaN};
    head_speed_net={NaN};
    
    %% ---analyze Tracks----
    
    for T= 1:length(Tracks)
        % only tracks which start after establishment of the gradient
        % (ca frame 500)
        if  Tracks(1,T).Frames(end)>min_time & Tracks(1,T).Frames(1)<max_time
            
            % --plot tracks--
            
            %             hold on
            %             plot( Tracks(1,T).SmoothX, Tracks(1,T).SmoothY,'Color',[0.5 0.5 0.5]);
            %             title(T);
            
            
            
            %% get indices of all behavior events for excluding them from the
            % trajectories:
            
            Rev1=[];
            Omegas=[];
            if isfield(Tracks(1,T),'polishedReversals')
                if size(Tracks(1,T).polishedReversals,1)>1
                    gi=find(Tracks(1,T).polishedReversals(:,4)==1);
                    for ii =1:length(gi)
                        if Tracks(1,T).polishedReversals(gi(ii),2)-Tracks(1,T).polishedReversals(gi(ii),1)>minrl*hz
                            Rev=(Tracks(1,T).polishedReversals(ii,1):Tracks(1,T).polishedReversals(ii,2));
                            Rev1=horzcat(Rev1,Rev);
                        end
                    end
                end
            else
                if size(Tracks(1,T).Reversals,1)>1
                    gi=find(Tracks(1,T).Reversals(:,4)==1);
                    
                    for ii =1:length(gi)
                        if Tracks(1,T).Reversals(gi(ii),2)-Tracks(1,T).Reversals(gi(ii),1)>minrl*hz
                            Rev=(Tracks(1,T).Reversals(ii,1):Tracks(1,T).Reversals(ii,2));
                            Rev1=horzcat(Rev1,Rev);
                        end
                    end
                else Rev=[];
                end
            end

        
        RunX=  (Tracks(1,T).SmoothX);
        RunY=  (Tracks(1,T).SmoothY);
        AngVelocity=abs(Tracks(1,T).AngSpeed);
        RunX(Rev1)=NaN;
        RunY(Rev1)=NaN;
        
        
        %----exclude those which occur close to border:----
        d=[];
        for i=1:length(RunY)
        d(i)=abs(sqrt(((center-RunX(i)).^2)+((center-RunY(i)).^2)));
    end
    bi=find(d>boundary);
    RunX(bi)=NaN;
    RunY(bi)=NaN;
    
    % find start and end points, avoid first and last 4 points at end of run,
    % they often belong already to reversal or omega:
    ci = isnan(RunX);
    RunEnd=find(diff(ci)>0)-5;
    RunStart=find(diff(ci)<0)+5;
    if length(ci)>0
        if ci(1)==0;
            RunStart=cat(2,1,RunStart);
        end
        if ci(end)==0;
            RunEnd=cat(2,RunEnd,length(RunX));
        end
    end
    
    
    % find runs longer than X
    X=3;
    longRuns=find((RunEnd-RunStart)>X);
    RunStart=RunStart(longRuns);
    RunEnd=RunEnd(longRuns);
    
    %kill runs which have unusual high angular speed
    %and don't cover 1 small worm travel distance (20 px):
    meanAngVelocity1=NaN(1,1);
    for ii= 1:length(RunStart)
        meanAngVelocity1(ii)=(nanmean(abs(Tracks(1,T).AngSpeed(RunStart(ii):RunEnd(ii)-1))));
        if meanAngVelocity1(ii)>30
            RunStart(ii)=NaN;
            RunEnd(ii)=NaN;
        else
            X=RunX(RunStart(ii):RunEnd(ii));
            Y=RunY(RunStart(ii):RunEnd(ii));
            d=NaN(1,1);
            ci=1;
            for iii=1:bin:length(Y)
                d(ci)=sqrt(((X(1)-X(iii)).^2)+((Y(1)-Y(iii)).^2));
                ci=ci+1;
            end
            if max(d)<20
                RunStart(ii)=NaN;
                RunEnd(ii)=NaN;
            end
        end
    end
    
    RunStart=RunStart(~isnan(RunStart));
    RunEnd=RunEnd(~isnan(RunEnd));
    RunX1={NaN};
    RunY1={NaN};
    for ii=1:length(RunStart)
        RunX1{ii}=RunX(RunStart(ii):bin:RunEnd(ii)); %downsample
        RunY1{ii}=RunY(RunStart(ii):bin:RunEnd(ii));
        RunFrames1{T,ii}=Tracks(1,T).Frames((RunStart(ii):bin:RunEnd(ii)));
    end
    
    %--plot--
    %             if length(RunStart)>1
    %                 for ii= 1:length(RunStart)
    %                     hold on
    %                     %plot([RunX(RunStart(ii)) RunX(RunStart(ii)+6)],[RunY(RunStart(ii)) RunY(RunStart(ii)+6)],'LineWidth',3);
    %                     plot(RunX(RunStart(ii):RunEnd(ii)),RunY(RunStart(ii):RunEnd(ii)),'k');
    %                     %scatter(RunX(RunStart(ii)),RunY(RunStart(ii)) ,3.5, HeadingStart1(ii),'filled');
    %
    %                     %colorbar
    %                 end
    %             end
    
    
    
    %--- get all bearing and heading angles for each left over run relative to
    %%  preferred 02 isocline:------
    if ~isnan(RunX1{1})
        
        for R=1:length(RunX1)
            
            x=RunX1{R}(1:end);
            y=RunY1{R}(1:end);
            s=Tracks(1,T).Speed(RunStart(R):RunEnd(R));
            s=s(1:bin:end);
            L=length(x);
            bearingangles=NaN(1,1);
            headingangles=NaN(1,1);
            beelineX=NaN;
            D=NaN;
 
            
            for ii=1:L-1
                beelineX=(CO2-(x(ii)));
                %                     if F==5
                %                          beelineX(ii)=(roi(1,1)-(y(ii)));
                %                     end
                beelineY=0;%(roi(2,2)+100-(x(ii)));  % alternatively bearing towards gas inlet
                headingX=(x(ii+1))-(x(ii));
                headingY=(y(ii+1))-(y(ii));
                nenner=((beelineX*headingX)+(beelineY*headingY));
                Brecher=(sqrt((beelineX).^2+(beelineY).^2))*(sqrt((headingX).^2+(headingY).^2));
                cosinus_angle=nenner/Brecher;
                bearingangles(ii)=(acos(cosinus_angle))*(360/(2*pi));
                headingangles(ii) = atan2(headingY,headingX);
                D(ii)=sqrt(((x(ii)-x(ii+1)).^2)+((y(ii)-y(ii+1)).^2));
            end
            
            cc=cc+1;
            
            Ctr_Bearing{T,R}=[bearingangles NaN ];
            db=[diff(bearingangles) NaN NaN ];
            db(db>10 | db<-10)=NaN;
            Ctr_dBearing{T,R}=[db];
            Ctr_Speed{T,R}=s;
            Ctr_Dist{T,R}=[D NaN];
            
            
        end % end Run loop
        
    end
    
    %% --plot tracks--
    if controlplots==1 & T>66 & F==7 && T<80 & length(Tracks(1,T).SmoothY)>10
        %
        hold on
        scatter(Tracks(1,T).SmoothX(1),Tracks(1,T).SmoothY(1),'dk')
        
        for R=1:length(RunX1)
            try
                scatter (RunX1{R}(1:end),RunY1{R}(1:end),30,(Ctr_dBearing{T,R}./(Ctr_Speed{T,R})),'filled')
                %                         scatter(HeadRunsX{T,R},HeadRunsY{T,R},13,head_speed{T,R},'filled')
            catch;
            end
        end
        plot( Tracks(1,T).SmoothX, Tracks(1,T).SmoothY,'Color',[0.5 0.5 0.5]);
        title([T,F])
        colorbar
        caxis([-10 10])
        
    end
    
end

end %end Tracks loop
%%

bearing{F}=Ctr_Bearing;
dBearing{F}=Ctr_dBearing;
speed_ctr{F}=Ctr_Speed;
RunFrames{F}=RunFrames1;
DistanceTrav{F}=Ctr_Dist;

end
%%

display('...save')
beep
als_date=datestr(now);

save(['runinfo_min' num2str(minrl) 'srevs_' dirname(cd) '_to' num2str(max_time/600) 'min'] , ...
    'bearing' ,'dBearing' ,'speed_ctr','RunFrames','als_date','DistanceTrav');



