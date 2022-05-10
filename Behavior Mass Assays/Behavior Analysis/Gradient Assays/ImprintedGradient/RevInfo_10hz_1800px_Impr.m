%(c) Julia Riedl 2019-21
%gathers data for reversals from Tracks.mat files to a single struct to be
%used by other analysis codes
%type of data: for each reversal:
% -Frames Nr
% -Seconds
% -bearing
% -diff(bearing)
% -speed
% -reversal onset indices
% -omega onset indices

%----
clear
warning off
%----parameters:---
bin=10; %binning factor
maxtime=50; %in minutes
hz=10; %recording rate
boundary=850;
center=1000;
reversals={NaN};

files=dir('*NI_als.mat'); %search string
if isempty(files)
    files=dir('*als.mat');
end

for F=1:length(files)
    
    clearvars -except F exp maxtime  files bearing dBearing speed_ctr...
        HeadSens HeadSpeed HeadSpeedNet reversals bin gradient sensPath boundary omegas seconds  center
    cc=0;
    disp(F)
    load (files(F).name);
    disp((files(F).name));
    Reversals={NaN};
    dCdT={NaN};
    Second={NaN};
    
    %% ---analyze Tracks----
    
    for T= 1:length(Tracks)

        if  Tracks(1,T).Frames(end)>100 & Tracks(1,T).Frames(1)<maxtime*10*60
            
            TX=  (Tracks(1,T).SmoothX);
            TY=  (Tracks(1,T).SmoothY);
            
            %----exclude those which occur close to border:----
            %remove the ones outside  boundaries;
            d=[];
            for i=1:length(TX)
                d(i)=abs(sqrt(((center-TX(i)).^2)+((center-TY(i)).^2)));
            end
            bi=find(d>boundary);
            TX(bi)=NaN;
            TY(bi)=NaN;
            
            %kill runs which have unusual high angular speed
            %and don't cover 1 small worm travel distance (20 px):
            
            meanAngVelocity=(nanmean(abs(Tracks(1,T).AngSpeed)));
            if meanAngVelocity>30
                continue
            else
                d=NaN(1,1);
                ci=1;
                for iii=1:3:length(TY)
                    d(ci)=sqrt(((TX(1)-TX(iii)).^2)+((TY(1)-TY(iii)).^2));
                    ci=ci+1;
                end
                if max(d)<20
                    disp('no displacement')
                    continue
                end
            end
            
            
            % reversal vector:
            Revs=zeros(1,ceil(length(TX)/bin));
            try
                Rev=(Tracks(1,T).polishedReversals);
            catch
                Rev=(Tracks(1,T).Reversals);
            end
            
            
            if ~isempty(Rev) & max(d)>=20;
                gi=find(Rev(:,4)==1);
                Rev=Rev(gi,:);
                Rev=Rev(:,1);
                Rev=round(Rev/bin);
                Revs(Rev)=1;
                
            end
            %remove the ones outside  boundaries;
            bi(bi<3)=[];
            Revs(round(bi/bin))=0;
            Reversals{T}=Revs;
            
            % omegas vector:
            omv=zeros(1,ceil(length(TX)/bin));
            if ~isempty(Tracks(1,T).OmegaTrans)
                os=Tracks(1,T).OmegaTrans(:,1);
                os=round(os/bin);
                omv(os)=1;
            end
            %remove the ones outside  boundaries;
            omv(round(bi/bin))=0;
            Omega{T}=omv;
            
            %seconds:
            Frame{T}=round(Tracks(T).Frames(1:bin:end));
            sec=round(Tracks(T).Frames(1:bin:end)/hz);
            sec(round(bi/bin))=NaN;
            %remove the ones outside  boundaries;
            Second{T}=sec;
            
            %--- get all bearing and heading angles for each left over run relative to
            %%  preferred C02 isocline:------
            if   max(d)>=20
                
                x=TX(1:bin:end);
                y=TY(1:bin:end);
                s=Tracks(1,T).Speed;
                s=s(1:bin:end);
                L=length(x);
                bearingangles=NaN(1,1);
                headingangles=NaN(1,1);
                beelineX=NaN;
                
                CO2=1;
                
                for ii=1:L-1
                    beelineX=(CO2-(x(ii)));
                    beelineY=0;%(roi(2,2)+100-(x(ii)));  % alternatively bearing towards gas inlet
                    headingX=(x(ii+1))-(x(ii));
                    headingY=(y(ii+1))-(y(ii));
                    nenner=((beelineX*headingX)+(beelineY*headingY));
                    Brecher=(sqrt((beelineX).^2+(beelineY).^2))*(sqrt((headingX).^2+(headingY).^2));
                    cosinus_angle=nenner/Brecher;
                    bearingangles(ii)=(acos(cosinus_angle))*(360/(2*pi));
                    headingangles(ii) = atan2(headingY,headingX);
                end
                
                cc=cc+1;
                
                Ctr_Bearing{T}=[bearingangles NaN ];
                Ctr_dBearing{T}=[diff(bearingangles) NaN NaN ];
                Ctr_Speed{T}=s;
                
            else
                
                Ctr_Bearing{T}=NaN(1,length(omv));
                Ctr_dBearing{T}=NaN(1,length(omv));
                Ctr_Speed{T}=NaN(1,length(omv));
                
                
            end
            
        end
        
    end %end Tracks loop
    %%
    
    bearing{F}=Ctr_Bearing;
    dBearing{F}=Ctr_dBearing;
    speed_ctr{F}=Ctr_Speed;
    reversals{F}=Reversals;
    omegas{F}=Omega;
    seconds{F}=Second;
    frames{F}=Frame;
    Filenames{F}=files(F).name;
    date_analyzed=datestr(now);
    
end
%%

home=cd;
nd=(cd);
d= strfind(cd, '\');
name=nd(d(end)+1:end);

display('...save')

save(['rev_info_10hz_JR_'   name] , 'date_analyzed','Filenames' ,'bearing' ,'dBearing' ,'speed_ctr','reversals','bin','omegas','seconds','frames');



