%(c) Julia Riedl 2018
% calculates % of tracks in reversal and first reversal latency of all tracks over time for various tracks files
%input: Tracks*als.mat files from Behavior Tracker CO2 stimulus
%characteristics
clear

files= dir('*als.mat');

%---to edit:---
hz=10;
CO2on=(2440:3000:16500); %CO2 stimulus times
up=60;% CO2 stimulus on time in sec

meanspeed_bin=NaN;
norm_revs=NaN;
rev_bin=NaN;

r_latency=NaN(length(files)*5,30);
r_latency_off=NaN(length(files)*5,30);
miri=NaN;

bc=1;
sc=1;
sco=1;

for F=1:length(files)
    tic
    disp(['loading..' ,files(F).name]);
    load(files(F).name);
    toc
    IRIall=NaN(300,10);
    IRIall_CO2off=NaN(300,10);
    ic=1;
    oc=1;
    %preallocate
    cc=1;
    speed_matrix=NaN(300,30*60*hz);
    rev_matrix=NaN(300,30*60*hz);
    animals=NaN(300,30*60*hz);
    tc=NaN;
    
    for T=1 :length(Tracks)
        
        TX=  (Tracks(1,T).SmoothX);
        TY=  (Tracks(1,T).SmoothY);
        %----exclude those which occur close to border:----
        TX(TX<20)=NaN;
        TY(TY<20)=NaN;
        TX(TX>1780)=NaN;
        TY(TY>1780)=NaN;
        
        %kill runs which have unusual high angular speed
        %and don't cover 1 small worm travel distance (20 px):
        meanAngVelocity=(nanmean(abs(Tracks(1,T).AngSpeed)));
        if meanAngVelocity>30 | Tracks(T).Frames(end)>Tracks(end).Frames(end);
            continue
        else
            d=NaN(1,1);
            ci=1;
            for iii=1:10:length(TY)
                d(ci)=sqrt(((TX(1)-TX(iii)).^2)+((TY(1)-TY(iii)).^2));
                ci=ci+1;
            end
            if max(d)<20
                disp('no displacement')
                continue
            end
        end
        %speed:
        if ~isempty(Tracks(T).Speed)
            speed_matrix(cc,1:30*60*hz)=NaN;
            speed_matrix(cc,Tracks(T).Frames(1):Tracks(T).Frames(end))=Tracks(T).Speed;
        end
        %reversals:
        RevV=zeros(1,length(Tracks(T).Frames));
        Rev=(Tracks(1,T).polishedReversals);
        frames=Tracks(T).Frames;
        
        if ~isempty(Rev) & max(d)>=20;% only 1 sec long reversals
            gi=find(Rev(:,4)>0);
            Rev=Rev(gi,:);
            for i=1:size(Rev,1)
                RevV(Rev(i):Rev(i)+(Rev(i,2)-Rev(i,1)))=1;
            end
        end
        
        rev_matrix(cc,frames(1):frames(end))=RevV;
        animals(cc,frames)=1;
        
        %remove reversals from speed
        ri=find(rev_matrix(cc,:)==1);
        speed_matrix(cc,ri)=NaN;
        speed_matrix(cc,ri-1)=NaN;
        tc(cc)=T;
        cc=cc+1;
        
    end
    %%
    num_animals=nansum(animals);
    meanspeed=medfilt1(nanmeanJ(speed_matrix),3);
    bi=find(speed_matrix==0);
    speed_matrix(bi)=NaN;
    rev_matrix(bi)=NaN;
    all_revs=nansum(rev_matrix,1);
    std_speed=medfilt1(nanstd(speed_matrix),3);
    sv=(1:length(speed_matrix))/hz;

    %  bin to 1s:    
    bin=10;
    c=1;
    for i=1:bin:length(all_revs)-bin
        meanspeed_bin(F,c)=nanmeanJ(meanspeed(i:i+bin));
        rev_bin(F,c)=ceil(nanmeanJ(all_revs(i:i+bin)))/nanmeanJ(num_animals(i:i+bin));
        c=c+1;
    end
    
    
    % shift-triggererd_averages:
    for jj=CO2on/bin-5;
        st_speed(bc,:)=meanspeed_bin(F,jj:jj+up*2);
        st_rev(bc,:)= rev_bin(F,jj:jj+up*2);
        bc=bc+1;
    end
    
    %% inter reversal intervals when CO2 up :

    for jj=CO2on
            lc=1;
        for i=1:cc-1
            if Tracks(tc(i)).Frames(1)<jj & Tracks(tc(i)).Frames(end)>jj+400
                dr=diff(rev_matrix(i,jj:jj+(up*10)));
                on=find(dr==1);
                off=find(dr==-1);
                % latency:
                if ~isempty(on)
                    r_latency(sc,lc)=on(1);
                    lc=lc+1;
                else
                    r_latency(sc,lc)=400;
                    lc=lc+1;
                end
                if  ~isempty(on) & ~isempty(off) & length(on)>1
                    if  on(1)<off(1)
                        fi=[];
                        if length(on)>length(off)                                                   
                        IRI=abs(off(1:end)-on(2:end));                       
                        else
                            IRI=abs(off(1:end-1)-on(2:end)); 
                            if off(end)<530  & ~isnan(dr(end))
                                IRI=[IRI up*10-off(end)];
                            end
                        end
                            
                    else
                        IRI=abs(on-off(1:length(on)));
                    end
                    IRIall(ic,1:length(IRI))=IRI;
                    ic=ic+1;
                elseif length(on)==1 & length(off)==1 & on(1)<off(1) & ~isnan(dr(end))
                    IRIall(ic,1)=(up*10)-off(1);
                    ic=ic+1;
                    
                end
            end
        end
        sc=sc+1;
    end
    %---CO2 down----
    for jj=CO2on+600
            lco=1;
        for i=1:cc-1
            if Tracks(tc(i)).Frames(1)<jj & Tracks(tc(i)).Frames(end)>jj+600
                dr=diff(rev_matrix(i,jj:jj+(up*30)));
                on=find(dr==1);
                off=find(dr==-1);
                % latency:
                if ~isempty(on)
                    r_latency_off(sco,lco)=on(1);
                    lco=lco+1;
                else
                     if Tracks(tc(i)).Frames(end)<jj+900
                    r_latency_off(sco,lco)=Tracks(tc(i)).Frames(end)-jj;
                     else
                         r_latency_off(sco,lco)=900;
                     end
                    lco=lco+1;
                end
                if  ~isempty(on) & ~isempty(off)
                    if  on(1)<off(1)
                        if length(on)>length(off)
                            on(end)=[];
                        end
                        IRI=abs(off(1:end-1)-on(2:end));
                    else
                        IRI=abs(on-off(1:length(on)));
                    end
                    IRIall_CO2off(oc,1:length(IRI))=IRI;
                    oc=oc+1;
                end
            end
        end
        sco=sco+1;
    end
    IRIall(IRIall==0)=NaN;
    miri(F)=nanmean(nanmean(IRIall(:,1:3),2));
end

%%
IRIall(IRIall==0)=NaN;
IRIall_CO2off(IRIall_CO2off==0)=NaN;
IRIall(IRIall<7)=NaN;
IRIall_CO2off(IRIall_CO2off<7)=NaN;

save miri miri
save st_rev st_rev
save st_speed st_speed
save  IRIall IRIall
save  IRIall_CO2off IRIall_CO2off
save r_latency r_latency
save r_latency_off r_latency_off
IRIall(IRIall==0)=NaN;
IRIall_CO2off(IRIall_CO2off==0)=NaN;
IRIall(IRIall<7)=NaN;
IRIall_CO2off(IRIall_CO2off<7)=NaN;
%%     plot speed & reversals
figure
subtightplot(2,2,1)
%     speed_cent=meanspeed_bin-nanmedian(meanspeed_bin);
%     normspeed=smoothn((speed_cent/nanstd(speed_cent)),2);
sem=nanstd(meanspeed_bin)./sqrt(F);
shadedErrorBar2([],nanmean(meanspeed_bin),sem,'b',10);

set(0,'DefaultTextInterpreter','none');
title('speed')
ylim([0.025 0.21])
ylabel('mean speed (mm/s)')


subtightplot(2,2,2)
sem=nanstd(rev_bin)./sqrt(F);
shadedErrorBar2([],nanmean(rev_bin),sem,'k',10);
ylim([0 0.75])
ylabel('# reversals/worm')
title('reversals')

subtightplot(2,2,3)
% plot shift triggered speed average:
p=patch([ 5 65 65 5 ], [0.035 0.035 0.2 0.2], [1 1 1]);
set(p,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9])
hold on
shadedErrorBar2([],nanmean(st_speed(:,1:end)),nanstd(st_speed)/sqrt(bc-2),'b',0)
hold on
% plot([10 10],[0.1 0.15],'--','color',[0.5 0.9 0.6])
title('shift onset average speed')
ylabel('speed (mm/s)')
xlabel('time(s)')
ylim([0.025 0.21])

subtightplot(2,2,4)
p=patch([ 5 65 65 5 ], [0.035 0.035 0.2 0.2], [1 1 1]);
set(p,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9])
hold on

shadedErrorBar2([],nanmean(st_rev(:,1:end)),nanstd(st_rev)/sqrt(bc-2),'r',0)
ylim([0 0.75])
hold on
% plot([10 10],[0 0.15],'--','color',[0.5 0.9 0.6])
title('shift onset average rev')
ylabel('rev/animal')
xlabel('time(s)')
suptitle(dirname(cd))
%%
plotdir
cd(pdir(1).name)
saveas(gcf, 'shift_rev_state&speed.fig')
cd ..\

%% inter reversal intervals and latency
figure
subtightplot(1,2,1)
n=NaN;
for i=1:7
    n(i)=sum(~isnan(IRIall_CO2off(:,i)));
end
n(n==0)=NaN;
shadedErrorBar2([1:7],nanmedian(IRIall_CO2off(:,1:7))/10,nanstd(IRIall_CO2off(:,1:7)/10)./sqrt(n),'k',0)
hold on
scatter([1:7],nanmedian(IRIall_CO2off(:,1:7))/10, 'k')
for i=1:7
    n(i)=sum(~isnan(IRIall(:,i)));
end
n(n==0)=NaN;
shadedErrorBar2([1:7],nanmedian(IRIall(:,1:7))/10,nanstd(IRIall(:,1:7)/10)./sqrt(n),'g',0);
scatter([1:7],nanmedian(IRIall(:,1:7))/10, 'g');
ylim([2 24])
ylabel('IRI (sec)')
xlabel('#reversal')


subtightplot(1,2,2)
L=NaN(10,2);
L(1:size(r_latency,1),2)=nanmedian(r_latency./10,2);
L(1:size(r_latency_off,1),1)=nanmedian(r_latency_off./10,2);
boxplot(L,{'0%CO2';'1.5% CO2'})
ylabel('reversal latency after CO2 upshift(sec)')
ylim([0 100])
suptitle({dirname(cd);'inter reversal intervals & 1st reversal latency in CO2 upshifts/downshifts'})
%%
cd(pdir(1).name)
saveas(gcf, 'iri&onset.fig')
cd ..\


