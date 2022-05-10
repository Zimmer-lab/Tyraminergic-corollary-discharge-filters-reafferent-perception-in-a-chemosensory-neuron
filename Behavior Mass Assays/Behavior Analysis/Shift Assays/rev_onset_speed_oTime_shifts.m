
% calculates % of tracks on reversal of all tracks over time for various tracks files
clearvars -except out
% close
files= dir('*als.mat');
hz=10;
CO2on=(2450:3000:16500);
% CO2on=(1250:600:11500);
pause=240;
up=60;


meanspeed_bin=NaN;
norm_revs=NaN;
rev_bin=NaN;

bc=1;

for F=1:length(files)
    tic
    disp(['loading..' ,files(F).name]);
    load(files(F).name);
    toc
    %preallocate
    cc=1;
    speed_matrix=NaN(300,30*60*hz);
    rev_matrix=NaN(300,30*60*hz);
    animals=NaN(300,30*60*hz);
    
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
        RevV=zeros(1,length(rev_matrix));
        Rev=(Tracks(1,T).polishedReversals);
        frames=Tracks(T).Frames;
        
        if ~isempty(Rev) & max(d)>=20;% only 1 sec long reversals
            gi=find(Rev(:,4)==1);
            Rev=Rev(gi,:);
            RevF=frames(Rev(:,1))';
            for i=1:size(Rev,1)
                RevV(RevF(i))=1;
            end
        end
        
        %correct for simulus onset delay on different setups:
        if F>2 & F<7 & strcmp(dirname2(cd),'tdc-1') & strcmp(dirname(cd),'WT')
            RevV=[RevV(50:end) NaN(1,49)];
        end
         if F>2 & F<=8 & strcmp(dirname2(cd),'tdc-1') & strcmp(dirname(cd),'tdc-1')
            RevV=[RevV(50:end) NaN(1,49)];
        end
        if F==3 & strcmp(dirname2(cd),'tdc-1') & strcmp(dirname(cd),'tdc-1')
            RevV=[RevV(120:end) NaN(1,119)];
        end
        
        if F==2 & strcmp(dirname2(cd),'tdc-1') & strcmp(dirname(cd),'tdc-1')
            RevV=[RevV(25:end) NaN(1,24)];
        end
        
        rev_matrix(cc,1:length(RevV))=RevV;
        animals(cc,frames)=1;
        
        %remove reversals from speed
        ri=find(rev_matrix(cc,:)==1);
        speed_matrix(cc,ri)=NaN;
        if ~isempty(ri) & ri(1)==1
            ri(1)=2
        end
        speed_matrix(cc,ri-1)=NaN;
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
        rev_bin(F,c)=ceil(nansum(all_revs(i:i+bin)))/nanmeanJ(num_animals(i:i+bin));
        c=c+1;
    end
      
     % shift-triggererd_averages:

     for jj=CO2on/bin-5;
             st_speed(bc,:)=meanspeed_bin(F,jj:jj+up*2);
             st_rev(bc,:)= rev_bin(F,jj:jj+up*2);
             bc=bc+1;
     end

     
end

save st_rev st_rev
save st_speed st_speed
save rev_bin rev_bin

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
ylim([0 0.35])
ylabel('# reversals/worm/s')
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
ylim([0 0.35])
hold on
% plot([10 10],[0 0.15],'--','color',[0.5 0.9 0.6])
title('shift onset average rev')
ylabel('rev/animal')
xlabel('time(s)')
suptitle(dirname(cd))
%%
plotdir
cd(pdir(1).name)
saveas(gcf, 'shift_rev_onset&speed.fig')
cd ..\


