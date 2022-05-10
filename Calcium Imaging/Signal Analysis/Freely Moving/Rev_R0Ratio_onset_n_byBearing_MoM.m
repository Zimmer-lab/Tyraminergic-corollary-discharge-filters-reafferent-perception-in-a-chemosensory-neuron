%reversal onset normalized R/R0 mean of mean
%(c) Julia  Riedl 2020
clear
load ImagingDataAll

%-- parameters:--
hz=30;
minlength=1;%sec
minlength=minlength*hz;

prewin=6; %in sec
prewin=prewin*hz;

postwin=6; %in sec
postwin=postwin*hz;

winsize=prewin+postwin;
%----

revAll_UG=NaN(20,length(ImagingData));
revAll_DG=NaN(20,length(ImagingData));

up_MoM=NaN(length(ImagingData),winsize);
down_MoM=NaN(length(ImagingData),winsize);


cc3=0;

for F=1:length(ImagingData)
    
    prerevUG=NaN(1,winsize);
    prerevDG=NaN(1,winsize);
    prerevUG_c=NaN(1,winsize);
    prerevDG_c=NaN(1,winsize);
    
    cc1=1;
    cc2=1;
    
    try
        cherry_bag=medfilt1(ImagingData{F}.cherry_bag,10);
        gcamp_bag=medfilt1(ImagingData{F}.gcamp_bag,10);
    catch
        cherry_bag=medfilt1(ImagingData{F}.cherry,10);
        gcamp_bag=medfilt1(ImagingData{F}.gcamp,10);
    end
    if min(gcamp_bag)<0
        gcamp_bag=gcamp_bag+abs(min(gcamp_bag));
    end
    
    ratio_bag=(gcamp_bag/nanmean(gcamp_bag))./(cherry_bag/nanmean(cherry_bag));
    
    %R/R0: divide by the lower percentile
    sr=sort(ratio_bag);
    R0_bag=nanmean(sr(1:length(ratio_bag)/10));
    if  R0_bag<0.15
        R0_bag=nanmean(sr(1:length(ratio_bag)/8));
    end
    ratioR0=ratio_bag./ R0_bag;
    
    %CO2
    sensIn=ImagingData{F}.CO2;
    nanIdx=find(isnan(sensIn));
%     SensNorm=sensIn-nanmedian(sensIn);
    CO2=smoothn(sensIn,10);
    CO2(nanIdx)=NaN;
    
    %reversals:
    if ~isfield(ImagingData{F},'RevFrames30hz')
        continue
    end
    RevON=ImagingData{F}.RevFrames30hz(1:2:end);
    RevEND=ImagingData{F}.RevFrames30hz(2:2:end);
    
    cc3=cc3+1; %number of animals
    
    for rev=1:length(RevEND)
        
        if RevON(rev)>length(ImagingData{F}.XY)-101 | RevON(rev)<2
            continue
        end
        if (RevEND(rev)-RevON(rev))<minlength
            continue
        end
        
        try
            XY=ImagingData{F}.XY(RevON(rev)-100:RevON(rev)+100,:);
        catch
            XY=ImagingData{F}.XY(RevON(rev)-1:RevON(rev)+100,:);
            XY=[NaN(99,2);XY];
        end
        
        [bearing]=getbearing_freelyMoving_global(XY(:,1),XY(:,2),1,10);
        
        %%%% reversals upgradient:
        if nanmean(bearing(100:end))>90
            
            if RevON(rev)>prewin & RevEND(rev)<length(ratioR0)
                
                prerevUG(cc1,:)=NaN(1,winsize);
                prerevUG_c(cc1,:)=NaN(1,winsize);
                revratio=ratioR0(RevON(rev)-prewin:RevEND(rev));
                prerevUG(cc1,1:length(revratio))=ratioR0(RevON(rev)-prewin:RevEND(rev));
                prerevUG_c(cc1,1:length(revratio))=CO2(RevON(rev)-prewin:RevEND(rev));
                cc1=cc1+1;
                % if reversal is too close to start of trace:
            elseif RevON(rev)<prewin
                prerevUG(cc1,:)=NaN(1,winsize);
                revratio=ratioR0(1:RevEND(rev));
                prerevUG(cc1,(prewin+1)-(RevON(rev)):prewin-(RevON(rev))+length(revratio))=revratio;
                cc1=cc1+1;
                %or to end
            elseif RevON(rev)+prewin>length(ratioR0)
                prerevUG(cc1,:)=NaN(1,winsize);
                revratio=ratioR0(RevON(rev)-winsize:end);
                prerevUG(cc1,1:length(revratio))=revratio;
                cc1=cc1+1;
            end
            
            
            %%%%%%% reversals downgradient:
        elseif nanmean(bearing(100:end))<90
            
            if RevON(rev)>prewin & RevEND(rev)<length(ratioR0)
                
                prerevDG(cc2,:)=NaN(1,winsize);
                revratio=ratioR0(RevON(rev)-prewin:RevEND(rev));
                prerevDG(cc2,1:length(revratio))=ratioR0(RevON(rev)-prewin:RevEND(rev));
                prerevDG_c(cc1,1:length(revratio))=CO2(RevON(rev)-prewin:RevEND(rev));
                cc2=cc2+1;
                % if reversal is too close to start of trace:
            elseif RevON(rev)<prewin
                
                prerevDG(cc2,:)=NaN(1,winsize);
                revratio=ratioR0(1:RevEND(rev));
                prerevDG(cc2,(prewin+1)-(RevON(rev)):prewin-(RevON(rev))+length(revratio))=revratio;
                cc2=cc2+1;
                %or to end
            elseif RevON(rev)+prewin>length(ratioR0)
                prerevDG(cc2,:)=NaN(1,winsize);
                revratio=ratioR0(RevON(rev)-winsize:end);
                prerevDG(cc1,1:length(revratio))=revratio;
                cc2=cc2+1;
            end
            
        end
        prerevUG=prerevUG(:,1:winsize);
        prerevDG=prerevDG(:,1:winsize);
        prerevUG_c=prerevUG_c(:,1:winsize);
        prerevDG_c=prerevDG_c(:,1:winsize);
    end % end reversal loop
    
    % normalize:
    %down gradient:
    prerev_n=NaN(size(prerevDG));
    
    for i=1:cc2-1
        mv=nanmean(prerevDG(i,prewin-5:prewin+5));
        prerev_n(i,:)=medfilt1(prerevDG(i,:),5)./mv;
    end
     down_MoM(F,1:size(prerev_n,2))=nanmean(prerev_n,1);
    
    %average each reversal over sec 1-3:
    ns1_3=nanmean(prerev_n(:,prewin+(hz):prewin+(hz*4)),2);
    revAll_DG(1:length(ns1_3),F)=ns1_3;
    
    % up-gradient:
    prerev_n=NaN(size(prerevUG));
    for i=1:cc1-1
        mv=nanmean(prerevUG(i,prewin-5:prewin+5));
        prerev_n(i,:)=medfilt1(prerevUG(i,:),5)./mv;
    end
    up_MoM(F,1:size(prerev_n,2))=nanmean(prerev_n,1);
    
    %average each reversal over sec 1-3:
    ns1_3=nanmean(prerev_n(:,prewin+(hz):prewin+(hz*4)),2);
    revAll_UG(1:length(ns1_3),F)=ns1_3;
    
   
    
end %end files

revAll_DG(revAll_DG==0)=NaN;
revAll_UG(revAll_UG==0)=NaN;

%% MoM boxplot up- vs. down-radient reversals:
figure
subtightplot(1,2,1)
sem=NaN;
for i=1:length(up_MoM)
    sem(i)=nanstd(up_MoM(:,i))/sqrt(numel(find(~isnan(up_MoM(:,i)))));
    nn(i)=numel(find(~isnan(up_MoM(:,i))));
end
h1(1)=shadedErrorBar2(-prewin/hz:0.0333334:postwin/hz,nanmean(up_MoM),sem,'r',0);
hold on
sem=NaN;
for i=1:length(down_MoM)
    sem(i)=nanstd(down_MoM(:,i))/sqrt(numel(find(~isnan(down_MoM(:,i)))));
    nn(i)=numel(find(~isnan(down_MoM(:,i))));
end
plot([-5 5],[1 1],'--k')
h1(2)=shadedErrorBar2(-prewin/hz:0.0333334:postwin/hz,nanmean(down_MoM),sem,'b',0);
plot([0 0],[min(nanmean(down_MoM)-0.2) max(nanmean(down_MoM)+0.2)],'r')

set(gca,'XTick',[-prewin/hz:2:postwin/hz])
xlabel('s')
ylabel('dF/F0')
title(['up/down gradient reversals (>100), n=' num2str(F)])
ylim([0.8 1.3])
xlim([-prewin/hz prewin/hz])


pdata=[nanmean(revAll_DG)' nanmean(revAll_UG)'];
subtightplot(1,2,2)
boxplot(pdata)
ylim([0.75 1.25])
set(gca,'xticklabels',{'down';'up'})
ylabel('mean dR/R0')
title([' reversals (>100), n=' num2str(length(ImagingData))])
suptitle(dirname(cd))

plotdir
cd(pdir(1).name)
saveas(gca, ['rev_ratioRO_n_MoM' dirname2(cd)])
cd ..\







