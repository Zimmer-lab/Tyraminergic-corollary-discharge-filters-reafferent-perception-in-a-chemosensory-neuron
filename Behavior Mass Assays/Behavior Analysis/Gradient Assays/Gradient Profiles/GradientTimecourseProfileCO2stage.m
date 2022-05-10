%(c) Julia Riedl 2017
% plots mean X dimension profiles of worm residency in % of all animals per
% bin

figure('visible','on')
set(0,'DefaultTextInterpreter','none');
clear
home=cd;
genotypes=cell(8);
genotypes(1:8)={''};

% to edit:
tb=12; %averaging timebin
nP=2; % number of profiles

 
% average over e experiments for t timebins:
file=dir('HMnorm_batch*');
load(file(1).name);
% load HMnorm_batch
dim=NaN(1,3);
for i=1:length(HMnorm_batch)
dim(i,:)=size(HMnorm_batch{i});
end
dim=max(dim);
% clearvars -except directory home dd HMnorm_batch SpeedM_batch angSpeedM_batch cc genotypes
c=1;

HMoE={NaN};
t=1;
for c=1:tb:length(HMnorm_batch{1})-tb
    % it starts with a structure containing the normalized (percentage of
    % animals at a given timepoint in a givenbin of the arena)
    mHM=NaN(dim(1),dim(2));
    
    for e= 1:length(HMnorm_batch)
        mHM1=squeeze(nanmeanJ(HMnorm_batch{1,e}(:,:,c:c+tb),3)); %averaging over tb 
        if (length(HMnorm_batch{1})-(c+tb))<tb
            mHM1=squeeze(nanmeanJ(HMnorm_batch{1,e}(:,:,c:end),3));
        end
        mHM(1:size(mHM1,1),1:size(mHM1,2),e)=mHM1;
    end

  %each cell: averaged data for e experiments and one timebin
  bi=find (mHM==Inf);
  mHM(bi)=NaN;
    HMoE{t}=mHM;
 
    t=t+1;
end
disp(['last bin: ' num2str(length(HMnorm_batch{1})-c) '...length of stack:' num2str(length(HMnorm_batch{1}))])
 nd=(cd);
 d= strfind(cd, '\');
name=nd(d(end)+1:end);


cd (home)

%% plot mean ocupancy for sectors over nP time bins
%now the data are organized by timepoint: each cell contains n arenas for
%one time point (values: %!!):
J=winter(length(HMoE));
cc=1;

for i=1:floor(length(HMoE)/nP):length(HMoE)
    
    sHM=nansum(HMoE{i},1); %!! summing along one dimension of arena(previous error:taking mean instead of the sum
    sHM=squeeze(sHM);
    if sum(isnan(sHM(5,:)))<length(HMnorm_batch)/2
        meanOE=nanmeanJ(sHM,2);% averaging over experiments
        sorted=sort(meanOE);
        semOE=nanstd(sHM')/sqrt(e);
        hold on
        
        l=shadedErrorBar2([],meanOE,semOE,J(i,:),0);
        sum(meanOE);
    
    if ~isempty(l) & ~isnan(nanmean(meanOE))
        try
        lh(cc)=l.mainLine;
        catch
            lh(cc)=l;
        end          
        lg{cc}=num2str(round((i/length(HMoE))*20));
        set(gca,'XTick',[1:3:length(meanOE)])
        set(gca,'XTickLabel',(round([1:(2.9/length(meanOE)*3.2):4]*10))/10);
        cc=cc+1;
    end
    end
end
disp(i)

ylim([0 16])
% xlim([1 30])
title([name ' , n=' num2str(size(sHM,2))])
legend(lh,lg)
xlabel('% CO2')
ylabel (' % worms')
set(gca, 'FontSize', 11)

