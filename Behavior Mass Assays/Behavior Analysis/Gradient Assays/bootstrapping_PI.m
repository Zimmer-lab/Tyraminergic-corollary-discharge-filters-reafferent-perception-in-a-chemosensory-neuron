% bootsrapping chemotaxis pref. index:

% draw 2 random samples from the pooled data( control and gradient e.g.)and
%calculate PI, iterate X times
clear

tic
home=cd;
pvalue=NaN;
datapool=[];
nn=NaN(1,2);
%-----paramteres to edit----
iterations=12000;  %how many bootstrapping iterations
timebin=20;  % How many time points averaged??? 10=every minute for 10 hz recordings
TP=0.5;    %=time point to compare(fraction of 20 min experiment)
cp=-0.011;   %input('control PI?');
ap=-0.593;   %input('gradient PI?');

for F=1:2
    
    %load data (control and gradient) and pool
    if F==1
        ddir=dir('0-5%*');
        cd(ddir(1).name)
        sdir1=dir('WT*');
        disp(F)
        disp(sdir1(1).name)
        cd(sdir1(1).name)
        
    elseif F==2
        ddir=dir('0 %C*');
        cd(ddir(1).name)
        sdir1=dir('Lite*');
        disp(F)
        disp(sdir1(1).name)
        cd(sdir1(1).name)
    end
    
    sdir2=dir('*profile_analysis');
    cd(sdir2(1).name)
    load HM_batch.mat
    c=1;
    for t=1:(length(HM_batch{1,1})/timebin)-1
        % it starts with a structure containing the normalized (percentage of
        % animals at a given timepoint in a givenbin of the arena)
        mHM=NaN(size(HM_batch{1,1}(:,:,1)));
        
        for e= 1:length(HM_batch)
            mHM1=squeeze(nanmean(HM_batch{1,e}(:,:,c:c+timebin),3)); %averaging over time
            if size(mHM(:,:,1))==size (mHM1)
                mHM(:,:,e)=mHM1;
            else
                if length (mHM(:,:,1))<size (mHM1)
                    mHM(:,:,e)=mHM1(1:size(mHM,1),1:size(mHM,2));
                elseif length (mHM(:,:,1))>size (mHM1)
                    mHM(1:size(mHM1,1),1:size(mHM1,2),e)=mHM1;
                end
            end
        end
        %each cell: averaged data for e experiments and one timebin
        HMoE{t}=mHM;
        c=c+timebin;
    end
    md=length(HMoE);
    i=round(length(HMoE)*TP); %length(HMoE)=timebins, in each cell 3rd dimension is experiments, so each i is one timepoint
    sHM=squeeze(nansum(HMoE{i},1)); %summing up all animals along the short axis of arena, keeping n experiments
    nn(F)=size(sHM,2);
    datapool=[datapool;sHM'];
    
    cd(home)
    
end
if nn(1)>nn(2)
    datapool=datapool(nn(1)-nn(2)+1:end,:);
end
PI=NaN;
for i= 1:iterations % calculate PI "iterations" times
    
    S1=NaN(1,1);
    ds=size(datapool,2);
    
    for bin= 1:ds %draw 1 sample for each position bin
        S1(bin) = randsample(datapool(:,bin),1,1);
    end
%     S1=S1/sum(S1);%normalize
    F1=nansum(S1(ceil(ds/2):ds));%this vector has length of number of experiments
    F2=nansum(S1(1:floor(ds/2)));
    PI(i)=(F1-F2)./(F2+F1); %calculate PI
   
end

toc
%% plot histogram:
%variance of distribution:
rf=100;%rounding factor
[b,v]=normhist(PI,20,1);
phi=roundd(mean(PI)+(abs(std(PI))*1.96),2); % 1.96 SD cover 95% of the data
phi2=round(mean(PI)-(abs(std(PI))*1.96),2);

hold on
name=dirname2(cd);
title([ name ' # of iterations: ' num2str(iterations)])
plot([phi phi ],[0 0.12],'--k')
plot([phi2 phi2 ],[0 0.12],'--k')
text(phi,0.10,['2sig=' num2str(round(phi*1000)/1000)])

if ap~=0
    l(1)=plot([ap ap ],[0 0.12],'r');
    text(ap,0.12,['gradient PI=' num2str(round(ap*1000)/1000)])
end

if cp~=0
    %%
    l(2)=plot([cp cp ],[0 0.12],'color', [0.5 0.5 0.5]);
    text(cp,0.12,['control PI=' num2str(round(cp*1000)/1000)])
end

%% p values:

if ap>mean(PI)
    p_gradient=length(find(PI>ap))/length(find(PI<ap));
else
    p_gradient=length(find(PI<ap))/length(find(PI>ap));
end

if cp>mean(PI)
    p_control=length(find(PI>cp))/length(find(PI<cp));
else
    p_control=length(find(PI<cp))/length(find(PI>cp));
end

legend(l,['p=' num2str(p_gradient)],['p=' num2str(p_control)])
