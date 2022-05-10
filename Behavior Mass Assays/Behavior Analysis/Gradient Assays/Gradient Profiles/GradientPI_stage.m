%(c) Julia Riedl 2016
%calculates chemotaxis preference PI over timepoints  and genotype folders and plots

clear
figure('visible','on')

%--- edit search strings if necessary for your genotype or whatever your folders naming is---
directory=dir('*_0*');
genotypes(1)={'*_0*'};

home=cd;

%----- to edit----
timebin=10; % How many time points??? 10=every minute for 10 hz recordings
TP=15; %timepoint for comparing indices between conditions

CI2=NaN(10,length(directory));
cc=1;
for dd=1:length(directory)
    clearvars -except dd cc home directory genotypes CI2 H TP timebin
    disp(directory(dd,1).name)
    if ~isdir(directory(dd,1).name)
        continue
    end
    cd(directory(dd,1).name)
    
    sd=dir('*analysis*');
    if isempty(sd)
        cd ..\
        continue
    end
    cd(sd(1).name)
    
    for cond=1%:2
        
        
        
        %     if mod(cc,2)==0
        %          cd('gradient')
        %     else
        %         cd('control')
        %
        %     end
        
        
        
        % average over e experiments for t timebins:
        file=dir('HMnorm_batch*');
        load(file(1).name);
        % clearvars -except directory home dd HMnorm_batch SpeedM_batch angSpeedM_batch cc genotypes
        c=1;
        
        
        
        for t=1:(length(HMnorm_batch{1,1})/timebin)-1
            % it starts with a structure containing the normalized (percentage of
            % animals at a given timepoint in a givenbin of the arena)
            mHM=NaN(size(HMnorm_batch{1,1}(:,:,1)));
            
            for e= 1:length(HMnorm_batch)
                mHM1=squeeze(nanmean(HMnorm_batch{1,e}(:,:,c:c+timebin),3)); %averaging over 100 sec
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
            disp(['last bin: ' num2str(length(HMnorm_batch{1})-c) '...length of stack:' num2str(length(HMnorm_batch{1}))])
            %each cell: averaged data for e experiments and one timebin
            HMoE{t}=mHM;
            c=c+timebin;
        end
        
        %% plot mean ocupancy for sectors over 4 time bins
        J=jet(length(HMoE));
        nd=(cd);
        d= strfind(cd, '\');
        name=nd(d(end)+1:end-9)
        title (name)
        % legend('sem','1-7.5','','7.5-15','','15-22.5','','22.5-30');
        % saveas(gcf, strcat('distribution_',name))
        
        
        %% plot Gradient Index:
        J=jet(4);
        jj=1;
        
        for i=1:length(HMoE)%length(HMoE)=timebins, in each cell 3rd dimension is experiments, so each i is one timepoint
            sHM=squeeze(nansum(HMoE{i},1)); %summing up all animals along the short axis of arena, keeping n experiments
            md=size(sHM,1);
            F1=nansum(sHM(ceil(md/2):md,:),1);%this vector has length of number of experiments
            F2=nansum(sHM(1:floor(md/2),:),1);
            CI(jj,1:size(sHM,2))=(F1-F2)./(F2+F1); %calculate PI
            jj=jj+1;
        end
        
        mCI=nanmean(CI,2);
        sem=nanstd(CI')/sqrt(e);
        subplot(1,2,1)
        hold on
        H(cc)=shadedErrorBar2([1:20/(length(mCI)+1):20],mCI,sem,{'-or','Color',J(cc,:)},0);
        xlabel(' time (min)')
        ylabel ('PI')
        % boxplot(PI');
        hold on
        genotypes{cc}= (name);
        cc=cc+1;
        CI2(1:length(CI(end,:)),cc-1)=CI(TP,:);
        save CI CI
        
        cd ..\
        
    end
    
    cd(home)
    
    
end
%% plot

plot([0 25],[0.0 0.0],'--k')

try
    legend(H(1:length(H)),genotypes);
end
ylim([-1 0.95])
xlim([0 24])
% title('CI over time')
%%
subplot(1,2,2)
cla
%notBoxplot(CI2,[]'whisker',1,'labels',{'N2 ctrl','N2 gradient','flp-1,nlp12 ctrl',' flp-1,nlp12 grad'})
CI2(CI2==0)=NaN;
notBoxPlot(CI2(:,1:cc-1))

try
    set(gca,'Xticklabel',{directory(1).name(1:10),directory(2).name(1:10),directory(3).name(1:12)},'FontSize',8)
catch
    try
        set(gca,'Xticklabel',{directory(1).name(1:15),directory(2).name(1:15)},'FontSize',8)
    catch
        set(gca,'Xticklabel',{directory(1).name,directory(2).name},'FontSize',8)
    end
end
ylim([-1 0.5])
title(['timepoints compared:' num2str(round(TP*(20/length(CI))))])

for i=1:2
    try
        a=CI2(:,i);
        a(isnan(a))=[];
        b=CI2(:,i+1);
        b=b(~isnan(b));
        [h(i) p(i)]=ttest2JR(a,b)
        text(i+0.2,0.15,['p=' num2str(round(p(i)*100000)/100000)])
    end
end


