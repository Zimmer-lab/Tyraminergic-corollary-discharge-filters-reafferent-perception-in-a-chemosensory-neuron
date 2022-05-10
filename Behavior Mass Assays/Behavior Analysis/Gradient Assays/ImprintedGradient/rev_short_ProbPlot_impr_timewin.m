%% reversal frequency modulation by bearing in specific time window of behavioral experiment
%(c) Julia Riedl 2019-21
% use on runinfo files
clear
warning off;
CC=1;
CM= jet(6);
bintype=1; % bin for dC/dT(2) or bearing(1)
delay=1;

%plot?
plotting=1;
files =dir('*rev_*info*');
load(files(1).name);
disp(files(1).name);

%parameters:
if bintype==1
    first=0;
    last=180;
    binsize=20;
else
    first=-0.13;
    last=0.13;
    binsize=0.02;
end
hz=10; %recording rate
mintime=5;


for win=25 %maximal time into exp (minutes)
    disp(['<' num2str(win) ' minutes'])
    % for each experiment: put all run data into one vector:
    mB= NaN(55,40);
    mC= NaN(55,40);
    mT= NaN(55,40);
    revN= NaN(55,40);
    revNnorm=NaN(55,40);
    c=0;
    
    for F=1:length(bearing)
        
        if ~isempty(bearing{F})
            
            c=c+1;
            bearingAll=[];
            dCdtAll=[];
            SpeedAll=[];
            revAll=[];
            
            
            for i=1:length(bearing{F})
                if ~isempty(bearing{F}{i}) & seconds{F}{i}(1)<win*60 & seconds{F}{i}(1)>mintime*60
                    bv=bearing{F}{i};
                    %                 cp=sensPath{F}{i};
                    %                 dC=[NaN NaN NaN cp(4:end)-cp(1:end-3) ];
                    sp=speed_ctr{F}{i};
                    rev=reversals{F}{i};
                    
                    bearingAll=cat(2,bearingAll,bv);
                    %                 dCdtAll=cat(2,dCdtAll,dC);
                    SpeedAll=cat(2,SpeedAll,sp);
                    revAll=cat(2,revAll,rev);
                end
            end
            
            
            bv=round(bearingAll*10)/10;
            
            
            %%%%%bin:
            cc=1;
            if ~isempty(bv)
                
                for i=first :binsize:(last-binsize)
                    
                    if bintype==1
                        bin_idx= find(bv>i & bv<=i+binsize);
                    else
                        bin_idx= find(dCdtAll>i & dCdtAll<=i+binsize);
                    end
                    bin_idx(bin_idx>length(bv)-delay)=[];
                    
                    if ~isempty(bin_idx)
                        %         dC=dCdtAll(bin_idx);
                        revV=revAll(bin_idx+delay);
                        binB=bv(bin_idx);
                        
                        
                        mB(c,cc)=nanmean(binB);
                        %         mC(c,cc)=nanmean(dC);
                        revN(c,cc)=nansum(revV);
                        revNnorm(c,cc)=nansum(revV)/length(bin_idx);
                    end
                    
                    cc=cc+1;
                    
                end
                
            end
        end
    end
    
    revNnorm=revNnorm(1:c,1:cc-1);
    % mC=mC(1:c,1:cc-1);
    revN=revN(1:c,1:cc-1);
    mB=mB(1:c,1:cc-1);
    
    %quantification: fold change low bearing to high bearing
    lb=nanmean(revNnorm(:,1:2),2)*3.3*60;
    hb=nanmean(revNnorm(:,end-1:end),2)*3.3*60;
    rq=hb./lb;
    [~, p]=ttest(revNnorm(:,1),revNnorm(:,end))
    save rq rq
    
    
    %% plot:(1) reversal rate binned by bearing
    name=[dirname2(cd) ': ' dirname(cd)];
    
    if plotting==1
        if CC==1
            figure
        end
        hold on
        sem=nanstd(revNnorm*3.3*60,1)/sqrt(c);
        hold on
        h(CC)=plot(nanmean(mB),nanmean(revNnorm*3.3*60,1),'color',CM(CC,:));
        errorb(nanmean(mB),nanmean(revNnorm*3.3*60,1),sem);
        %         h(CC)=shadedErrorBar2(round(nanmedian(mB,1)*1)/1,nanmean(revNnorm*3.3*60,1),sem,{'-or','color',CM(CC,:)},0);
        xlabel('bearing')
        ylabel('reversal frequency (rev/min)')
        %ylim([0.012 0.025])
        title ([name ])
        
    end
    leg{CC}=[num2str(mintime) '-' num2str(win) 'min'];
    CC=CC+1;
end % end win loop
%%
legend(leg)
ylim([0.4 1.9])
%save:
plotdir=dir('*plots*');
if isempty(plotdir)
    name=dirname(cd);
    mkdir([name ' plots'])
    plotdir=dir('*plots*');
end

cd(plotdir(1).name)
saveas(gca, ['reversal_freq_modulation_' num2str(mintime) '-' num2str(win) '.fig'])
cd ..\
% save turningbias mT
y=chirp(1:0.001:1.5,30);
sound(y)


%% binary plot:
figure
scatter([1 2],[nanmean(lb) nanmean(hb)],'b')
hold on
errorbar([1 2],[nanmean(lb) nanmean(hb)],[std(lb)/sqrt(length(lb)) std(hb)/sqrt(length(hb))], 'b')
ylim([0.4 1.9])
set(gca, 'xtick',[1 2])
set(gca,'xticklabel',{'bearing<30'; 'bearing<150'})
xlim([0.5 2.5])
title(dirname(cd))
cd(plotdir(1).name)
saveas(gca, ['reversal_freq_modulation_binary' num2str(mintime) '-' num2str(win) '.fig'])
cd ..\




