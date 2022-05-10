% RIM state change (plateau to fall/plateau to rise)triggered ratio averages
%Julia Riedl 2017

clear
load ImagingDataAll

%  close

%parameters:
hz=10; % recording rate
time=8.1; %recording time in minutes
prewin=5*hz;%time before reversal rise_MoM
minp=10; %minimal plateau length;

revAll_rise=NaN(20,length(ImagingData));
rise_MoM=NaN(length(ImagingData),prewin*2);

plots=0;

cutoff_o=0.0031;

if plots==1
    figure
end

for F=1:length(ImagingData)
    
    if isfield(ImagingData{F}, 'cherry_rim')==1
        fall_bag=NaN(10,215);
        fall_rim=NaN(10,215);
        rise_bag=NaN(10,215);
        rise_rim=NaN(10,215);
        CO2_fall=NaN;
        CO2_rise=NaN(10,215);
        
        CC=1;
        CC2=1;
        
        endflag=0;
        
        %ratio bag:
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
        
        ratio_bag=gcamp_bag./cherry_bag;
        
        %R/R0: divide by the lower percentile
        sr=sort(ratio_bag);
        R0_bag=nanmean(sr(1:length(ratio_bag)/10));
        if  R0_bag<0.15
            R0_bag=nanmean(sr(1:length(ratio_bag)/8));
        end
        ratioR0=ratio_bag./ R0_bag;
        
        %ratio RIM
        cherry=medfilt1(ImagingData{F}.cherry_rim,10);
        gcamp=medfilt1(ImagingData{F}.gcamp_rim,10);
        ratio_rim=(gcamp/nanmean(gcamp))./(cherry/nanmean(cherry));
        ratio_rim=smoothn(ratio_rim,10);
        ratio_rim=ratio_rim-nanmedian(ratio_rim);
        
        %CO2
        if isfield(ImagingData{1}, 'stimulus')
            CO2=ImagingData{1}.stimulus;
            try
                CO2=[NaN(330,1) ; CO2(:,2)];
            catch
                CO2=[NaN(20,1) ; CO2(:,1)];
            end
        end
        
        ratioR0=ratioR0(1:length(CO2));
        if length(ratio_rim)>=length(CO2)
            ratio_rim=ratio_rim(1:length(CO2));
        else
            continue
        end
        
        %bleach correct
        [~, P]=linfit(1:length(ratio_rim(~isnan(ratio_rim))),ratio_rim(~isnan(ratio_rim)));
        rfit(:,1) = P(1)*(1:length(ratio_rim))+P(2);
        rfit=rfit+abs(min(rfit))+1;
        %bl=abs(min(0.7*rfit))+0.5;
        rim=(ratio_rim-rfit)+nanmedian(rfit);
        rim(1:250)=NaN; %remove data pre _stimulus onset
                %gcamp bleach correction:
        ratioR0_bc=(ratioR0-rfit)+nanmedian(rfit);
        
        %find rise and falls via derivative
        rf=0.0006;
        if nanstd(rim(1500:end))>0.8
            dRIM=derivReg(rim,0.0005,10);
            cutoff=cutoff_o*1.4;
            rf=0.0017;
        elseif nanstd(rim(1500:end))>0.37
            dRIM=derivReg(rim,0.00015,8);
            cutoff=cutoff_o*4;
            rf=0.0009;
        elseif nanstd(rim(1500:end))>0.25
            dRIM=derivReg(rim,0.0004,10);
            cutoff=cutoff_o*0.9;
        elseif nanstd(rim(1500:end))>0.16
            dRIM=derivReg(rim,0.00015,5);
            cutoff=cutoff_o*0.8;
            rf=0.0006;
        elseif nanstd(rim(1500:end))>0.1
            dRIM=derivReg(rim,0.0001,5);
            cutoff=cutoff_o*0.6;
        elseif nanstd(rim(1500:end))<0.1 & nanstd(rim(1500:end))>=0.064
            dRIM=derivReg(rim,0.000009,5);
            cutoff=cutoff_o*0.3;
        elseif nanstd(rim(1500:end))<0.064
            disp(['skipped F#' num2str(F)])
            continue
        end
        rc=1;
        while nanstd(dRIM(500:end))<cutoff/2 & rc<6
            rf =rf*0.5;
            dRIM=derivReg(rim,rf,7);
            rc=rc+1;
        end
        
        dRIM=medfilt1(dRIM,10);
        % rim=medfilt1(rim,20);
        %         rim=smoothn(rim,20);
        %         dRIM=[diff(rim); NaN];
        
        %start and endpoints of plateau phase:
        platvec=zeros(1,length(rim));
        platvec(dRIM>-cutoff*1 & dRIM<(cutoff*0.9))=1;
        
        tp=find(diff(platvec)~=0); %remove plateaus which are very short
        shortp=find(diff(tp)<minp);
        for i=1:length(shortp)
            platvec(tp(shortp(i)):tp(shortp(i)+1)+1)=1;
        end
        
        pstart=find(diff(platvec)==1)+1;
        pend=find(diff(platvec)==-1)-1;
        
        if length(pend)>1
            
            
            platvec(platvec==0)=2;
            if pend(1)<pstart(1) & pend(1)>prewin;
                endflag=pend(1);
                pend(1)=[];
            elseif pend(1)<pstart(1)
                pend(1)=[];
            end
            
            sc=0;
            low_start=NaN;
            while isnan(low_start) & sc<=3
                
                sc=sc+1;
                high_start=NaN;
                high_end=NaN;
                low_start=NaN;
                low_end=NaN;
                
                cc=1;
                cc2=1;
                
                if sc==1
                    mco=nanmean(rim);
                else
                    mco=mco*1.25;
                end
                
                for i=1:length(pstart) %go through plateau starts and asssign high or low state
                    
                    if pstart(i)>prewin & i==1 & nanmean(dRIM(pstart(i)-prewin:pstart(i)))>0.0015
                        high_start(cc)=pstart(i)+3;
                            if i<=length(pend)
                                high_end(cc)=pend(i)-5;
                                platvec(pstart(i):pend(i)+1)=3; %3=high state
                            end
                            cc=cc+1;
                            
                    elseif pstart(i)>prewin & nanmean(platvec(pstart(i)-19:pstart(i)))>1.8
                        if nanmean(dRIM(pstart(i)-prewin:pstart(i)))>0.0015
                            high_start(cc)=pstart(i);
                            if i<=length(pend)
                                high_end(cc)=pend(i)-5;
                                platvec(pstart(i):pend(i)+1)=3; %3=high state
                            end
                            cc=cc+1;
                             
                        elseif nanmean(dRIM(pstart(i)-prewin:pstart(i)))<-0.0015
                            
                            if i>1 & cc>1 & rim(pstart(i))>(0.2+nanstd(rim(1500:end))) % if an intermdediate up-plateau-end was alreday asssigned, overwrite
                                if i<=length(pend)
                                    high_end(cc-1)=pend(i);
                                else
                                    high_end(cc-1)=pstart(i);
                                end
                                
                            else
                                %low starts:
                                if rim(pstart(i))<mco
                                    low_start(cc2)=pstart(i)+5;
                                    
                                    if i<=length(pend)
                                        low_end(cc2)=pend(i);
                                        platvec(pstart(i):pend(i))=1;
                                    end
                                    cc2=cc2+1;
                                end
                            end
                        end
                        
                    end
                end
                
            end
            if endflag~=0 & low_start(1)<high_start(1)
                high_end=[endflag high_end];
            elseif endflag~=0 & low_start(1)>high_start(1) & ~isnan(low_end)
                low_end=[endflag low_end];
            end
            
            if low_start(end)>length(rim)-2
                low_start(end)=[];
            end
            
            %----stack rim and bag ratios around change points:----
            high_end(high_end>4690)=[];
            
            
            %rim fall episodes
            for i=1:length(high_end)
                if plots==1 &i==1
                    cla
                    scatter(1:length(rim),rim,10,platvec,'filled')
                    hold on
                    plot(ratioR0-1,'b')
                    plot(CO2*20,'color',[0.5 0.5 0.5])
                    title([F,nanstd(rim(1500:end))])  
                end
                %                 if i<=length(low_start)
                fall_bag(CC,1:size(fall_bag,2))=NaN;
                fall_rim(CC,1:size(fall_bag,2))=NaN;
                if high_end(i)>prewin
                    %find closest low start:
                    ds=high_end(i)-low_start;
                    if numel(find(ds<0))>1
                        ds(ds>0)=min(ds);
                    end
                    [~,idx]=min(abs(ds));
                    if ~isnan(low_start(idx)+10) & low_start(idx)+10<length(ratioR0_bc)
                        b=ratioR0_bc(high_end(i)-prewin:low_start(idx)+10);
                        if length(b)<215
                            if high_end(i)<1300 | nanmean(CO2(high_end(i)-10:low_start(idx)+10))>=0.022 | isnan(nanmean(diff(CO2(high_end(i)-2:low_start(idx)+10)))) ...
                                    | nanmean(CO2(high_end(i)-10:low_start(idx)+10))<=0.016
                                r=rim(high_end(i)-prewin:low_start(idx)+10);
                                if nanmean(diff(r(prewin:end-10)))<-cutoff*0.5
                                    %                                         if plots==1
                                    %                                             scatter(high_end(i),rim(high_end(i)),'o', 'markeredgecolor',[0.8 0.0 0.7],'SizeData',70);
                                    %                                             scatter(low_start(idx),rim(low_start(idx)),'p','filled', 'markerfacecolor',[0 1 0.8]);
                                    %                                         end
                                    fall_bag(CC,1:length(b))=b;
                                    fall_rim(CC,1:length(r))=r;
                                    CO2_fall(CC)=nanmean(diff(CO2(high_end(i)-10:high_end(i)+10)));
                                    CC=CC+1;
                                end
                            end
                        end
                    end
                end
                %                 end
            end
            
            %---RIM rise episodes:---
            if high_start(1)<low_end(1)
                high_start(1)=[];
            end
            
            for i=1:length(low_end)
                if ~isnan(high_start) & low_end(i)>350
                    rise_bag(CC2,1:size(rise_bag,2))=NaN;
                    rise_rim(CC2,1:size(rise_bag,2))=NaN;
                    if low_end(i)>prewin
                        %find closest low start:
                        ds=low_end(i)-high_start;
                        if numel(find(ds<0))>1
                            ds(ds>0)=min(ds);
                        end
                        [~,idx]=min(abs(ds));
                        if high_start(idx)>length(ratioR0_bc)-11
                            continue
                        end
                        b=ratioR0_bc(low_end(i)-prewin:high_start(idx)+10);
                        
                        if length(b)<215 % & nanmedian(diff(b(prewin:end)))<0.005
                            r=rim(low_end(i)-prewin:high_start(idx)+10);
                            if nanmean(diff(r(prewin:end-10)))>cutoff*0.5 & nanmean(diff(CO2(low_end(i)-10:high_start(idx)+1)))>=0 & nanmean(diff(CO2(low_end(i)-10:high_start(idx)+1)))<max(diff(CO2))/2% only CO2 stable or rising
                                
                                if plots==1
                                    scatter(low_end(i),rim(low_end(i)),'p','filled', 'markerfacecolor',[0 0.2 1],'SizeData',90);
                                    scatter(high_start(idx),rim(high_start(idx)),'p','filled', 'markerfacecolor',[1 0 0.5]);
                                pause(0.1)
                                end
                                rise_bag(CC2,1:length(b))=b;
                                rise_rim(CC2,1:length(r))=r;
                                CO2_rise(CC2,1:length(b))=CO2(low_end(i)-prewin:high_start(idx)+10);
                                CC2=CC2+1;
                                
                            end
                        end
                    end
                end
            end
            
        end
    end
    
    %normalise:
    norm_bag=NaN(size(rise_bag));
    for i=1:size(rise_bag,1)
        norm_bag(i,:)=rise_bag(i,:)./mean(rise_bag(i,prewin:prewin+5));
    end
    
    revAll_rise(1:size(rise_bag),F)=nanmean(norm_bag(:,prewin+hz:prewin+(4*hz)),2);
    rise_MoM(F,1:size(norm_bag,2))=nanmean(norm_bag,1);
    
end % end imaging files loop

rise_MoM(rise_MoM==0)=NaN;
%quantification:
nls_4=nanmean(rise_MoM(:,prewin+hz:prewin+(4*hz)),2);

%% MoM RIM High onset normalized reversals:
figure
sem=NaN;
tv=(-prewin:prewin*1.2)/hz;
for i=1:length(rise_MoM)
    sem(i)=nanstd(rise_MoM(:,i))/sqrt(numel(find(~isnan(rise_MoM(:,i)))));
    nn(i)=numel(find(~isnan(rise_MoM(:,i))));
end
h1(1)=shadedErrorBar2(tv,nanmean(rise_MoM(:,1:length(tv))),sem(1:length(tv)),'b',0);
hold on

plot([0 0],[min(nanmean(rise_MoM)-0.2) max(nanmean(rise_MoM)+0.2)],'r')
plot([-5 5],[1 1],'--k')

% set(gca,'XTick',[-rise_MoM/hz:2:rise_MoM/hz])
ylabel('norm dR')
xlabel('time (s)')
title([dirname(cd) ': BAG inhib, n= ' num2str(F)])
ylim([0.8 1.2])

plotdir
try
cd(pdir(1).name)
saveas(gca, ['RIM_trig_BAG_n_MoM.fig'])
cd ..\
catch
    saveas(gca, ['RIM_trig_BAG_n_MoM.fig'])
end




