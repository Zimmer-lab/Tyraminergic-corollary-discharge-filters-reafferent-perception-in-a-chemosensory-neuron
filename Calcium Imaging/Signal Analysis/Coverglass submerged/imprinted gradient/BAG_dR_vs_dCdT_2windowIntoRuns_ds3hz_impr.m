%correlates ratio with dC during runs within specified window into run
clear % -except ImagingData
load ImagingDataAll

% close all
CC=0;
home=cd;
cc=1;

% win_ons=1; %minimal

bin=10; %downsample factor;
hz=30; %original recording rate;
minrl=3;% minimal reversal length
omega_co=1;

warning off

%----control plots?---
plots=0;

if plots==1
    figure
end
%-----------------------
wstep=5;
winl=wstep;

for win_ons=3:wstep:wstep*2  % time window into run  (in sec) )
    
    
    CC=CC+1;
    if CC==2
        winl=350;
    end
    %preallocate:
    mBAG=NaN(1, bin);
    mSens=NaN(1, bin);
    
    ratioAll=[];
    sensInAll=[];
    dCdTAll=[];
    dRdTAll=[];
    RunsAll=[];
    
    for F=1:length(ImagingData)
        
        if ~isfield(ImagingData{F}, 'RevFrames30hz')
            continue
        end
        %ratio: normalized and smoothed
        cherry=medfilt1(ImagingData{F}.cherry_bag,5);
        gcamp=medfilt1(ImagingData{F}.gcamp_bag,5);
        nanIdx=find(isnan(cherry));
        
        
        %         ratio_o=(gcamp/nanmean(gcamp))./(cherry/nanmean(cherry));
        eq_ratio=(gcamp)./(cherry);
        ratio_o=ratio_R0(eq_ratio);
        
        %CO2:
        sensIn=ImagingData{F}.CO2;
        
        
        %downsample to X hz:
        ratio=bindata(ratio_o, bin);
        nanidx=find(isnan(ratio));
        ratio=smoothn(ratio,10);
        ratio(nanidx)=NaN;
        sensIn=bindata(sensIn, bin);
        nanidx=find(isnan(sensIn));
        sensIn=smoothn(sensIn,10);
        sensIn(nanidx)=NaN;
        
        %reversals:
        RevON=round(ImagingData{F}.RevFrames30hz(1:2:end)/bin);
        RevEND=round(ImagingData{F}.RevFrames30hz(2:2:end)/bin);
        RevON(RevON<1)=1;
        RevEND(RevEND>length(ratio))=[];
        RevON(RevON>length(ratio))=[];
        
        
        
        %(1)take only  winl s of runs after reversal:
        Runs=zeros(length(ratio),1);
        if length(RevEND)<length(RevON) %Reversal lengths
            revl=abs(RevON(1:end-1)-RevEND);
        else
            revl=abs(RevON-RevEND);
        end
        
        for i=1:length(RevEND) %set run indices within timewindow to 1 after reversal...
            
            if i<length(RevEND)
                if   revl(i)>=minrl*(hz/bin) &  RevEND(i+1)>RevEND(i)+((win_ons+winl)*hz/bin)
                    Runs(RevEND(i)+(win_ons*(hz/bin)):RevEND(i)+((win_ons+winl)*hz/bin))=1;
                elseif revl(i)>=minrl*(hz/bin)
                    Runs(RevEND(i)+(win_ons*(hz/bin)):RevON(i+1))=1;
                end
            elseif revl(i)>=minrl*(hz/bin) & length(ratio) > RevEND(i)+((win_ons+winl)*hz/bin)
                
                Runs(RevEND(i)+(win_ons*(hz/bin)):RevEND(i)+((win_ons+winl)*hz/bin))=1;
                
            end
        end
        Runs=Runs(1:length(ratio));
        
        % plot: proofreading plot
        if plots==1
            cla
            hold on
            plot([1:length(ratio)],ratio,'k')
            h=scatter([1:length(ratio)],ratio,5,Runs);
            scatter(RevEND,ratio(RevEND),'r')
            scatter(RevON,ratio(RevON),'pg')
            plot((sensIn/15)-0.05)
            set(gca,'Xgrid','on');
            %         caxis([0 2])
            title(F)
        end
        
        
        try
            sensIn=sensIn(1:length(ratio));
            
        catch
            sensIn=[sensIn NaN(1,length(ratio)-length(sensIn))];
            sensIn=sensIn(1:length(ratio));
        end
        
        ratioAll=[ratioAll; ratio(:)];
        RunsAll=[RunsAll; Runs(:)];
        %derivatives:
        sensInAll=[sensInAll sensIn] ;
        dRdT=[diff(ratio)'; NaN];
        dRdTAll=[dRdTAll ; dRdT];
        
        dCdT=[diff(sensIn) NaN];
        dCdTAll=[dCdTAll dCdT];
        
        
        %  end
        if length(dRdTAll)~=length(dCdTAll)
            disp('inequal vector length')
            break
        end
    end
    
    
    %% %%%% bin & plot %%%%%%%
    
    subratio=dRdTAll;
    subdC=dCdTAll;
    gi=find(RunsAll==1); %  timepoints within run window only
    subratio=dRdTAll(gi);
    subdC=dCdTAll(gi);
    
    mBAG=NaN(5000,19);
    mSens=NaN(5000,19);
    mDC=NaN(5000,19);
    mDR=NaN(5000,19);
    binsize=0.0032;
    maxDC=nanmean(dCdTAll)+(3*nanstd(dCdTAll)); %0.019;
    minDC= -(nanmean(dCdTAll)+(3*nanstd(dCdTAll))); %-0.019;
    binsize=(maxDC-minDC)/10;
    disp(['data range 4 SD: +- ' num2str(nanmean(dCdTAll)+(3*nanstd(dCdTAll))) '; edge=' num2str(minDC) '; binsize=' num2str(binsize)])
    
    
    cc=1;
    for i=minDC :binsize:(maxDC-binsize)
        
        bin_idx= find(subdC>i & subdC<=i+binsize);
        if i >= maxDC-(binsize*1.9)
            bin_idx= find(subdC>i);
        elseif i== minDC
            bin_idx= find(subdC<i+binsize);
        end
        
        
        
        if length(bin_idx)>length(find(~isnan(subratio)))/300
            
            bV=ratioAll(bin_idx);
            sV=sensInAll(bin_idx);
            dCV=subdC(bin_idx);
            dRV=subratio(bin_idx);
            
            rd(cc)=length(find(bV(~isnan(bV))));
            
            mBAG(1:length(bV),cc)=(bV);
            mSens(1:length(sV),cc)=(sV);
            mDC(1:length(dCV),cc)=(dCV);
            mDR(1:length(dRV),cc)=(dRV);
            
        end
        cc=cc+1;
    end
    %% %%%% plot %%%%%%%
    mDR=mDR(:,1:cc-1);
    mDC=mDC(:,1:cc-1);
    % % boxplot/scatter:
    if CC==1
        broadfig;
    end
    subtightplot(1,2,CC)
    boxplot(mDR*(hz/bin),'symbol','.w','whisker',0.72)
    h=findobj(gca,'tag','Outliers');
    delete(h)
    xl=nanmedian(mDC);
    set(gca,'XTick',[1:1:length(xl)])
    set(gca,'XTickLabel',round(xl(1:1:end)*(hz/bin*1000))/1000);
    if CC==1
        ylim([-0.55 0.45])
    elseif CC==2
        ylim([-0.55 0.45])
    end
    
    %     yticks([-0.15:0.05:-0.15])
    xlabel(' CO2 change (%/sec)')
    ylabel(' dR/sec')
    hold on
    plot([0 21],[0 0],'--k')
    
    %correlations:
    mp=floor((cc-1)/2);
    [r(CC,1),p(CC,1)]=corr(nanmean(mDR(:,1:mp))',nanmean(mDC(:,1:mp))');
    [r(CC,2),p(CC,2)]=corr(nanmean(mDR(:,mp+1:mp*2))',nanmean(mDC(:,mp+1:mp*2))');
    [r(CC,3),p(CC,3)]=corr(nanmean(mDR)',nanmean(mDC)');
    if isnan(r(CC,2))
        [r(CC,2),p(CC,2)]=corr(nanmean(mDR(:,mp+1:(mp*2)-1))',nanmean(mDC(:,mp+1:(mp*2)-1))');
        if isnan(r(CC,2))
            [r(CC,2),p(CC,2)]=corr(nanmean(mDR(:,mp+1:(mp*2)-2))',nanmean(mDC(:,mp+1:(mp*2)-2))');
        end
        
    end
    if isnan(r(CC,1))
        [r(CC,1),p(CC,1)]=corr(nanmean(mDR(:,2:mp))',nanmean(mDC(:,2:mp))');
    end
    gi=find(~isnan(nanmean(mDC)));
    gi2=find(~isnan(nanmean(mDR)));
    gi=intersect(gi,gi2);
    [r(CC,3),p(CC,3)]=corr(nanmean(mDR(:,gi))',nanmean(mDC(:,gi))');
    
    r=roundd(r,2);
    p=roundd(p,3);
    
    plot([1 mp],[0.7 0.7],'-k')
    text(0.9,0.66,['r=' num2str(r(CC,1)) '; p=' num2str(p(CC,1))])
    plot([mp+1 mp*2],[0.7 0.7],'-k')
    text(mp+1,0.74,['r=' num2str(r(CC,2)) '; p=' num2str(p(CC,2))])
    plot([1 mp*2],[-1 -1],'-k')
    text(mp-3,-1.05,['r=' num2str(r(CC,3)) '; p=' num2str(p(CC,3))])
    
    title([ num2str(round(win_ons)) '-'  num2str(round(win_ons+winl)) ' s into run,n=' num2str(length(subratio)) 'dp'])
    suptitle([dirname(cd) ': '  num2str(hz/bin) 'hz, min rev l=' num2str(minrl) 's'])
    %     %increase window start point:
    %     win_ons=win_ons+wstep;
    
end
%save
pdir=plotdir;
cd(pdir(1).name)
saveas(gca, ['coding_2xrun_win' num2str(wstep) 'win_minrl_' num2str(minrl) '_' num2str(round(hz/bin)) 'hz.fig'])
cd ..\


