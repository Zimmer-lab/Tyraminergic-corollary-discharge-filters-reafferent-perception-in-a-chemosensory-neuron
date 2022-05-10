%correlates ratio with dC during runs within specified window into run
clear % -except ImagingData
load ImagingDataAll

% close all
CC=0;
home=cd;
cc=1;

win_ons=1; %time into run window starts (in sec)

bin=15; %downsample factor;
hz=30; %original recording rate;
minrl=2;% minimal reversal length
omega_co=1;

warning off

mBAG=NaN(2000,12);
mSens=NaN(2000,12);
mDC=NaN(2000,12);
mDR=NaN(2000,12);

%----control plots?---
plots=0;

if plots==1
    figure
end
%-----------------------
wstep=2;
winl=2;

for win_ons=0:wstep:wstep*5  %WINDOW WIDTH (IN S )
    
    
    CC=CC+1;
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
        
%         eq_ratio=(gcamp)./(cherry);
%         ratio_o=ratio_R0(eq_ratio);
ratio_o=(gcamp/nanmean(gcamp))./(cherry/nanmean(cherry));
        
        %CO2:
        sensIn=ImagingData{F}.CO2;
        
        
        %downsample to 3 hz:
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
                if   revl(i)>minrl*(hz/bin) &  RevEND(i+1)>RevEND(i)+((win_ons+winl)*hz/bin)
                    Runs(RevEND(i)+(win_ons*(hz/bin)):RevEND(i)+((win_ons+winl)*hz/bin))=1;
                elseif revl(i)>minrl*(hz/bin)
                    Runs(RevEND(i)+(win_ons*(hz/bin)):RevON(i+1))=1;
                end
            elseif revl(i)>minrl*(hz/bin) & length(ratio) > RevEND(i)+((win_ons+winl)*hz/bin)
                
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
    xl{CC}=[num2str(win_ons) '-' num2str(win_ons+winl) 's'] ;
    subratio=dRdTAll;
    subdC=dCdTAll;
    gi=find(RunsAll==1); %  timepoints within run window only
    subratio=dRdTAll(gi);
    subdC=dCdTAll(gi);
    
    binsize=0.0095;
    minDC=0.001;
    
    bin_idx= find(subdC>minDC & subdC<=i+binsize);
    
    if length(bin_idx)>length(find(~isnan(subratio)))/200
        
        bV=ratioAll(bin_idx);
        sV=sensInAll(bin_idx);
        dCV=subdC(bin_idx);
        dRV=subratio(bin_idx);
        
        rd(cc)=length(find(bV(~isnan(bV))));
        
        mBAG(1:length(bV),CC)=(bV);
        mSens(1:length(sV),CC)=(sV);
        mDC(1:length(dCV),CC)=(dCV);
        mDR(1:length(dRV),CC)=(dRV);
        
    end
end
%% %%%% plot %%%%%%%
% % boxplot/scatter:
figure
subtightplot(1,2,1)
boxplot(mDR(:,1:CC)*(hz/bin))
h=findobj(gca,'tag','Outliers');
delete(h)
set(gca,'XTick',[1:length(xl)])
set(gca,'XTickLabel',xl);
hold on
plot([0 CC],[0 0],'--k')
ylim([-0.2 0.2])
xlabel('window into run (sec)')
ylabel('dR/sec')
title([dirname(cd) ': '  num2str(hz/bin) 'hz, min rev l=' num2str(minrl) 's'])

subtightplot(1,2,2)
shadedErrorBar2([0:wstep:wstep*5 ],nanmedian(mDR(:,1:CC)),nanstd(mDR(:,1:CC))./sqrt(length(bin_idx)),'b',0);


ylim([-0.05 0.05])
% yticks([-1.2:0.3:1.1])
xlabel('window into run (sec)')
ylabel('dR/sec')
hold on
plot([0 wstep*6 ],[0 0],'--k')

%stats: one sample ttest versus mean 0
% [H, P]=ttest(mDR);
% for tt=1:CC
% text(tt,0.25,['p~=0:' num2str(roundd(P(tt),4))])
% end

title([dirname(cd) ': '  num2str(hz/bin) 'hz, min rev l=' num2str(minrl) 's'])
%     %increase window start point:
%     win_ons=win_ons+wstep;


%% save
pdir=plotdir;
cd(pdir(1).name)
saveas(gca, ['coding_diff_zero_pos' num2str(winl) 'run_wins_minrl_' num2str(minrl) '.fig'])
cd ..\


