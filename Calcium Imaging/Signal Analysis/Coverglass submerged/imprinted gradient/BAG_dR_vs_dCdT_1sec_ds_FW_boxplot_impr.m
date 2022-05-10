%ratio vs CO2 derivative bin plot

clear
load ImagingDataAll.mat

CC=1;

warning off

bagAll=[];
afdAll=[];
sensInAll=[];
dCdTAll=[];
dRdTAll=[];

binning=30;

plots=0;

for F=1:length(ImagingData)
    
    %if isfield(ImagingData{F}, 'cherry_AFD')==1
    
    %ratio:
    
    cherry=medfilt1(ImagingData{F}.cherry_bag,15);
    gcamp=medfilt1(ImagingData{F}.gcamp_bag,15);
    nanIdx=find(isnan(cherry));
    eq_ratio=(gcamp/nanmean(gcamp))./(cherry/nanmean(cherry));
    
    %ratio=smoothn((gcamp./cherry),10);
    ratio=smoothn(eq_ratio,20);
    ratio(nanIdx)=NaN;
    %bin for 1s:
    ratio_ds =NaN;
    bc=1;
    for i=1:binning:length(ratio)-binning
        ratio_ds(bc)=nanmedian(ratio(i:i+binning));
        bc=bc+1;
    end
    nanIdx=find(isnan(ratio_ds));
    ratio_ds=smoothn(ratio_ds,2);
    ratio_ds(nanIdx)=NaN;
    
    %CO2
    sensIn=ImagingData{F}.CO2;
    sensIn=smoothn(sensIn,binning);
    sensIn(nanIdx)=NaN;
    
    %bin for 1s:
    sensIn_ds =NaN;
    bc=1;
    for i=1:binning:length(sensIn)-binning
        sensIn_ds(bc)=nanmedian(sensIn(i:i+binning));
        bc=bc+1;
    end
    sensIn=smoothn(sensIn,2);
    sensIn(nanIdx)=NaN;
    try
        ratio_ds=ratio_ds(1:length(sensIn_ds));
    catch
        sensIn_ds=sensIn_ds(1:length(ratio_ds));
    end
    
    %reversals:
    if ~isfield(ImagingData{F},'RevFrames30hz')
        F
        disp('no reversals annotated!')
        continue
    end
    RevON=ImagingData{F}.RevFrames30hz(1:2:end);
    
    RevEND=ImagingData{F}.RevFrames30hz(2:2:end);
    Revs=zeros(length(ratio),1);
    RevON(RevON<1)=1;
    if ~isnan(RevON)
    for i=1:length(RevEND)
        if RevON(i)>2
            Revs(RevON(i)-2:RevEND(i)+2,1)=1;
        else
            Revs(RevON(i):RevEND(i)+2,1)=1;
        end
        
    end
    end
    Revs=Revs(1:length(ratio));
    %     revAll=[revAll;Revs];
    
    %----set  forward  episodes to NaN: ----
    ri=find(Revs==1);
    %downsample:
    ri=round(ri/binning);
    ri(ri==0)=[];
    
    %remove reversal or forward episodes:
    ratio_ds(ri)=NaN;
    sensIn_ds(ri)=NaN;
    
    %     %adjust vector length:
    %     bearingGlob=bearingGlob(1:length(eq_ratio_norm));
    
    % plot:
    if plots==1
        figure(f1)
        cla
        %subtightplot(4,4,CC)
        hold on
        plot(sensIn,'k')
        plot(ratio,'b')
        %     plot(eq_ratio_norm_AFD,'color',[1 0 0.5])
        try
            h=scatter(ri,sensIn(ri),'marker','.','SizeData',1);
            set(h,'MarkerEdgeColor',[0.5 0.5 0.5]);
        end
        set(gca,'Xgrid','on');
    end
    
    %% concat data:
    %     try
    %         afd=eq_ratio_norm_AFD;
    %         bag=eq_ratio_norm(1:length(afd));
    %     catch
    %         bag=eq_ratio_norm;
    %         afd=eq_ratio_norm_AFD(1:length(bag));
    %     end
    bag=ratio_ds';
    
    try
        sensIn=sensIn(1:length(bag));
        
    catch
        sensIn_ds=[sensIn_ds NaN(1,length(bag)-length(sensIn))];
        sensIn_ds=sensIn_ds(1:length(bag));
    end
    
    bagAll=[bagAll; bag];
    %   afdAll=[afdAll; afd];
    
    %derivatives:
    sensInAll=[sensInAll sensIn_ds] ;
    dRdT=[diff(bag); NaN];
    dRdTAll=[dRdTAll ; dRdT];
    
    dCdT=[diff(sensIn_ds) NaN];
    dCdTAll=[dCdTAll dCdT];
    
    % plot derivatives:
    if plots==1
        
        figure(f2)
        cla
        %subtightplot(4,4,CC)
        hold on
        plot(dRdT,'color',[0.3 0.3 1])
        plot(dCdT,'k')
        %     plot(eq_ratio_norm_AFD,'color',[1 0 0.5])
        try
            h=scatter(ri,dCdT(ri),'marker','.','SizeData',1);
            set(h,'MarkerEdgeColor',[0.5 0.5 0.5]);
        end
        set(gca,'Xgrid','on');
        
    end
    
    np=num2str(CC);
    if CC<10
        np= [num2str(0) num2str(CC)];
    end
    
    CC=CC+1;
    
    %  end
    if length(dRdTAll)~=length(dCdTAll)
        disp('inequal vector length')
        break
    end
end

%%%%%% bin & plot %%%%%%%

%%
dCdTAll(dCdTAll== min(dCdTAll))=NaN;
dCdTAll(dCdTAll== 0)=NaN;
% % bin for dCdt

disp('derivative')
mBAG=NaN(10000,16);
mAFD=NaN(10000,16);
mSens=NaN(10000,16);
mDC=NaN(10000,16);
mDR=NaN(10000,16);
binsize=0.0022;

cc=1;
for i=-0.016 :binsize:(0.016-binsize)
    
    bin_idx= find(dCdTAll>i & dCdTAll<=i+binsize);
    
    if length(bin_idx)>length(find(~isnan(dCdTAll)))/300
        
        bV=bagAll(bin_idx);
        %         aV=afdAll(bin_idx);
        sV=sensInAll(bin_idx);
        dCV=dCdTAll(bin_idx);
        dRV=dRdTAll(bin_idx);
        
        rd(cc)=length(find(bV(~isnan(bV))));
        
        mBAG(1:length(bV),cc)=(bV);
        %         mAFD(1,cc)=nanmedian(aV);
        mSens(1:length(sV),cc)=(sV);
        mDC(1:length(dCV),cc)=(dCV);
        mDR(1:length(dRV),cc)=(dRV);
        
    end
    cc=cc+1;
end


%% %%%% plot %%%%%%%

% boxplot/scatter:

figure
boxplot(mDR,'symbol','.w','whisker',0.72)
h=findobj(gca,'tag','Outliers');
delete(h)
xl=nanmedian(mDC);
set(gca,'XTick',[1:2:length(xl)])
set(gca,'XTickLabel',round(xl(1:2:end)*10000)/10000);
ylim([-0.10 0.10])
xlabel(' CO2 change (%/sec)')
ylabel(' dR/dt')
hold on
plot([0 21],[0 0],'--k')
title([dirname(cd) ': fw'])

ni=~isnan(xl);
yl=nanmedian(mDR);
yl=yl(ni);
xl=xl(ni);
xl=xl(:);yl=yl(:);
[r , p]=corr(yl,xl) %correlation
text(1,0.07,['r=' num2str(r) ,'; p=' num2str(p)])
%split correlations for up and downgradient:
[rp , pp]=corr(yl(xl>=0),xl(xl>=0));
[rn , pn]=corr(yl(xl<=0),xl(xl<=0));

%% save
pdir=plotdir;
cd(pdir(1).name)
saveas(gca, ['dR_vs_dCO2_FW' dirname2(cd)])
cd ..\






