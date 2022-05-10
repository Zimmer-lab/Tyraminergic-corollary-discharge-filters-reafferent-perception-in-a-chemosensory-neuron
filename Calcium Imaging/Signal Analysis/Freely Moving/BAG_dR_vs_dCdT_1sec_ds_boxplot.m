%((c) Julia Riedl 2018
%correlates ratio and CO2 derivatives in bins and returns a boxplot plot.
%set exR variable to exclude: exR=0 : nothing, exR=1 reversals, exR=2 forward
% downsamples data to 1hz

%Input: ImagingData All file created with data gathering code (CaData...)

clearvars -except ImagingData
%---set to exclude reversals or forward episodes:
exR=1;
%---
hz=30; % data rate
CC=1;
warning off

bagAll=[];
sensInAll=[];
dCdTAll=[];
dRdTAll=[];

plots=0;

for F=1:length(ImagingData)
    
    %ratio:
    
    cherry=medfilt1(ImagingData{F}.cherry,15);
    gcamp=medfilt1(ImagingData{F}.gcamp,15);
    nanIdx=find(isnan(cherry));
    eq_ratio=(gcamp/nanmean(gcamp))./(cherry/nanmean(cherry));
    
    %ratio=smoothn((gcamp./cherry),10);
    ratio=smoothn(eq_ratio,20);
    ratio(nanIdx)=NaN;
    %bin for 1s:
    ratio_ds=bindata(ratio,hz);
    nanIdx=find(isnan(ratio_ds));
    ratio_ds=smoothn(ratio_ds,2);
    ratio_ds(nanIdx)=NaN;
    
    %CO2
    sensIn=ImagingData{F}.CO2;
    sensIn=smoothn(sensIn,10);
    sensIn(nanIdx)=NaN;
    
    %bin for 1s:
    sensIn_ds=bindata(sensIn,hz);
    try
        ratio_ds=ratio_ds(1:length(sensIn_ds));
    catch
        sensIn_ds=sensIn_ds(1:length(ratio_ds));
    end
    
    %reversals:
    if ~isfield(ImagingData{F},'RevFrames30hz')
        F;
        disp('no reversals annotated!')
        continue
    end
    RevON=ImagingData{F}.RevFrames30hz(1:2:end);
    RevEND=ImagingData{F}.RevFrames30hz(2:2:end);
    Revs=zeros(length(ratio),1);
    for i=1:length(RevEND)
        if RevON>2
            Revs(RevON(i)-2:RevEND(i)+2,1)=1;
        else
            Revs(RevON(i):RevEND(i)+2,1)=1;
        end
        
    end
    Revs=Revs(1:length(ratio));
    %     revAll=[revAll;Revs];
    
    %----set  reversal or forward  episodes to NaN: ----
    if exR==1
        ri=find(Revs==1);
    elseif exR==2
        ri=find(Revs==0);
    end
    
    %downsample:
    ri=round(ri/hz);
    ri(ri==0)=[];
    
    if exR>0
        %remove reversal or forward episodes:
        ratio_ds(ri)=NaN;
        sensIn_ds(ri)=NaN;
    end
    
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

    bag=ratio_ds';
    
    try
        sensIn=sensIn(1:length(bag));
        
    catch
        sensIn_ds=[sensIn_ds NaN(1,length(bag)-length(sensIn))];
        sensIn_ds=sensIn_ds(1:length(bag));
    end
    
    bagAll=[bagAll; bag];
    
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
    
    if length(dRdTAll)~=length(dCdTAll)
        disp('inequal vector length')
        break
    end
end

%% %%%%%% bin & plot %%%%%%%

% % bin for dCdt

mBAG=NaN(10000,16);
mSens=NaN(10000,16);
mDC=NaN(10000,16);
mDR=NaN(10000,16);
binsize=0.0025;

cc=1;
for i=-0.02 :binsize:(0.018-binsize)
    
    bin_idx= find(dCdTAll>i & dCdTAll<=i+binsize);
    
    if length(bin_idx)>length(find(~isnan(dCdTAll)))/500
        
        bV=bagAll(bin_idx);
        sV=sensInAll(bin_idx);
        dCV=dCdTAll(bin_idx);
        dRV=dRdTAll(bin_idx);
        
        rd(cc)=length(find(bV(~isnan(bV))));
        
        mBAG(1:length(bV),cc)=(bV);
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
set(gca,'XTickLabel',round(xl(1:2:end)*1000)/1000);
ylim([-0.10 0.10])
xlabel(' CO2 change (%/sec)')
ylabel(' dR/dt')
hold on
plot([0 21],[0 0],'--k')
if exR==0
title([dirname(cd) ': all data'])
elseif exR==1
title([dirname(cd) ': reversals only'])
elseif exR==2
title([dirname(cd) ': forward only'])
end





