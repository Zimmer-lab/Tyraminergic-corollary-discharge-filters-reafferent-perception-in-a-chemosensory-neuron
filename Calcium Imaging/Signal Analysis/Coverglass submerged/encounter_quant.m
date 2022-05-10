%encounter position-triggered R/R0 vs.distance
% input: hand annotated encounter timepoints as cell array.
% e.g. for 1 experiment:
%encounter_timepoints =
%   5×1 cell array
% 
%     {'W1_1_R1_'  }    
%     {[    1.5280]}    
%     {[    1.4020]}    
%     {[2.7504e+03]}    
%     {[2.5236e+03]}   

close
clear

load 'ImagingDataAll'
load 'encounter_timepoints'
offset=0; % offset between video and visiview in sec

ne=size(encounter_timepoints,1);
hz=30; %acquisition rate
prewin=15; %in sec
prewin=ceil(prewin*hz);
enm=NaN(length(ImagingData)*3,prewin+11);
ed=NaN(length(ImagingData)*3,prewin+11);

plots=0;
if plots==1
    f1=figure;
end

cc=1;

for F=1:length(ImagingData)
    
    sr=sort(ImagingData{F}.ratio);
    R0=nanmean(sr(1:length(ImagingData{F}.ratio)/10));
    rr0=ImagingData{F}.ratio./R0;
    clear tp

    tp=encounter_timepoints{3,F}; %30hz
    tp=tp(tp~=0);
    tp=round(tp-(offset*hz));
    
    
    for jj=1:length(tp)
        if isfield(ImagingData{F},'XY') %calculate distance for event
            X=ImagingData{F}.XY(:,1);
            Y=ImagingData{F}.XY(:,2);
            
            if plots==1
                figure(f1)
                tv=[1:length(rr0)]/(30*60);
                subtightplot(1,2,1)
                if jj==1
                    cla
                    plot(tv,rr0,'b')
                    hold on
                    scatter(tv(tp),rr0(tp),'r')
                    xlabel('min')
                    title(F)
                end
                subtightplot(1,2,2)
                cla
                if tp(jj)>900
                    scatter(X(tp(jj)-900:tp(jj)+1),Y(tp(jj)-900:tp(jj)+1),10,rr0(tp(jj)-900:tp(jj)+1),'filled')
                    title(F)
                else
                    scatter(X(1:tp(jj)+10),Y(1:tp(jj)+10),10,rr0(1:tp(jj)+10),'filled')
                end
                hold on
                scatter(X(tp(jj)),Y(tp(jj)),'pk')
                axis equal
            end
            
            %distance:
            cj=1;
            if tp(jj)>prewin
                for i=tp(jj)-prewin:tp(jj)+10
                    ed(cc,cj)=edist([X(tp(jj)),Y(tp(jj))],[X(i),Y(i)]);
                    cj=cj+1;
                end
            else
                for i=1:tp(jj)+10
                    ed(cc,cj)=edist([X(tp(jj)),Y(tp(jj))],[X(i),Y(i)]);
                    cj=cj+1;
                end
            end
        end
        
        
        if tp(jj)>prewin
            enm(cc,:)=rr0(tp(jj)-prewin:tp(jj)+10);
        else
            cud=rr0(1:tp(jj)+10);
            enm(cc,end-length(cud)+1:end)=cud;
        end
        cc=cc+1;
    end
    
end

%% normalize:

enm_n=NaN(size(enm));
for i=1:size(enm,1)
    if ~isnan(nanmedian(enm(i,:)))
        sre=sort(enm(i,:));
        R0=nanmean(sre(1:length(find(~isnan(enm(i,:))))/10));
        enm_n(i,:)=enm(i,:)./R0;
        
    end
end

%% bin data:
edbin=NaN(size(ed,1),20);
erbin=NaN(size(ed,1),20);

binwin=18;

for F=1:size(ed,1)
    ccb=1;
    for i=1:binwin:790
        
        binidx=find(ed(F,:)>i & ed(F,:)<i+binwin);
        if ~isempty(binidx)
            edbin(F,ccb)=nanmedian(ed(F,binidx));
            erbin(F,ccb)=nanmedian(enm_n(F,binidx));
        else
            edbin(F,ccb)=NaN;
            erbin(F,ccb)=NaN;
        end
        ccb=ccb+1;
        
    end
    
end

%% plot

figure
subtightplot(1,2,1)
sem=nansem(enm);
shadedErrorBar2((1:length(enm))/hz,nanmean(enm),sem,[0 0.5 0.5],0);
xlabel('sec')
ylabel('R/R0')

subtightplot(1,2,2)
sem=nansem(erbin);
shadedErrorBar2(nanmean(edbin)/1000,nanmean(erbin),sem,[0 0.5 0.5],0);
xlabel('sec')
set ( gca, 'xdir', 'reverse' )
xlabel('distance in mm')
ylabel('R/R0')

suptitle(['worm-worm encounters, n=' num2str(cc-1)] )

%quantification:
rqt=nanmean(enm(:,end-30:end),2);
rqd=nanmean(erbin(:,1:3),2);

%%save
plotdir
cd(pdir(1).name)
saveas(gcf, 'encounter_quant.fig')
cd ..\




