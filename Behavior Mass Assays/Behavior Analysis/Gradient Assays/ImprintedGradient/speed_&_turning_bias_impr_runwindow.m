%% weathervaning hist
%plot?
clear
plots=1;
files =dir('runinfo_min3*');
maxt=27;% maximum time into run in 3hz
mint=9;% minimum time into run in 3hz
hz=3;

% if isempty(strfind(files(1).name,'rev'))
load(files(1).name);
disp(files(1).name);
% else
%     load(files(2).name);
%     disp(files(2).name);
% end
% for each experiment: put all run data into one vector:
mB= NaN(length(bearing),18);
mT= NaN(length(bearing),18);
mS= NaN(length(bearing),18);
mdS= NaN(length(bearing),18);
c=0;

for F=1:length(bearing)
    
    disp(F)
    if ~isempty(bearing{F})
        
        c=c+1;
        bearingAll=[];
        dBearingAll=[];
        SpeedAll=[];
        diffSpeedAll=[];
        DAll=[];
        
        
        for i=1:size(bearing{F},1)

            if ~isempty(dBearing{F}{i,1})
                for rr=1:size(bearing{F}(i,:),2)
                    if  ~isempty(speed_ctr{F}{i,rr})
                        bv=bearing{F}{i,rr};
                        tr=dBearing{F}{i,rr};
                        sp=smoothn(speed_ctr{F}{i,rr},10);
                        di=DistanceTrav{F}{i,rr};
                        fr=RunFrames{F}{i,rr};
                        if length(bv)>maxt %& fr(1)<mexp*3*60
                            bv1=bv(mint:maxt);
                            tr1=tr(mint:maxt);
                            dsp1=[diff(sp(mint:maxt)) NaN];
                            sp1=sp(mint:maxt);
                            di1=di(mint:maxt);
                        elseif length(bv)>mint & length(bv)<maxt %& fr(1)<mexp*3*60
                            bv1=bv(mint:end);
                            tr1=tr(mint:end);
                            dsp1=[diff(sp(mint:end)) NaN];
                            sp1=sp(mint:end);
                            di1=di(mint:end);
                        else
                            bv1=NaN;
                            tr1=NaN;
                            dsp1= NaN;
                            sp1=NaN;
                            di1=NaN;
                            
                        end
                    else
                        bv1=NaN;
                        tr1=NaN;
                        dsp1= NaN;
                        sp1=NaN;
                        di1=NaN;
                    end
                    
                    bearingAll=cat(2,bearingAll,bv1);
                    dBearingAll=cat(2,dBearingAll,tr1);
                    SpeedAll=cat(2,SpeedAll,sp1);
                    diffSpeedAll=cat(2,diffSpeedAll,dsp1);
                    DAll=cat(2,DAll,di1);
                    
                end
            end
        end
        
        
        bv=(bearingAll);
        tr=(dBearingAll);
        sp=(SpeedAll);
        spd=(diffSpeedAll);
        
        %kill extremely high turning rates (short reversals etc)
        
        bi=tr>19 | tr<-19;
        tr(bi)=NaN;
        bv(bi)=NaN;
        sp(bi)=NaN;
        spd(bi)=NaN;
        tnorm=tr./DAll;
        bi=tnorm>15 | tnorm<-15;
        tnorm(bi)=NaN;
        tnorm=deg2rad(tnorm);  %conversion from deg/pix-->rad/px 1/(67*(360/(2*pi)))
        tnorm=(tnorm)/0.025;  %conversion from rad/pix-->deg/mm
        
        %      figure(fig2)
        %      try
        %      subplot(1,5,F)
        %      hist(bv,20)
        %      end
        
        if ~isempty(bv)
            
            [X,v]=hist(bv,10);
            X1(c,:)=X./sum(X);
            cc=1;
            %bin -->for what?
            a=15;
            
            for i=10 :a:165
                bin_idx= bv>i & bv<=i+a;
                bin_idx=find(bin_idx);
                pb1=bv(bin_idx);
                pt=tnorm(bin_idx);
                st=sp(bin_idx);
                dst=spd(bin_idx);
                
                mB(c,cc)=nanmean(pb1);
                mT(c,cc)=nanmean(round(pt*10)/10);
                mS(c,cc)=nanmean(st);
                mdS(c,cc)=nanmean(dst);
                
                cc=cc+1;
            end
            
        end
    end
end
save mT mT

%% quantification:
%1: turning bias
quant.turning45_135=nanmean(mT(:,4:10),2);
%peed modulation
quant.speed_foldchange=nanmean(mS(:,1:3),2)./nanmean(mS(:,cc-3:cc-1),2);
save(['beh_quant_' num2str(round(mint/hz))], 'quant');
%% plot:(1) turning rate binned by bearing
if plots==1
    
    figure
    subtightplot(1,3,1) %weathervaning
    sem=nanstd(mT(:,1:cc-1),1)/sqrt(c);
    hold on
    bar(nanmean(mT(:,1:cc-1),1));
    errorb(nanmean(mT(:,1:cc-1),1),sem)
    set(gca,'XTick',[1:cc-1])
    set(gca,'XTickLabel',round(nanmedian(mB,1)));
    xlabel('bearing')
    ylabel('turning rate (rad/mm)')
    ylim([-0.14 0.05])
    
    subtightplot(1,3,2) %speed
    sem=nanstd(mS(:,1:cc-1),1)/sqrt(c);
    hold on
    plot(nanmedian(mB(:,1:cc-1),1),nanmean(mS(:,1:cc-1),1));
    errorb(nanmedian(mB(:,1:cc-1),1),nanmean(mS(:,1:cc-1),1),sem)
    % set(gca,'XTick',[1:cc-1])
    % set(gca,'XTickLabel',round(nanmedian(mB,1)));
    xlabel('bearing')
    ylabel('speed(mm/sec)')
    ylim([0.02 0.1])
    
    subtightplot(1,3,3) %acceleration
    sem=nanstd(mdS(:,1:cc-1),1)*hz/sqrt(c);
    hold on
    plot(round(nanmedian(mB(:,1:cc-1),1)),nanmean(mdS(:,1:cc-1),1)*hz);
    errorb(round(nanmedian(mB(:,1:cc-1),1)),nanmean(mdS(:,1:cc-1),1)*hz,sem)
    % set(gca,'XTick',[1:cc-1])
    % set(gca,'XTickLabel',round(nanmedian(mB,1)));
    xlabel('bearing')
    ylabel('acceleration(mm/sec^2)')
    % xlim([1:cc-1])
    ylim([-0.005 0.005])
    
end % end files loop

%quantification
TB(:,1)=nanmean(mT);

nd=(cd);
d= strfind(cd, '\');
name=nd(d(end)+1:end);
suptitle([name(1:end) ':' num2str(round(mint/3.3)) 's-' num2str(round(maxt/3.3)) ' into run; speed_bias'])

%% save:
plotdir=dir('*plots*');
if isempty(plotdir)
    mkdir([name ' plots'])
    plotdir=dir('*plots*');
end
cd(plotdir(1).name)


saveas(gca, ['speed_bias_maxrunl' num2str(round(maxt/3.3)) '.fig'])
cd ..\




