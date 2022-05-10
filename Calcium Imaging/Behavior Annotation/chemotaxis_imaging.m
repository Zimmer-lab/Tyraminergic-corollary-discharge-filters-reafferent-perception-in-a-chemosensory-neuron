%(c) Julia Riedl 2019
%calculates a preference PI over time based on position in CO2 gradient of various
%single worm imaging experiments
clear

try
    load ImagingDataAll
catch
    load ImagingData
end
home=cd;
cc=1;

warning off


timepoint=0:2:9.5; %in minutes
timepoint=timepoint*(60);
bins=0.5:0.2:2.6;

hist_data=NaN(length(bins),2);

plots=0;
ctl_o='123';

CO2All=[];
timeAll=[];

for F=1:length(ImagingData)
    
    if isfield(ImagingData{F}, 'RevFrames30hz')==1
        ctl=ImagingData{F}.TrialLabel;
        try
        sc=ctl(1:length(ctl_o))==ctl_o;     
        catch
            sc=ctl==ctl_o(1:length(ctl));   
        end
        
        %CO2
        if ctl(1:3)==ctl_o(1:3) & ctl(1)=='W' %concatenate data if they belong to a split trila dataset
            CO2=[ImagingData{F-1}.CO2 ImagingData{F}.CO2];
            CO2All=CO2All(1:end-length(ImagingData{F-1}.CO2)/10);
            timeAll=timeAll(1:end-length(ImagingData{F-1}.CO2)/10);
        elseif sc(1:end-4)==1 & sc(end)==0 & length(ctl)==length(ctl_o)
            CO2=[ImagingData{F-1}.CO2 ImagingData{F}.CO2];
            CO2All=CO2All(1:end-length(ImagingData{F-1}.CO2)/10);
            timeAll=timeAll(1:end-length(ImagingData{F-1}.CO2)/10);
        else
            CO2=ImagingData{F}.CO2;
        end
        
        
        nanIdx=find(isnan(CO2));
        CO2=smoothn(CO2,30);
        CO2(nanIdx)=NaN;
        
        %downsample to 3 hz:
        CO2=bindata(CO2,10);
        CO2All=[CO2All;CO2(:)];
        
        %assign time (sec)
        time=[1:length(CO2)]/3;
        timeAll=[timeAll ; time(:)];
                
        ctl_o=ctl;
             
    end %if
 
end

for T=1:length(timepoint)-1 %time points in minutes
    si=find(timeAll>=timepoint(T) & timeAll<timepoint(T+1));
    sCO2=CO2All(si);
    
    %hist:
    c=hist(sCO2,bins)/sum(~isnan(sCO2));
    hist_data(:,T)=c;
    lt{T}=[num2str(timepoint(T)/60) 'min'];
end
%%
%%%%%% plot %%%%%%%

figure
subtightplot(1,2,1)
plot(bins,hist_data)
title(dirname(cd))
xlabel('CO2')
legend(lt)

%INDEX (F1-F2)./(F2+F1)
bn=floor(length(bins)/2);
ci=sum(hist_data(bn:end,:))-sum(hist_data(1:bn,:));
subtightplot(1,2,2)
plot(ci)
xlabel('min')
ylabel('PI')
ylim([-0.6 0.6])

%save
plotdir
cd(pdir(1).name)
saveas(gca, ['imaging_PI_' dirname2(cd)])
cd ..\






