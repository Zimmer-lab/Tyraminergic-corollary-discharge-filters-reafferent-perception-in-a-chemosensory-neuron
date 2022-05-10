
clear
warning off
%parameters:
minB=1; %in sec, bin lag after exp start for quantif

xlabels={NaN};

genot={'lite*','WT*','lgc*','BAG_r*','AVB*','BAG_A*'};
quant_var='rq.mat';

data=NaN(70,ceil(length(genot)*1.3));
data_length=NaN;
home=cd;

try
    CC=1;
    for gt=1:length(genot)
        gt_directory=dir(genot{gt});
        disp((genot{gt}))
        
        for dd=1:length(gt_directory)
            
            cd(gt_directory(dd,1).name)
            
            sd=dir(genot{gt});
            load(quant_var)
            data(1:size(rq,1),CC)=rq;
            data_length(CC)=size(rq,1);
            xlabels{CC}=dirname(cd);
            CC=CC+1;
            cd ..\
        end % end subdirectories
    end
    
    cd(home)
    
catch ME1
    beep
    cd(home)
    disp(ME1.message)
end

data=data(1:max(data_length),1:CC-1);
save('revModulation_quant','data','xlabels')


%% plot:
figure
boxplot(data,'labels',xlabels)
title(['Reversal Frequency Modulation:' dirname(cd)])
ylabel('fold change rev freq')
xtickangle(25)
saveas(gca, ['reversal_modulation_multi.fig'])


