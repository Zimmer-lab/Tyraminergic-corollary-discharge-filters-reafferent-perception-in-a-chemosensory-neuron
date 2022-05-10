%(c) Julia Riedl 2021
%quantifies turning and speed from impr gradient behavior
clear
warning off
% global rois
% cd('mat')
%quantifies weathervaning (average turning bias from bearing values 45-135)
clear

xlabels={NaN};

genot={'lite1*','lgc*','WT*','BAG*'};
speed=NaN(20,length(genot)*1);
turn=NaN(20,length(genot)*1);

home=cd;
try
    
    CC=1;
    for gt=1:length(genot)
        gt_directory=dir(genot{gt});
        disp((genot{gt}))
        clearvars -except CC home gt_directory genot gt xlabels speed turn binsize
        
        for dd=1:length(gt_directory)
            
            cd(gt_directory(dd,1).name)
            
            sd=dir(genot{gt});
           qf= dir('beh_quant*');
            load(qf(1).name)
            turn(1:size(quant.turning45_135,1),CC)=quant.turning45_135;
            speed(1:size(quant.speed_foldchange,1),CC)=quant.speed_foldchange;
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
            ME1.stack
    end
    
 speed=speed(:,1:CC-1);
 turn=turn(:,1:CC-1);
%     data(data==0)=NaN;
 save('beh_quant_multi','speed','turn','xlabels')
    

%% plot:
figure
subtightplot(1,2,1)
boxplot(speed,'labels',xlabels)
title(['speed modulation:' dirname(cd)])
ylabel ('fold change')
subtightplot(1,2,2)
boxplot(turn,'labels',xlabels)
title(['turning bias:' dirname(cd)])
ylabel('turning(rad/mm)')
saveas(gca, ['beh_quant_multi.fig'])


