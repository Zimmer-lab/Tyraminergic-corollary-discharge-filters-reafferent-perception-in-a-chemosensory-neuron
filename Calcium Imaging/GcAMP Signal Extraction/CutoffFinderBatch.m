% code finds threshold from cherry stack for getting neuron region
%in various folders:
close
clear
thresholds=[3000 6000];

home=cd;
cc=1;
warning off
worms=dir('*W*');
start=input('start folder?');


for W=start:length(worms)
    
%     try
if ~isdir(worms(W).name)
    continue
end
%         
        cd(worms(W).name)
        disp(['W: ' num2str(W)])
        disp(worms(W).name)
        rfiles=dir('*R1*');
        if ~isempty(rfiles)
            disp('already analyzed')
            cd(home)
            continue
        end
        % find cherry stk:
        stks= dir(['*.stk']);
        stknames={stks.name};
        ch_i=find( cellfun('isempty', strfind( stknames,'_gc.')));
        % SN=input('which stk?');
        ch_stk=stks(ch_i(1)).name;
        disp(ch_stk)
        
        if exist('f1')==0
            f1=figure('Name',ch_stk,'Position' ,[730 500 540 400]);
        end
        set(f1,'Name',ch_stk);
 
        pause(0.5)
        
        for tr= 1:length(thresholds)
            
            threshold=thresholds(tr);
            disp(threshold)
            
            for F= 1:900:7250
                try
                    figure(f1)
                    currentframe = tiffread2(ch_stk,F);
                catch
                    disp(F)
                    continue
                end
                if ~isempty(currentframe)
                    smoothframe=smoothn(double(currentframe(1).data),3);
                    subplot(1,2,1)
                    imagesc(smoothframe)
                    %    caxis([50 600])
                    title([ch_stk(1:10) ': # '  num2str(F)])
                    if F==1
                        pause(1)
                    end
                    %%
                    threshedframe = smoothframe>threshold;
                    subplot(1,2,2)
%                     if F==1
%                         suptitle(ch_stk)
%                     end
                    imagesc(threshedframe)
                    title(['W' num2str(W) ':  '  num2str(threshold)])
                    
                    pause(0.5)
                end
                
            end

            pause(0.5)
        end
        
        %background:
       try
        bg=mean(mean(smoothframe(450:500,1:50)))
       catch
          bg=mean(mean(smoothframe(240:256,1:50)))
       end
        
        cd(home)
        pause(3)
%     catch
%         disp('error')
%         cd(home)
%     end
end



