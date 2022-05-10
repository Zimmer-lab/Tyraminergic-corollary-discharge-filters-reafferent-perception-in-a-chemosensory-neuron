%(c) Julia Riedl 2017
% finds 2nd largest thresholded region in cherry channel

function [thiswormdata2,imgnum] =getregionvalues_jr1_region2(ch_stk,gc_stk,foldername,...
    currentthreshold2,numframes,stack,thiswormdata2)

close all
warning off

bg_pos='left';

if length(thiswormdata2)<2
    thiswormdata = NaN(numframes,4);
    firstframe=1;
else
    [ns ne]=nanedge(thiswormdata2(:,1));
    if ~isempty(ne)
    if ne(end)==length(thiswormdata2)
        firstframe=ns(end);
    else
        firstframe=length(thiswormdata2);
    end
    else
        firstframe=length(thiswomrdata2)-1;
    end
end
%visualize:
scr=get(0,'ScreenSize');
figure('Name','Regions','Position', [scr(3)/1.4 scr(4)/1.5 scr(3)/4 scr(4)/4]);
loadtime=NaN(1,numframes);
tic

proofpics={};

%offset correction:
%get offset values:
corr_offset=channel_offset(ch_stk,gc_stk);
X_off=round(nanmedian(corr_offset(:,1)));
Y_off=round(nanmedian(corr_offset(:,2)));
disp(['Xoff:' num2str(X_off)])
disp(['Yoff:' num2str(Y_off)])
CC=1;

for imgnum = firstframe:numframes
    
    %------- % (1) mCHERRY: -----
    try
        tStart=tic;
        currentframe = tiffread2(ch_stk,imgnum,imgnum);
        smoothframe=smoothn(double(currentframe.data),3);
        loadtime(imgnum)=toc(tStart);
    catch
        
        save([foldername '_' num2str(stack)  '_R2_log'],'thiswormdata2');
        break
    end
    
    %% (a)threshold
    threshedframe = smoothframe>currentthreshold2;
    connectedparts = bwconncomp(threshedframe);
    for n = 1:size(connectedparts.PixelIdxList,2)
        partsize = size(connectedparts.PixelIdxList{n},1);
        if partsize < 50 % minimum 50 pixels for a region
            for m = 1:partsize
                badpixels = connectedparts.PixelIdxList{n};
                threshedframe(badpixels(m,1)) = 0;
            end
        end
    end
    connectedpartsleft = bwconncomp(threshedframe);
    
    if connectedpartsleft.NumObjects > 0
        %the following is to only select the largest region
        partsizes = NaN(1,connectedpartsleft.NumObjects);
        for n = 1:connectedpartsleft.NumObjects
            partsizes(1,n) = size(connectedpartsleft.PixelIdxList{n},1);
        end
        
        % find 2 largest objects:
        [~, idx]=sort(partsizes,'descend');
        %largest region:
        R1=zeros(size(threshedframe));
        R1(connectedpartsleft.PixelIdxList{idx(1)})=1;
        
        
        %%   (b) find second biggest region using same threshold:
        
        if connectedpartsleft.NumObjects >1
                       
            %if there's more than 1 additional object, remove very
            %elongated objects
            if connectedpartsleft.NumObjects >2
            ecc=regionprops(threshedframe, 'eccentricity','area','PixelIdxList');
            ev=[ecc.Eccentricity];
            bi=find(ev>0.98);
            for i=1:length(bi)
            badpixels = ecc(bi).PixelIdxList;
            threshedframe(badpixels) = 0;
            disp(imgnum)
            disp('too elongated')
            end
            end
            
            %2nd largest region:
            connectedpartsleft = bwconncomp(threshedframe);
            R2=zeros(size(threshedframe));
            R2(connectedpartsleft.PixelIdxList{idx(2)})=1;
            
            
            %check if region is touching the border or way too large:
            if size(threshedframe,2)==size(threshedframe,1)
                ims=size(threshedframe,1);
            end
            %points along the border:
            bp=[1:ims*2,(ims+1):ims:ims^2,(ims+2):ims:ims^2,(ims-1):ims:ims^2,(ims):ims:ims^2,numel(threshedframe)-ims*2:ims^2];
            rc=1;
            if numel(connectedpartsleft.PixelIdxList{idx(rc+1)})>10000
                rc=rc+1;
            end
            
            if length(idx)>2 % increase threshold when 2nd area is too large 
                %and there is a third object( e.g. gut flourescence)
                if numel(connectedpartsleft.PixelIdxList{idx(rc+1)})>10000
                    rc=rc+1;
                    currentthreshold2=currentthreshold2*0.95;
                else
                    R2=zeros(size(threshedframe));
                    R2(connectedpartsleft.PixelIdxList{idx(rc+1)})=1;
                    
                    while ~isempty(intersect(connectedpartsleft.PixelIdxList{idx(rc+1)},bp)) & rc<4
                        
                        if length(idx)<3
                            thiswormdata2(imgnum,1:3)=NaN(1,3);
                            currentthreshold2=currentthreshold2*0.95;
                            continue
                        end
                        try
                        R2=zeros(size(threshedframe));
                        R2(connectedpartsleft.PixelIdxList{idx(rc+2)})=1;
                        rc=rc+1;
                        catch
                            thiswormdata2(imgnum,1:3)=NaN(1,3);
                            break
                        end
                    end
                end
            else
                thiswormdata2(imgnum,1:3)=NaN(1,3);
            end
            
            %background:  takes a region of similar size as neuron, 8 neuron diameters
            %up                     
            bbR2=regionprops(R2,'Boundingbox');           
            if ~isempty(bbR2)
                bbR2=round(bbR2.BoundingBox);
                %cherry backgraound:
                if bbR2(2)+bbR2(4)*8<size(smoothframe,1)
                bckg_cherry=nanmean(nanmedian(smoothframe((bbR2(2)+bbR2(4)*7):(bbR2(2)+bbR2(4)*8),bbR2(1):bbR2(1)+15)));                               
                if imgnum==1
                    cf=(smoothframe);
                    cf((bbR2(2)+bbR2(4)*7):(bbR2(2)+bbR2(4)*8),bbR2(1):bbR2(1)+15)=15000;
                    figure
                    imagesc(cf)
                    title(dirname(cd))
                    saveas(gcf, 'background_area.fig')                  
                end
                else
                    bckg_cherry=50;
                end
                else
                thiswormdata2(imgnum,1:3)=NaN(1,3);
                continue
            end
            
            
        else % if no second region continue to next frame:
            thiswormdata2(imgnum,1:3)=NaN(1,3);
            
            continue
        end
        
        %now get values for region in cherry and gcamp channel       
        numpixelsthreshed = sum(sum(R2));
        thresholdpixvalskeptmc = smoothframe.*R2; %thresholdpixvalskept retains the pixel brightness values but loses any below threshold; this is what is saved
        %cherry intensity sum of largest area & area size
        thiswormdata2(imgnum,2) = (sum(sum(thresholdpixvalskeptmc))) -(bckg_cherry*numpixelsthreshed);
        thiswormdata2(imgnum,3) = numpixelsthreshed;
        %cherry background 1 & 2:
        thiswormdata2(imgnum,5) = bckg_cherry*numpixelsthreshed;
        
        % ---GCAMP:---
        currentframe1 = tiffread2(gc_stk,imgnum);
        currentframe1=double(currentframe1.data);
        %---correct offset---
        if X_off<0 & Y_off>=0
            currentframe2=currentframe1(1:end-Y_off,1-X_off:end);
            currentframe2=[currentframe2 zeros(size(currentframe2,1),abs(X_off))];
            currentframe2=[zeros(abs(Y_off),size(currentframe2,2));currentframe2 ];
        elseif X_off>0 & Y_off<=0
            currentframe2=currentframe1(1+abs(Y_off):end,1:end-X_off);
            currentframe2=[ zeros(size(currentframe2,1),abs(X_off)) currentframe2];
            currentframe2=[currentframe2;zeros(abs(Y_off),size(currentframe2,2)) ];
        elseif X_off>0 & Y_off>=0
            currentframe2=currentframe1(1:end-Y_off,1:end-X_off);
            currentframe2=[ zeros(size(currentframe2,1),abs(X_off)) currentframe2];
            currentframe2=[zeros(abs(Y_off),size(currentframe2,2));currentframe2 ];
        elseif X_off<0 & Y_off<=0
            currentframe2=currentframe1(1+abs(Y_off):end,1-X_off:end);
            currentframe2=[currentframe2 zeros(size(currentframe2,1),abs(X_off))];
            currentframe2=[currentframe2;zeros(abs(Y_off),size(currentframe2,2)) ];
        else
            currentframe2=currentframe1;
        end
 
        %gcamp background:
                if bbR2(2)+bbR2(4)*8<size(currentframe2,1)
                bckg_gcamp=nanmean(nanmedian(currentframe2((bbR2(2)+bbR2(4)*7):(bbR2(2)+bbR2(4)*8),bbR2(1):bbR2(1)+15)));                               
                else
                    bckg_gcamp=10;
                end
        
        thresholdpixvalskeptgc = currentframe2.*R2;
        thiswormdata2(imgnum,1) = (sum(sum(thresholdpixvalskeptgc)));
        %backgrounds:
        thiswormdata2(imgnum,3) =bckg_cherry;
        thiswormdata2(imgnum,4)=bckg_gcamp;
        
        
    end
    
    
    if mod(imgnum,100)==0
        display(['frame#' num2str(imgnum)])
        % show 2 regions
        
        ipic=zeros(size(threshedframe));
        ipic(find(R1==1))=1;
        ipic(find(R2==1))=3;
        subtightplot(1,2,1)
        imagesc(thresholdpixvalskeptgc);
        subtightplot(1,2,2)
        imagesc(ipic);
        title(num2str(imgnum))
        pause(0.01)
        assignin('base', 'thiswormdata2', thiswormdata2);
        proofpics{CC}=ipic;
        CC=CC+1;
        
    end
    
    if mod(imgnum,1000)==0
        display(num2str(imgnum))
        save([foldername '_' num2str(stack)  '_R2_log'],'thiswormdata2');
        save proofpics proofpics;
    end
    
end

%save final:
save([foldername '_' num2str(stack)  '_R2_log'],'thiswormdata2');


end