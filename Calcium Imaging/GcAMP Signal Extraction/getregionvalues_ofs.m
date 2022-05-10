%(c) Julia Riedl 2017
% finds largest thresholded region and corrects for offset bewtween gcam
% and cherry channel


function [thiswormdata,thiswormdata2,imgnum] = getregionvalues_ofs(ch_stk,gc_stk,foldername,trackinstruction,...
    thistrackthreshold,currentthreshold2,thistrackbckgd,numframes,stack,thiswormdata)

pp=1; %if pp=1 proof pics of thresholded gcamp channel are stored

%offset correction:
%get offset values:
disp(datestr(now))


corr_offset=channel_offset(ch_stk,gc_stk);
correctionflag=0;

if std(corr_offset(:,1))>15 | std(corr_offset(:,2))>15
    disp ('huge variance in offset, calculate for each frame!')
    %     correctionflag=1;
    proof_pics={};
else
    correctionflag=0;
end

X_off=round(nanmedian(corr_offset(:,1)));
Y_off=round(nanmedian(corr_offset(:,2)));
if isnan(Y_off)
    Y_off=0;
end
if abs(X_off)>13 | abs(Y_off)>13
    disp('offset too large! set to 1')
    X_off=1;
    Y_off=1;
end
disp(['Xoff:' num2str(X_off)])
disp(['Yoff:' num2str(Y_off)])

close all

show=1;
thiswormdata2=NaN(1,3);

if length(thiswormdata)<2
    thiswormdata = NaN(numframes/2,4);
    firstframe=1;
else
    [ns ne]=nanedge(thiswormdata(:,1));
    if ~isempty(ne) & ne(end)==length(thiswormdata)
        firstframe=ns(end);
    else
        firstframe=length(thiswormdata);
    end
end

tic
eflag=1;
cc=1;
for imgnum = firstframe:numframes
    
    assignin('base','imgnum',imgnum)
    
    %------- % (1) mCHERRY: -----
    try
        currentframe = tiffread2(ch_stk,imgnum,imgnum);
        if ~isempty(currentframe)
            currentframe=(currentframe.data);
        else
            eflag=eflag+1;
            if eflag>20
                save ([foldername '_' num2str(stack) '_R1_log'],'thiswormdata');
                break
            else
                continue
            end
        end
        
    catch ME
        save ([foldername '_' num2str(stack) '_R1_log'],'thiswormdata');
        disp(ME.stack(1))
        continue
    end
    
    threshedframe = currentframe>thistrackthreshold-1; %have to subtract 1 from the metamorph threshold to get this to be the same
    connectedparts = bwconncomp(threshedframe);
    for n = 1:size(connectedparts.PixelIdxList,2)
        partsize = size(connectedparts.PixelIdxList{n},1);
        if partsize < 30
            for m = 1:partsize
                badpixels = connectedparts.PixelIdxList{n};
                threshedframe(badpixels(m,1)) = 0;
            end
        end
    end
    connectedpartsleft = bwconncomp(threshedframe);
    
    %check for low values of eccentrixity (round objects):
    ecc=regionprops(threshedframe,'Eccentricity','PixelIdxList');
    for jj=1:length(ecc)
        if ecc(jj).Eccentricity>0.86
            threshedframe(ecc(jj).PixelIdxList)=0;
            connectedpartsleft = bwconncomp(threshedframe);
        end
    end
    
    if connectedpartsleft.NumObjects > 0
        if trackinstruction == 0
            %the following is to only select the largest region
            partsizes = NaN(1,connectedpartsleft.NumObjects);
            for n = 1:connectedpartsleft.NumObjects
                partsizes(1,n) = size(connectedpartsleft.PixelIdxList{n},1);
            end
            maxpartindex = find(partsizes == max(partsizes));
            if size(maxpartindex,2) > 1
                maxpartindex = maxpartindex(1,1);
            end
            
            
            partsizes(maxpartindex) = 0;
            connectedpartsleft2 = connectedpartsleft;
            threshedframe2 = threshedframe;
            connectedpartsleft.PixelIdxList{maxpartindex} = [];
            for n = 1:size(connectedpartsleft.PixelIdxList,2)
                partsize = size(connectedpartsleft.PixelIdxList{n},1);
                for m = 1:partsize
                    badpixels = connectedpartsleft.PixelIdxList{n};
                    if size(connectedpartsleft.PixelIdxList{n},1)>0
                        threshedframe(badpixels(m,1)) = 0; %removes all pixels which are too small or unconnected
                    end
                end
            end
            
        end
        
        
        threshedframe = uint16(threshedframe);  %threshedframe is logical matrix (turned to uint16 to allow next line) only giving the pixel locations
        numpixelsthreshed = sum(sum(threshedframe));
        thresholdpixvalskeptmc = currentframe.*threshedframe; %thresholdpixvalskept retains the pixel brightness values but loses any below threshold; this is what is saved
        %cherry intensity sum of largest area & area size
        %dynamic background:
        bg_pix= currentframe<thistrackbckgd(1)*1.2;
        background=double(currentframe).*double(bg_pix);
        background=nanmean(nanmean(double(currentframe).*double(bg_pix)));
        if isnan(background)
            background=thistrackbckgd(1);
        end
        thiswormdata(imgnum,2) = (sum(sum(thresholdpixvalskeptmc))) -(background(1,1)*numpixelsthreshed);
        thiswormdata(imgnum,3) = numpixelsthreshed;
        
        
        
        % ---GCAMP:---
        
        currentframe1 = tiffread2(gc_stk,imgnum,imgnum);
        if isempty(currentframe1)
            thiswormdata(imgnum,1)=NaN;
            continue
        end
        currentframe1=currentframe1.data;
        if X_off==0
            X_off=1;
        end
        %---correct offset---
        if correctionflag==1
            if mod(imgnum,20)==0
                [max_ch, imax] = max(abs(currentframe(:)));
                [y, x] = ind2sub(size(currentframe),imax(1));
                sub_ch=NaN;
                sub_gc=NaN;
                try
                    for margin=40:50:90
                        sub_ch=medfilt1(currentframe(y-margin:y+margin,x-margin:x+margin));
                        sub_gc=medfilt1(currentframe1(y-margin:y+margin,x-margin:x+margin));
                    end
                catch
                    if isnan(sub_ch) | isnan(sub_gc)
                        disp('no valid max')
                        continue
                    end
                end
                cr = normxcorr2(sub_gc(:,:,1),sub_ch(:,:,1));
                [max_c, imax] = max(abs(cr(:)));
                [ypeak, xpeak] = ind2sub(size(cr),imax(1));
                corr_offset = [(xpeak-size(sub_gc,2)) (ypeak-size(sub_gc,1))];
                X_off=round(corr_offset(1));
                Y_off=round(corr_offset(2));
            end
        end
        
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
        
        
        thresholdpixvalskeptgc = currentframe2.*threshedframe;
        thiswormdata(imgnum,1) = (sum(sum(thresholdpixvalskeptgc))) - (thistrackbckgd(2,1)*numpixelsthreshed);
        if isnan(thiswormdata(imgnum,1))
            imgnum
        end
        if imgnum==1 | imgnum==firstframe
            figure
            imagesc(currentframe)
            title('cherry')
            figure
            imagesc(currentframe1)
            title('gc')
            figure
            imagesc(currentframe2)
            title('corrected')
        end
        
        if pp==1 | correctionflag==1
            if mod(imgnum,20)==0
                proof_pics{cc}=thresholdpixvalskeptgc;
                cc=cc+1;
            end
        end
        
        if mod(imgnum,100)==0
            
            display(['frame#' num2str(imgnum)])
            if show==1
                if exist('f1')==0
                    f1=figure('Name',ch_stk,'Position' ,[730 500 540 400]);
                end
                imagesc(thresholdpixvalskeptgc);
                hold on
                title(imgnum)
                pause(0.01)
                %                 disp(corr_offset)
                save ([foldername '_' num2str(stack) '_R1_log'],'thiswormdata');
                if correctionflag==1 | pp==1;
                    save  proof_pics  proof_pics
                end
            end
            
        end
        
        if mod(imgnum,1000)==0
            toc
            display(num2str(imgnum))
            save ([foldername '_' num2str(stack) '_R1_log'],'thiswormdata');
            save (['Trackparameters_' date] , 'thistrackthreshold', 'thistrackbckgd','X_off','Y_off' );
            tic
        end
    else
        thiswormdata(imgnum,:)=NaN;
        if mod(imgnum,5)==0 | mod(imgnum,10)==0
            disp([num2str(imgnum) ': no object found..'])
        end
    end
    
    %save final:
    if imgnum==numframes
        save ([foldername '_' num2str(stack) '_R1_log'],'thiswormdata');
        %save parameters
        
        save (['Trackparameters_' date] , 'thistrackthreshold', 'thistrackbckgd','X_off','Y_off' );
        disp('done')
    end
    
    
end