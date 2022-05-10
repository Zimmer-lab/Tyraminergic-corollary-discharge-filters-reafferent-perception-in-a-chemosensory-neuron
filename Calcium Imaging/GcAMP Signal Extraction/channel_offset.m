% calculates position offset between gcamp and cherry channel
% find cherry & gc stk:
function corr_offset=channel_offset(ch_stk,gc_stk)
warning off
stks= dir(['*.stk']);
stknames={stks.name};

if ch_stk==0 & gc_stk==0
    ch_i=find( cellfun('isempty', strfind( stknames,'gc')));
    % SN=input('which stk?');
    ch_stk=stks(ch_i(1)+0).name;
    
    if length(stks)==2
        gc_stk=stks(ch_i(1)+1).name;
    else
        gc_stk=stks(ch_i(1)+2).name;
    end
end

c=1;
corr_offset=NaN(1,2);
for ff= 1:500:7480
    try
        gc= tiffread2(gc_stk,ff);
        ch = tiffread2(ch_stk,ff);
        [max_ch, imax] = max(abs(ch.data(:)));
        [y, x] = ind2sub(size(ch.data),imax(1));
        sub_ch=NaN;
        try
            for margin=40:20:80
            sub_ch=medfilt1(ch.data(y-margin:y+margin,x-margin:x+margin));
            sub_gc=medfilt1(gc.data(y-margin:y+margin,x-margin:x+margin));
            end
        catch
            if isnan(sub_ch)
            disp('no valid max')
            continue
            end
        end
        cr = normxcorr2(sub_gc(:,:,1),sub_ch(:,:,1));
        [max_c, imax] = max(abs(cr(:)));
        [ypeak, xpeak] = ind2sub(size(cr),imax(1));
        corr_offset(c,:) = [(xpeak-size(sub_gc,2)) (ypeak-size(sub_gc,1))];
        c=c+1;
    catch ME
    end
end
    
disp(corr_offset)