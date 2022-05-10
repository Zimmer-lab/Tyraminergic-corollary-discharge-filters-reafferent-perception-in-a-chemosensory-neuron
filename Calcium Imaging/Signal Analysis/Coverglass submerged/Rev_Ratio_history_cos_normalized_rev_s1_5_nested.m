%stacks the ratios pre and post all reversals and saves quantification in
%nested format for statistics
clear

load ImagingDataAll

%-- parameters:--
hz=30;
minlength=2;%sec
minlength=minlength*hz;

prewin=10; %in sec
prewin=prewin*hz;

postwin=10; %in sec
postwin=postwin*hz;

winsize=prewin+postwin;

revAll=NaN(20,length(ImagingData));
prerev_n_MoM=NaN(1,winsize);
perev_d_MoM=NaN(1,winsize);

for F=1:length(ImagingData)
    
    prerev=NaN(1,1800);
    cc=1;
    
    %-----ratio BAG:-------------
    try
        cherry_bag=medfilt1(ImagingData{F}.cherry_bag,5);
        gcamp_bag=medfilt1(ImagingData{F}.gcamp_bag,5);
    catch
        cherry_bag=medfilt1(ImagingData{F}.cherry,5);
        gcamp_bag=medfilt1(ImagingData{F}.gcamp,5);
    end
    if min(gcamp_bag)<0
        gcamp_bag=gcamp_bag+abs(min(gcamp_bag));
    end
    
    ratio_bag=(gcamp_bag/nanmean(gcamp_bag))./(cherry_bag/nanmean(cherry_bag));
    %     eq_ratio_cent=eq_ratio-nanmedian(eq_ratio);
    %     eq_ratio_norm=smoothn(eq_ratio_cent/nanstd(eq_ratio_cent),10);
    
    %R/R0: divide by the lower percentile
    sr=sort(ratio_bag);
    R0_bag=nanmean(sr(1:length(ratio_bag)/10));
    if R0_bag<0.15
        R0_bag=nanmean(sr(1:length(ratio_bag)/8));
    end
    
    ratio=ratio_bag./R0_bag;
    
    
    %reversals:
    if isfield(ImagingData{F},'RevFrames30hz')==1
        
        RevON=ImagingData{F}.RevFrames30hz(1:2:end);
        RevEND=ImagingData{F}.RevFrames30hz(2:2:end);
        if isempty(RevEND) | length(RevEND)<2
            prerev_n_MoM(F,:)=NaN;
            continue
        end
        if RevON(1)<0
            RevON(1)=1;
        end
        
        if RevEND(end)<RevON(end)
            RevEND=[RevEND(:) ; length(ratio)];
        end
        if RevEND(end)>length(ratio)
            RevEND(end)=length(ratio);
        end
        
        for rev=1:length(RevON)
            if RevEND(rev)-RevON(rev)>minlength
                
                if RevON(rev)>prewin+1 & RevON(rev)+300<length(ratio)
                    prerev(cc,1:1800)=NaN(1,1800);
                    cr=ratio(RevON(rev)-prewin:RevEND(rev));
                    prerev(cc,1:length(cr))=cr;
                    cc=cc+1;
                    % if reversal is too close to start of trace:
                elseif RevON(rev)<prewin
                    
                    gap=abs(RevON(rev)-prewin+1);
                    prerev(cc,1:1800)=NaN(1,1800);
                    cr=ratio(1:RevEND(rev));
                    prerev(cc,gap:gap+length(cr)-1)=cr;
                    cc=cc+1;
                    
                end
            end
            
        end % end reversal loop
        
        %% % dC/dT:
        
        prerev=prerev(:,1:winsize);
        
        %NaN episodes of strong rise at the end of reversal (->self-touch)
        for i=1:size(prerev,1)
            [ns ~]=nanedge(prerev(i,:));
            if isempty(ns)
                ns=size(prerev,2)-3;
            end
            ns=ns(end)-1;
            if ns>prewin+1+(30*6)
                sp=diff(diff(medfilt1(prerev(i,:),50)));
                if find(sp(ns-90:ns)>0.04)
                    prerev(i,ns-90:ns)=NaN;
                end
                
            end
        end
        
        
        % normalize:
        prerev_n=NaN(size(prerev));
        for i=1:cc-1
            mv=nanmean(prerev(i,prewin-hz-5:prewin-5));
            prerev_n(i,:)=medfilt1(prerev(i,:),5)./mv;
        end
        
        %average each reversal over sec 1-4:
        ns1_4=nanmean(prerev_n(:,prewin+(hz/5):prewin+(hz*4)),2);
        
        revAll(1:length(ns1_4),F)=ns1_4;
        
    end
end

revAll(revAll==0)=NaN;

RevDataNested.Data=revAll;
RevDataNested.Date=datestr(now);
st = dbstack;
RevDataNested.Code=st.name;
save RevDataNested RevDataNested