%(c) Julia Riedl 2017 
%separates stage position data retrieved from Metamorph into separate files
% for each imaging experiment.
%Input: data=csvimport_stagepos('Preleth_br_20_ch_1stagepos.csv');
clearvars -except data
cc=1;
str1=data{cc,1}(1:end);
XY=NaN(1,2);

for i=1:length(data)
    
    str2=data{i,1}(1:end);
    
    if strcmp(str1,str2)
        
        X=data{i,3};
        
        if length(X)>1
            if X(end-1)==','
                X(end-1)='.';
            end
        end
        if ~isempty(data{i,end})
            Y=data{i,4};
        else
            XY(cc,:)=NaN;
            cc=cc+1;
            continue
            
        end
        
        if length(Y)>1 & Y(end-1)==','
            Y(end-1)='.';
        end
        if ischar(X)
            bi=strfind(X,'"');
            X(bi)=[];
            X=str2num(X);
            if ~isnan(X)
                X=X(1);
            end
            
            bi=strfind(Y,'"');
            Y(bi)=[];
            Y=str2num(Y);
            if ~isnan(Y)
            Y=Y(1);
            end
        end
        XY(cc,:)=[X Y];
        if XY(cc,1)==XY(cc,2)
            i;
        end
        cc=cc+1;
    else
        disp(str1)
        ci=strfind(str1,'ch');
        if ~isempty(ci)
           save([str1(3:end-1) '_stagepos.mat'], 'XY');
        end
        
        str1=data{i,1}(1:end);
        cc=1;
        XY=NaN(1,2);
%         %plot
%         plot(XY(:,1))
%         hold on
%         plot(XY(:,2),'r');

    end
    if mod(i,7000)==0
        i
    end
    
end



