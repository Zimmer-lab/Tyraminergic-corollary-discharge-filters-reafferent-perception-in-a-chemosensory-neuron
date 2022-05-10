%(c) Julia Riedl 2018
%gather all imaging Data (*log.mat files) into one struct file "ImagingDataAll":

clear
files = dir('*BAG*');
DayFolders=files(find(vertcat(files.isdir)));
load CO2fitUS %load gradient fit vector
CaData=cell(1,1);

%gradient:
gradient=vertcat(CO2fitUS,CO2fitUS,CO2fitUS);
gradient1=gradient;
for i=1:12
    gradient=vertcat(gradient,gradient(1:end/1.05,:));   
end


home=cd;
cc=1;
cc2=0; %number of worms

for D=1:length(DayFolders)-0
    
    if isempty(strfind(DayFolders(D).name,'o'))
        
        cd(DayFolders(D).name)
        disp(DayFolders(D).name)
        
        worms=dir('W*');
        
        
        for W=1:length(worms)
            
            cd(worms(W).name)
            cc2=cc2+1;
            
            revfiles=dir('*SF*');
            rev_end_files=dir('*rev_end*');
            
            if ~isempty(revfiles)% & ~isempty(rev_end_files)
                rev=1;
            else
                rev=0;
            end
            
            %ratio:
            files=dir('*R1*_log.mat');
            if isempty(files)
                disp([worms(W).name ': no log file'])
                cd ..\
                continue
            end

            % go through various movies for this worm:
            for movie = 1:length(files)
                
                disp(files(movie).name)
                load (files(movie).name);
                
                gcamp=thiswormdata(:,1);
                cherry=thiswormdata(:,2);
                
                ratio=(thiswormdata(:,1)./thiswormdata(:,2));
                % find very rapid changes in size of cherry area and remove this data:
                areaR1=thiswormdata(:,3);
                areaJumps1=find(diff(medfilt1(areaR1))>300 | diff(medfilt1(areaR1))<-300);
                areaJumps1=vertcat(areaJumps1, find(areaR1>2*nanmean(areaR1) | areaR1<0.25*nanmean(areaR1)));
                ratio(areaJumps1+1)=NaN;
                ratio(ratio>2.5| ratio<-0.9)=NaN;
                
                
                %XY:
                posfiles=dir('*_stagepos*');
                disp(posfiles(movie).name);
                XYpos=load(posfiles(movie).name);
                if isstruct(XYpos)
                    XYpos=XYpos.XY;
                end
                X=1-XYpos(:,1);
                Y=XYpos(:,2);
                
                %normalize and round:
                X=round((X+33500)/10);
                Y=round((Y+23500)/10);
                if max(X)>size(gradient,2)
                    X=X-(max(X)-size(gradient,2))-10;
                end
                    
                %Map XY on gradient:                
                X(X>length(gradient))=NaN;
                C_ind=sub2ind(size(gradient),Y,X);
                
                %get CO2 concentration at this position of path:
                nan_idx=isnan(C_ind);
                C_ind=C_ind(nan_idx~=1);
                
                %get CO2 concentration at this position of path and reinsert Nans:
                gradV=(gradient(round(C_ind)));
                gradVal=NaN(1,length(nan_idx));
                gradVal(nan_idx~=1)=gradV;
                
                %annotated reversals:
                if rev==1 & length(revfiles)==length(files)
                    load (revfiles(movie).name)
                    if movie==1
                        rev1=RevFrames30hz;
                    elseif movie==2
                        rev1=[rev1 RevFrames30hz+length(gcamp_o)];
                    end
                     CaData{cc}.RevFrames30hz=rev1;
                   
                end
                %check for data length:
                 [ns,~]=nanedge(thiswormdata(:,1));
                 ns=ns-1;
                 if ~isempty(ns) & ns(end)>length(X)
                     disp(W)
                     disp('inconsistent XY and R length!')
                 end
                 if length(gcamp)<length(X)
                     X=X(1:length(gcamp));
                     Y=Y(1:length(gcamp));
                 end
                 
                 if movie==1
                     gcamp_o=gcamp(1:length(X));
                     cherry_o=cherry(1:length(X));
                     gradVal_1=gradVal;
                     X1=X;
                     Y1=Y;
                 elseif movie==2
                     gcamp_o=[gcamp_o;gcamp(1:length(X))];
                     cherry_o=[cherry_o;cherry(1:length(X))];
                     gradVal_1=[gradVal_1(:); gradVal(:)];
                     X1=[X1;X];
                     Y1=[Y1;Y];
                 end
                 
            end
            
            CaData{cc}.cherry=cherry_o;
            CaData{cc}.gcamp=gcamp_o;
            CaData{cc}.CO2=gradVal_1;
            CaData{cc}.XY=[X1 Y1];
            CaData{cc}.TrialLabel=posfiles(movie).name(1:end-16);
            
            cc=cc+1;
            
            cd ..\
            
        end
        
    end % end Worm Loop
    
    cd(home)
    
end

%% save:
ImagingData=CaData;
save ImagingDataAll ImagingData






