%(c) Julia Riedl 2017
% saves plots into a "plot" directory

pdir=dir('*plots*');
if isempty(pdir)
    disp('no plots folder')
end

try
if  exist('trackfiles','var') 
    pfiles=trackfiles;
else
    pfiles=[];
end
wfiles=dir('*W*');
afiles=dir('*als.mat');
if ~isempty(pfiles) | ~isempty(afiles)
    fileflag=1;
elseif ~isempty(wfiles)
    fileflag=1;
else
    fileflag=0;
end
    
if isempty(pdir) & fileflag==1
    
    name=dirname(cd);
    mkdir([name '_plots'])
    pdir=dir('*plots*');
    disp('..plot into folder')
end

catch ME

end