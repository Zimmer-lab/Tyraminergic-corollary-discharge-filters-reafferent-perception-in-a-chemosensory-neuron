function background = getbackgroundMZ64 (RawMovieName,frameinterval)

%this function adds every "frameinterval" frames and averages to obtain the background
%image
%-----------------------------------------------------------------------
 
FileInfo = VideoReader(RawMovieName);
m = FileInfo.Width;
n = FileInfo.Height;;

if strcmp(FileInfo.VideoFormat,'RGB24')
    
    cdatasum = zeros(n,m,3,'double'); %for 24-bit movies
    
else
    cdatasum = zeros(n,m,'double');   %for 8-bit movies
    
end;

% Mov = aviread(RawMovieName, 1);
 %Movcolormap = Mov.colormap;

for Frame = 1:frameinterval:FileInfo.NumberOfFrames
    Mov = read(FileInfo, Frame);
    MovX64 = double(Mov)/255;
    cdatasum = cdatasum + MovX64;  
end

cdataaverage = cdatasum./round(FileInfo.NumberOfFrames/frameinterval);
background = uint8(round(cdataaverage*255));

