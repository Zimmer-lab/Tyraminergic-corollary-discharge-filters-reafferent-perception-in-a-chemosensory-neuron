function background = getbackgroundMZ64_MSv101 (FileInfo)

%-- this function adds every "frameinterval" frames and averages to obtain 
%-- the background image
%-----------------------------------------------------------------------
frameinterval = 10;
m = FileInfo.Width;
n = FileInfo.Height;

if FileInfo.Channels == 3
    cdatasum = zeros(n,m,3,'double'); %for 24-bit movies
else
    cdatasum = zeros(n,m,'double');   %for 8-bit movies
end;
numFrames =200% FileInfo.NumFrames;

FileInfo.nextFrame();
for Frame = 1:frameinterval:numFrames
    MovX64 = double(FileInfo.getFrameUInt8())/255;
    cdatasum = cdatasum + MovX64;
    FileInfo + frameinterval;
end

cdataaverage = cdatasum./round(numFrames/frameinterval);
background = uint8(round(cdataaverage*255));
