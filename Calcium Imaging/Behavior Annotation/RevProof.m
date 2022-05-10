%(c) Julia Riedl 2017

% plots stageposition and annotated reversals corrected for delay and returns RevFrames in 10 and 30hz format
%reversals onsets are marked in red, reversal ends in green.

%delay input introduces a shift for 30 hz data to correct for delays in
%annotation or consistent video onset shift.


function [RevFrames30hz, recordname]=RevProof(F,RevFrames,delay)


global fn click handles

start=10;

click=0;
warning off
%close all
d=NaN;

%XY:
posfiles=dir('*stagepos*');
XYpos=load(posfiles(F).name);
try
    XY=XYpos.XY;
end
X=1-XY(1:5:end,1);
Y=XY(1:5:end,2);
Ya=XY(:,1);
Xa=1-XY(:,2);

%ratio:
rfiles=dir(['*_log.mat']);
load (rfiles(F).name);

%name:
stacks=dir('*ch*.stk'); stackname=stacks(F).name;
us= strfind(stackname, '.');
recordname=[stackname(1:us(end)-1) '_rev_SF_1'];

screen=get(0,'screensize');
mf=figure('Position', [screen(3)/3,screen(4)/3,screen(3)/1.8,screen(4)/1.9],'Name',stackname);

plot(Xa,Ya,'k');
axis equal
title(delay)
hold on
RevFrames30hz=round((RevFrames-delay)*3);

RevON=RevFrames30hz(1:2:end);
RevEND=RevFrames30hz(2:2:end);
% delete reversals which happened after recording end:
os=find(RevON>length(Xa));
RevON(os)=[];
os=find(RevEND>length(Xa));
RevEND(os)=[];

scatter(Xa(RevON(1)),Ya(RevON(1)),'bp','filled');
try
    scatter(Xa(RevON(2:end)),Ya(RevON(2:end)),'rp','filled');
    scatter(Xa(RevON(end)),Ya(RevON(end)),'mh','filled');
catch
    scatter(Xa(RevON(2:end-1)),Ya(RevON(2:end-1)),'rp','filled');
end


try
    scatter(Xa(RevEND(1:end)),Ya(RevEND(1:end)),'gp','filled')
catch
    scatter(Xa(RevEND(1:end-1)),Ya(RevEND(1:end-1)),'gp','filled')
end


end


