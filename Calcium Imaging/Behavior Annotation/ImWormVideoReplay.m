

% (c) Julia Riedl 2016
%replays worm infrared movie for annotation of reversals
%requires video software suite(VideoUtils v1.2.4; Copyright (c) 2012, Marc Vivet)
% and avi files as input
%movienumber= avi movie in current location, rate: the higher the faster is
%the replay (6-9 is average), start is start of for replay
%Replays the infrared worm video and provides a button for annotating reversal starts
%and ends which is output as 'RevFrames' variable as one vector. 
function fn=ImWormVideoReplay(movienumber,rate,start)

global fn
warning off

mnames=dir('*.avi');
moviename=mnames(movienumber).name;

if rate>10   
   vp = VideoPlayer(moviename, 'Verbose', false, 'ShowTime', false,'InitialFrame',start,'StepInFrames',2);
 
else
vp = VideoPlayer(moviename, 'Verbose', false, 'ShowTime', false,'InitialFrame',start);
end
screen=get(0,'screensize');
setPlotPosition(vp, [screen(3)/2.3,screen(4)/3.2,screen(3)/2,screen(4)/1.8]);

%annotation button:
AnnotateFrameRev(moviename,start)

while ( true )

        plot( vp );
        drawnow;
        fn=vp.FrameNum;
        pause(0.1/rate)

        if ( ~vp.nextFrame)
            break;
        end
     
end

end