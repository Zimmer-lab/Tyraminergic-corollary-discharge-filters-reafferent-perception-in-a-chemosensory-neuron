%(c) Julia Riedl 2017; code finds staring & end points of stretches of nan in a vector

function [NaNStart, NaNEnd]=nanedge(data)


%(1) create avector with 1=NaN
nanvector=zeros(1,length(data));
nanvector(isnan(data))=1;

%(2)find onsets of NaN episodes:
NaNStart= find(diff(nanvector)==1);
NaNEnd= find(diff(nanvector)==-1);

%add first point as start and last point as end if NaN:
if nanvector(1)==1
    NaNStart=[1 NaNStart];
end
if nanvector(end)==1
    NaNEnd=[ NaNEnd length(nanvector)];
end

end
