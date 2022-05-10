%(c) Julia Riedl 2016
%returns the name of the current top folder
%as a string
function name=dirname2(cd)

        nd=(cd);
        d= strfind(nd, '\');
        name=nd(d(end-1)+1:d(end)-1);
end