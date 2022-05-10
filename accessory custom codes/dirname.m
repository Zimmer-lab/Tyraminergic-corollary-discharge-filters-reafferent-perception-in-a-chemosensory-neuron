%(c) Julia Riedl 2016
%returns the name of the current folder
% as a string

function name=dirname(cd)

        nd=(cd);
        d= strfind(nd, '\');
        name=nd(d(end)+1:end);
end