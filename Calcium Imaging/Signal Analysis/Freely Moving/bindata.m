
%bins data by averaging and downsampling:

function B=bindata(X,binsize)

a=mod(length(X),binsize);
Xr=(reshape(X(1:end-a),binsize,[]));
B = nanmean(Xr);

end