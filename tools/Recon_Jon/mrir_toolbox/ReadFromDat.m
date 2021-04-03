function data = ReadFromDat(fname, dims)

if nargin == 1
    dims = 10;
end

fid = fopen(fname,'r');
data = single(fread(fid,'single'));
fclose(fid);
s = data(1:dims);
data = reshape(data(1+dims:dims+(end-dims)/2) + 1i*data(1+dims+(end-dims)/2:end),s');

end

