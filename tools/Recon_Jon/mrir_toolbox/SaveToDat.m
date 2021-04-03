
function SaveToDat(data, fname)

    fid = fopen(fname,'w');
    fwrite(fid,size(data),'single');
    fclose(fid);
    
    fid = fopen(fname,'a+');
    fwrite(fid,single(real(data)),'single');
    fwrite(fid,single(imag(data)),'single');
    fclose(fid);
    
end