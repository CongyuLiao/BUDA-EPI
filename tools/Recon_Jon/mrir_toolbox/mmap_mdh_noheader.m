function m = mmap_mdh_noheader(filename)
%
% mmap_mdh
%
% 2009 by Thomas Witzel
% Massachusetts Institute of Technology
%
% This file contains proprietary information of Siemens Healthcare
% Solutions, DO NOT DISTRIBUTE !!!!
%
fid = fopen(filename,'r');

% first read the offset into the file
offset = fread(fid,1,'uint32');
disp(sprintf('offset = %d',offset));

% then figure out how many samples per readout there should be
fseek(fid,offset,'bof');
fseek(fid,28,'cof');
NSAMP = fread(fid,1,'uint16');
NChannels = fread(fid,1,'uint16');
disp(sprintf('nsamples in k-space line = %d',NSAMP));

% now figure out how many k-space lines there are 
fseek(fid,0,'eof');
endpos = ftell(fid);
dsize = endpos-offset;
% this calculation assumes that the scan was completed and there are
% ADCEND data for each channel (16 samples + header)
nlinesinfile = (dsize-NChannels*(16*8+128))/(NSAMP*8+128);
disp(sprintf('nlines in k-space data = %ld',nlinesinfile));
nrepeats = nlinesinfile;
fclose(fid);

m = memmapfile(filename,'Offset', offset,        ...  
     'Format', {'single' ,[1 32+NSAMP*2], 'kdata'}, ...
     'Repeat',nrepeats);
