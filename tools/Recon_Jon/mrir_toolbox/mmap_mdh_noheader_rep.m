function m = mmap_mdh_noheader_rep(filename,linesperrep,nRep)
%
% mmap_mdh
%
% 2009 by Thomas Witzel
% Massachusetts Institute of Technology
%
% 29 Oct 2009 Modified by Thomas Witzel and Kawin Setsompop
% 
% This file contains proprietary information of Siemens Healthcare
% Solutions, DO NOT DISTRIBUTE !!!!
%
fid = fopen(filename,'r');

% first read the offset into the file
offset = fread(fid,1,'uint32');
%disp(sprintf('offset = %d',offset));

% then figure out how many samples per readout there should be
fseek(fid,offset,'bof');
fseek(fid,28,'cof');
NSAMP = fread(fid,1,'uint16');
NChannels = fread(fid,1,'uint16');
disp(sprintf('nsamples in k-space line = %d',NSAMP));

% now figure out the end of the file 
fseek(fid,0,'eof');
endpos = ftell(fid);
fclose(fid);

offset = endpos - (linesperrep*(nRep)*(NSAMP*8+128) + NChannels*(16*8+128)); % ADCEND data for each channel at the end (16 samples + header)

m = memmapfile(filename,'Offset', offset,        ...
    'Format', {'single' ,[1 linesperrep*(32+NSAMP*2)], 'kdata'}, ...
    'Repeat',nRep);




% 
% % now figure out how many k-space lines there are
% fseek(fid,0,'eof');
% endpos = ftell(fid);
% dsize = endpos-offset;
% fclose(fid);
% 
% % this calculation assumes that the scan was completed and there are
% % ADCEND data for each channel (16 samples + header)
% 
% nlines = (dsize-NChannels*(16*8+128))/(NSAMP*8+128)
% nrepsinfile = floor(nlines/linesperrep)
% nlines_firstrep = linesperrep + rem(nlines,linesperrep)
% 
% offset = offset+nlines_firstrep*(NSAMP*8+128)
% %disp(sprintf('nreps in k-space data = %ld',nrepsinfile));
% 
% offset2 = endpos - (linesperrep*(nRep-1)*(NSAMP*8+128) + NChannels*(16*8+128));
% 
% m = memmapfile(filename,'Offset', offset,        ...
%     'Format', {'single' ,[1 linesperrep*(32+NSAMP*2)], 'kdata'}, ...
%     'Repeat',nRep-1);
