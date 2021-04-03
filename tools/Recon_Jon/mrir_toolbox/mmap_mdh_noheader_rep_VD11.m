function m = mmap_mdh_noheader_rep_VD11(filename,linesperrep,nRep,NSAMP,NChannels,data_begin)
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

isvd = 1

if isvd == 0  % VB17 version
    
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
 elseif isvd == 1 % VD version 
    
    % first read the offset into the file
    
%     fseek(fid, 16, 'bof');
%     header_start = fread(fid, 1, 'uint32');
%     fseek(fid, header_start, 'bof');
%     offset = header_start + fread(fid, 1, 'uint32');
%     
    offset = data_begin 
    mdh_length_float32 = 48;
    mdh_ch_length_float32 = 8;
    
    largeHeader = mdh_length_float32*4;
    smallHeader = mdh_ch_length_float32*4;
    
     
    
    m = memmapfile(filename,'Offset', offset,        ...
        'Format', {'single' ,[1 linesperrep*(mdh_length_float32/NChannels + mdh_ch_length_float32+ NSAMP*2)], 'kdata'}, ...
        'Repeat',nRep);  
elseif isvd == 2  % non working VD 11
    
    % first read the offset into the file
    
    fseek(fid, 16, 'bof');
    header_start = fread(fid, 1, 'uint32');
    fseek(fid, header_start, 'bof');
    offset = header_start + fread(fid, 1, 'uint32');
    
    % then figure out how many samples per readout there should be
    %fseek(fid,offset,'bof');
    %fseek(fid,124,'cof');
    %NSAMP = fread(fid,1,'uint16');
    %NChannels = fread(fid,1,'uint16');
%     NSAMP = 168;
%     NChannels = 64;
%     
    disp(sprintf('nsamples in k-space line = %d',NSAMP));
    
    % now figure out the end of the file
    fseek(fid,0,'eof');
    endpos = ftell(fid);
    fclose(fid);
    
    mdh_length_float32 = 48;
    mdh_ch_length_float32 = 8;
    
    largeHeader = mdh_length_float32*4;
    smallHeader = mdh_ch_length_float32*4;
    
    %chunk_of_crap_at_end = NChannels*(16*8+largeHeader);
    chunk_of_crap_at_end = largeHeader + NChannels*8;%  %NChannels*(16*8+ smallHeader)+largeHeader- 4*5*2*(168+24);
       
    % 5.5. for 2 rps  1x invivo
    % 8 for DSI invivo
        
% 6.5 for 2 reps or more data
%     10 for 1 rep data
    
   
%     CurrentEndPad = (32*2+mdh_ch_length_float32)*64 +mdh_length_float32
%     tooMuchby = (21+7*168)*2 +48*8
%     Desired = CurrentEndPad-tooMuchby
%     
    %chunk_of_crap_at_end = 1878*2 -(168+24)*3*4 - (95+24+38+40)*4;
    
    offset = endpos - ( linesperrep*(nRep)*(NSAMP*8 + largeHeader/NChannels + smallHeader) + chunk_of_crap_at_end ); % ADQEND data for each channel at the end (16 samples + header)
    
    rem(( linesperrep*(nRep)*(NSAMP*8 + largeHeader/NChannels + smallHeader) + chunk_of_crap_at_end )/4096,1)
    m = memmapfile(filename,'Offset', offset,        ...
        'Format', {'single' ,[1 linesperrep*(mdh_length_float32/NChannels + mdh_ch_length_float32+ NSAMP*2)], 'kdata'}, ...
        'Repeat',nRep);
end




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
