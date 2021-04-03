function m = mmap_mdh_noheader_rep_v2(filename,linesperrep,nRep,NSAMP,NChannels,data_begin,isvd)
%
% mmap_mdh
%
% 2009 by Thomas Witzel
% Massachusetts Institute of Technology
%
% 29 Oct 2009 Modified by Thomas Witzel and Kawin Setsompop
% 
% 20 Jan 2012 Totally rewriten and simplify everything to look for data begin so can read incomplete files 
%
% This file contains proprietary information of Siemens Healthcare
% Solutions, DO NOT DISTRIBUTE !!!!
%

if isvd == 1
    mdh_length_float32 = 48;
    mdh_ch_length_float32 = 8;
else
    mdh_length_float32 = 0;
    mdh_ch_length_float32 = 32;
end

m = memmapfile(filename,'Offset', data_begin,        ...
    'Format', {'single' ,[1 linesperrep*(mdh_length_float32/NChannels + mdh_ch_length_float32+ NSAMP*2)], 'kdata'}, ...
    'Repeat',nRep);

