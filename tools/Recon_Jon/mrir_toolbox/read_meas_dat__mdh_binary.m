function [mdh, scan_num] = read_meas_dat__mdh_binary(binary)
%READ_MEAS_DAT__MDH_BINARY  quickly parse MDH extracted as one binary vector
%
% mdh_struct = read_meas_dat__mdh_binary(binary)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/mar/24
% $Id: read_meas_dat__mdh_binary.m,v 1.1 2009/01/09 20:30:16 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % constants defined in <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>
  MDH_NUMBEROFEVALINFOMASK   = 2;
  MDH_NUMBEROFICEPROGRAMPARA = 4;

  MDH_FREEHDRPARA = 4;


  %==-----------------------------------------------------------------------==%

  ind = 1;

  ulDMALength = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;

  ulFlags1    = uint8(typecast(binary(ind:ind+0), 'uint8')); ind = ind + 1;
  ulFlags2    = uint8(typecast(binary(ind:ind+0), 'uint8')); ind = ind + 1;

  mdh.ulFlagsAndDMALength        = uint32(double(ulFlags2) * 2^24) + uint32(double(ulFlags1) * 2^16) + uint32(ulDMALength);

  mdh.lMeasUID                   = int32(typecast(binary(ind:ind+3), 'int32')); ind = ind + 4;
  mdh.ulScanCounter              = uint32(typecast(binary(ind:ind+3), 'uint32')); ind = ind + 4;
  mdh.ulTimeStamp                = uint32(typecast(binary(ind:ind+3), 'uint32')*2.5); ind = ind + 4; % milliseconds
  mdh.ulPMUTimeStamp             = uint32(typecast(binary(ind:ind+3), 'int32')*2.5); ind = ind + 4; % milliseconds

  mdh.aulEvalInfoMask(1:MDH_NUMBEROFEVALINFOMASK) = uint32(typecast(binary(ind:ind+4*MDH_NUMBEROFEVALINFOMASK-1), 'uint32')); ind = ind + 4*MDH_NUMBEROFEVALINFOMASK;

  % build 64-bit mask from two 32-bit integers
  mask = uint64(double(...
      bitshift(uint64(mdh.aulEvalInfoMask(2)), 32)) ...
                + double(mdh.aulEvalInfoMask(1)));

  mdh.ushSamplesInScan           = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
  mdh.ushUsedChannels            = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;

  if (1),
    %  sLoopCounter
    mdh.ushLine                    = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushAcquisition             = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushSlice                   = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushPartition               = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushEcho                    = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushPhase                   = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushRepetition              = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushSet                     = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushSeg                     = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushIda                     = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushIdb                     = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushIdc                     = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushIdd                     = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushIde                     = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
  end;

  if (1),
    % sCutOffData
    mdh.ushPre                     = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
    mdh.ushPost                    = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
  end;

  mdh.ushKSpaceCentreColumn      = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
  mdh.ushDummy                   = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
  mdh.fReadOutOffcentre          = single(typecast(binary(ind:ind+3), 'single')); ind = ind + 4;
  mdh.ulTimeSinceLastRF          = uint32(typecast(binary(ind:ind+3), 'uint32')); ind = ind + 4;
  mdh.ushKSpaceCentreLineNo      = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;
  mdh.ushKSpaceCentrePartitionNo = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind + 2;

  mdh.aushIceProgramPara(1:MDH_NUMBEROFICEPROGRAMPARA) = uint16(typecast(binary(ind:ind+2*MDH_NUMBEROFICEPROGRAMPARA-1), 'uint16')); ind = ind+2*MDH_NUMBEROFICEPROGRAMPARA;
  mdh.aushFreePara(1:MDH_FREEHDRPARA)                  = uint16(typecast(binary(ind:ind+2*MDH_FREEHDRPARA-1), 'uint16')); ind = ind+2*MDH_FREEHDRPARA;

  if (1),
    % sSliceData
    if (1),
      % sVector
      mdh.flSag            = single(typecast(binary(ind:ind+3), 'single')); ind = ind+4;
      mdh.flCor            = single(typecast(binary(ind:ind+3), 'single')); ind = ind+4;
      mdh.flTra            = single(typecast(binary(ind:ind+3), 'single')); ind = ind+4;
    end;
    mdh.aflQuaternion(1:4) = single(typecast(binary(ind:ind+(4*4)-1), 'single')); ind = ind+(4*4);
  end;

  mdh.ushChannelId  = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind+2;
  mdh.ushPTABPosNeg = uint16(typecast(binary(ind:ind+1), 'uint16')); ind = ind+2;


  %==-----------------------------------------------------------------------==%

  scan_num = double(mdh.ulScanCounter);
  
  
  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_meas_dat__mdh_binary.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

