function mdh = read_meas_dat__mdh(fp)
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/mar/24
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % constants defined in <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>
  MDH_NUMBEROFEVALINFOMASK   = 2;
  MDH_NUMBEROFICEPROGRAMPARA = 4;

  MDH_FREEHDRPARA = 4;

  
  %==-----------------------------------------------------------------------==%
  
  
  ulDMALength = uint16(fread(fp, 1, 'uint16'));
  
  ulFlags1    = uint8(fread(fp, 1, 'uint8'));
  ulFlags2    = uint8(fread(fp, 1, 'uint8'));
  
  mdh.ulFlagsAndDMALength        = uint32(double(ulFlags2) * 2^24) + uint32(double(ulFlags1) * 2^16) + uint32(ulDMALength);

  mdh.lMeasUID                   = int32(fread(fp, 1, 'int32'));
  mdh.ulScanCounter              = uint32(fread(fp, 1, 'uint32'));
  mdh.ulTimeStamp                = uint32(fread(fp, 1, 'uint32')*2.5); % milliseconds
  mdh.ulPMUTimeStamp             = uint32(fread(fp, 1, 'uint32')*2.5); % milliseconds

  mdh.aulEvalInfoMask(1:MDH_NUMBEROFEVALINFOMASK) = uint32(fread(fp, MDH_NUMBEROFEVALINFOMASK, 'uint32'));

  % build 64-bit mask from two 32-bit integers
  mask = uint64(double(...
      bitshift(uint64(mdh.aulEvalInfoMask(2)), 32)) ...
		+ double(mdh.aulEvalInfoMask(1)));

  mdh.ushSamplesInScan           = uint16(fread(fp, 1, 'uint16'));
  mdh.ushUsedChannels            = uint16(fread(fp, 1, 'uint16'));

  if (1),
    %  sLoopCounter
    mdh.ushLine                    = uint16(fread(fp, 1, 'uint16'));
    mdh.ushAcquisition             = uint16(fread(fp, 1, 'uint16'));  % note: acquisition is same as average
    mdh.ushSlice                   = uint16(fread(fp, 1, 'uint16'));
    mdh.ushPartition               = uint16(fread(fp, 1, 'uint16'));
    mdh.ushEcho                    = uint16(fread(fp, 1, 'uint16'));
    mdh.ushPhase                   = uint16(fread(fp, 1, 'uint16'));
    mdh.ushRepetition              = uint16(fread(fp, 1, 'uint16'));
    mdh.ushSet                     = uint16(fread(fp, 1, 'uint16'));
    mdh.ushSeg                     = uint16(fread(fp, 1, 'uint16'));
    mdh.ushIda                     = uint16(fread(fp, 1, 'uint16'));
    mdh.ushIdb                     = uint16(fread(fp, 1, 'uint16'));
    mdh.ushIdc                     = uint16(fread(fp, 1, 'uint16'));
    mdh.ushIdd                     = uint16(fread(fp, 1, 'uint16'));
    mdh.ushIde                     = uint16(fread(fp, 1, 'uint16'));
  end;

  if (1),
    % sCutOffData
    mdh.ushPre                     = uint16(fread(fp, 1, 'uint16'));
    mdh.ushPost                    = uint16(fread(fp, 1, 'uint16'));
  end;

  mdh.ushKSpaceCentreColumn      = uint16(fread(fp, 1, 'uint16'));
  mdh.ushDummy                   = uint16(fread(fp, 1, 'uint16'));
  mdh.fReadOutOffcentre          = single(fread(fp, 1, 'float32'));
  mdh.ulTimeSinceLastRF          = uint32(fread(fp, 1, 'uint32'));
  mdh.ushKSpaceCentreLineNo      = uint16(fread(fp, 1, 'uint16'));
  mdh.ushKSpaceCentrePartitionNo = uint16(fread(fp, 1, 'uint16'));

  mdh.aushIceProgramPara(1:MDH_NUMBEROFICEPROGRAMPARA) = uint16(fread(fp, MDH_NUMBEROFICEPROGRAMPARA, 'uint16'));
  mdh.aushFreePara(1:MDH_FREEHDRPARA)                  = uint16(fread(fp, MDH_FREEHDRPARA, 'uint16'));

  if (1),
    % sSliceData
    if (1),
      % sVector
      mdh.flSag            = single(fread(fp, 1, 'float32'));
      mdh.flCor            = single(fread(fp, 1, 'float32'));
      mdh.flTra            = single(fread(fp, 1, 'float32'));
    end;
    mdh.aflQuaternion(1:4) = single(fread(fp, 4, 'float32'));
  end;

  mdh.ushChannelId  = uint16(fread(fp, 1, 'uint16'));
  mdh.ushPTABPosNeg = uint16(fread(fp, 1, 'uint16'));

  
  return;

  

  %************************************************************************%
  %%% $Source: /home/jonnyreb/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

