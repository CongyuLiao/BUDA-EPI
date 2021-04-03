function [mdh, scan_num] = read_meas_dat__mdh_binary(binary, varargin)
%READ_MEAS_DAT__MDH_BINARY  quickly parse MDH extracted as one binary vector
%
% mdh_struct = read_meas_dat__mdh_binary(binary)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/mar/24
% $Id: read_meas_dat__mdh_binary__alt.m,v 1.1 2010/12/07 14:14:12 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  request_full = 0;
  if ( nargin >= 2 ), 
    request_full = varargin{1};
  end;

  isvd = (length(binary) == 192);	% Is this a VD11 data file?

  if (isvd)				% VD11 or newer MDH structure
    %==----------------------------------------------------------------------==%

    % constants defined in <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>
    MDH_NUMBEROFEVALINFOMASK   = 2;
    MDH_NUMBEROFICEPROGRAMPARA = 24;

    MDH_FREEHDRPARA = 4;


    %==----------------------------------------------------------------------==%

    ind = 1;

    ulDMALength = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    ulFlags1    = typecastc(binary(ind:ind+0), 'uint8');
    ind = ind + 1;
    ulFlags2    = typecastc(binary(ind:ind+0), 'uint8');
    ind = ind + 1;

    mdh.ulFlagsAndDMALength        = double(ulFlags2) * 2^24 + ...
                                     uint32(double(ulFlags1) * 2^16) + ...
                                     uint32(ulDMALength);

    mdh.lMeasUID                   = typecastc(binary(ind:ind+3), 'int32');
    ind = ind + 4;
    mdh.ulScanCounter              = typecastc(binary(ind:ind+3), 'uint32');
    ind = ind + 4;
    mdh.ulTimeStamp                = typecastc(binary(ind:ind+3), 'uint32')*2.5;
    ind = ind + 4; % milliseconds
    mdh.ulPMUTimeStamp             = typecastc(binary(ind:ind+3), 'int32')*2.5;
    ind = ind + 4; % milliseconds
    mdh.ushSystemType              = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if ( request_full),
      mdh.ushPTABPosDelay          = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.lPTABPosX                = typecastc(binary(ind:ind+3), 'int32');
      ind = ind + 4;
      mdh.lPTABPosY                = typecastc(binary(ind:ind+3), 'int32');
      ind = ind + 4;
      mdh.lPTABPosZ                = typecastc(binary(ind:ind+3), 'int32');
      ind = ind + 4;
      mdh.ulReserved1              = typecastc(binary(ind:ind+3), 'uint32');
      ind = ind + 4;
    else
      ind = ind + 18;
    end;

    mdh.aulEvalInfoMask(1:MDH_NUMBEROFEVALINFOMASK) = typecastc(binary(ind:ind+4*MDH_NUMBEROFEVALINFOMASK-1), 'uint32');
    ind = ind + 4*MDH_NUMBEROFEVALINFOMASK;

    % build 64-bit mask from two 32-bit integers
    mask = double(...
        bitshift(uint64(mdh.aulEvalInfoMask(2)), 32)) ...
                  + double(mdh.aulEvalInfoMask(1));

    mdh.ushSamplesInScan           = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushUsedChannels            = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if (1),
      %  sLoopCounter
      mdh.ushLine                  = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushAcquisition           = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSlice                 = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPartition             = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushEcho                  = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPhase                 = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushRepetition            = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSet                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSeg                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIda                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdb                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdc                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdd                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIde                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
    end;

    if (1),
      % sCutOffData
      mdh.ushPre                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPost                  = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
    end;

    mdh.ushKSpaceCentreColumn      = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushCoilSelect              = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.fReadOutOffcentre          = typecastc(binary(ind:ind+3), 'single');
    ind = ind + 4;
    mdh.ulTimeSinceLastRF          = typecastc(binary(ind:ind+3), 'uint32');
    ind = ind + 4;
    mdh.ushKSpaceCentreLineNo      = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushKSpaceCentrePartitionNo = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if (1),
      % sSliceData
      if (1),
        % sVector
        mdh.flSag                  = typecastc(binary(ind:ind+3), 'single');
        ind = ind + 4;
        mdh.flCor                  = typecastc(binary(ind:ind+3), 'single');
        ind = ind + 4;
        mdh.flTra                  = typecastc(binary(ind:ind+3), 'single');
        ind = ind + 4;
      end;
      if ( request_full), 
        mdh.aflQuaternion(1:4)     = typecastc(binary(ind:ind+15), 'single'); 
      end;
      ind = ind + 16;
    end;

    if ( request_full), 
      mdh.aushIceProgramPara(1:MDH_NUMBEROFICEPROGRAMPARA) = ...
        typecastc(binary(ind:ind+2*MDH_NUMBEROFICEPROGRAMPARA-1), 'uint16'); 
      ind = ind + 2*MDH_NUMBEROFICEPROGRAMPARA;
  
      mdh.aushFreePara(1:MDH_FREEHDRPARA) = ...
        typecastc(binary(ind:ind+2*MDH_FREEHDRPARA-1), 'uint16');
      ind = ind + 2*MDH_FREEHDRPARA;
    else
      ind = ind + 2*(MDH_NUMBEROFICEPROGRAMPARA + MDH_FREEHDRPARA);
    end;

    mdh.ushApplicationCounter      = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushApplicationMask         = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ulCRC                      = typecastc(binary(ind:ind+3), 'uint32');
    ind = ind + 4;

    mdh.ushChannelId               = [];	% For backwards compatibility

    %==----------------------------------------------------------------------==%
  else					% VB17 or older MDH structure
    %==----------------------------------------------------------------------==%

    % constants defined in <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>
    MDH_NUMBEROFEVALINFOMASK   = 2;
    MDH_NUMBEROFICEPROGRAMPARA = 4;

    MDH_FREEHDRPARA = 4;


    %==----------------------------------------------------------------------==%

    ind = 1;

    ulDMALength = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    ulFlags1    = typecastc(binary(ind:ind+0), 'uint8');
    ind = ind + 1;
    ulFlags2    = typecastc(binary(ind:ind+0), 'uint8');
    ind = ind + 1;

    mdh.ulFlagsAndDMALength        = double(ulFlags2) * 2^24 + ...
                                     uint32(double(ulFlags1) * 2^16) + ...
                                     uint32(ulDMALength);

    mdh.lMeasUID                   = typecastc(binary(ind:ind+3), 'int32');
    ind = ind + 4;
    mdh.ulScanCounter              = typecastc(binary(ind:ind+3), 'uint32');
    ind = ind + 4;
    mdh.ulTimeStamp                = typecastc(binary(ind:ind+3), 'uint32')*2.5;
    ind = ind + 4; % milliseconds
    mdh.ulPMUTimeStamp             = typecastc(binary(ind:ind+3), 'int32')*2.5;
    ind = ind + 4; % milliseconds

    mdh.aulEvalInfoMask(1:MDH_NUMBEROFEVALINFOMASK) = typecastc(binary(ind:ind+4*MDH_NUMBEROFEVALINFOMASK-1), 'uint32');
    ind = ind + 4*MDH_NUMBEROFEVALINFOMASK;

    % build 64-bit mask from two 32-bit integers
    mask = double(...
        bitshift(uint64(mdh.aulEvalInfoMask(2)), 32)) ...
                  + double(mdh.aulEvalInfoMask(1));

    mdh.ushSamplesInScan           = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushUsedChannels            = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if (1),
      %  sLoopCounter
      mdh.ushLine                  = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushAcquisition           = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSlice                 = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPartition             = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushEcho                  = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPhase                 = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushRepetition            = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSet                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSeg                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIda                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdb                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdc                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdd                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIde                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
    end;

    if (1),
      % sCutOffData
      mdh.ushPre                   = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPost                  = typecastc(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
    end;

    mdh.ushKSpaceCentreColumn      = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushDummy                   = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.fReadOutOffcentre          = typecastc(binary(ind:ind+3), 'single');
    ind = ind + 4;
    mdh.ulTimeSinceLastRF          = typecastc(binary(ind:ind+3), 'uint32');
    ind = ind + 4;
    mdh.ushKSpaceCentreLineNo      = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushKSpaceCentrePartitionNo = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if ( request_full), 
      mdh.aushIceProgramPara(1:MDH_NUMBEROFICEPROGRAMPARA) = ...
        typecastc(binary(ind:ind+2*MDH_NUMBEROFICEPROGRAMPARA-1), 'uint16'); 
      ind = ind + 2*MDH_NUMBEROFICEPROGRAMPARA;
  
      mdh.aushFreePara(1:MDH_FREEHDRPARA) = ...
        typecastc(binary(ind:ind+2*MDH_FREEHDRPARA-1), 'uint16');
      ind = ind + 2*MDH_FREEHDRPARA;
    else
      ind = ind + 2*(MDH_NUMBEROFICEPROGRAMPARA + MDH_FREEHDRPARA);
    end;

    if (1),
      % sSliceData
      if (1),
        % sVector
        mdh.flSag                  = typecastc(binary(ind:ind+3), 'single');
        ind = ind + 4;
        mdh.flCor                  = typecastc(binary(ind:ind+3), 'single');
        ind = ind + 4;
        mdh.flTra                  = typecastc(binary(ind:ind+3), 'single');
        ind = ind + 4;
      end;
      if ( request_full), 
        mdh.aflQuaternion(1:4)     = typecastc(binary(ind:ind+15), 'single'); 
      end;
      ind = ind + 16;
    end;

    mdh.ushChannelId               = typecastc(binary(ind:ind+1), 'uint16');
    ind = ind+2;

    if ( request_full),
      mdh.ushPTABPosNeg            = typecastc(binary(ind:ind+1), 'uint16');
    end;
    ind = ind+2;

    %==----------------------------------------------------------------------==%
  end

  scan_num = double(mdh.ulScanCounter);
  
  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_meas_dat__mdh_binary__alt.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

