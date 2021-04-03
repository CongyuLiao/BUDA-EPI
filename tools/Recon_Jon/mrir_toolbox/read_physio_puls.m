function [puls, t, physio] = read_physio_puls(filename, varargin)
%READ_PHYSIO_PULS
%
% [puls, t, physio] = read_physio_puls(filename)
%
%
% example:
%
%  [puls, t, physio] = read_physio_puls('2008_02_08__ge_functionals__MID408.puls');
%
%
%  figure; axis; hold on; box on;
%  stem(t, physio.trigger_on*max(puls),  'g', 'Marker', 'none');
%  stem(t, physio.trigger_off*max(puls), 'r', 'Marker', 'none');
%  plot(t, puls, 'k.-');
%  xlabel('time (s)');
%  trig_interval = round(mean(diff(t(find(physio.trigger_on)))));
%  title(sprintf('PF: %d / min, PP: %d ms', trig_interval*60, trig_interval*1000));

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/jan/23
% $Id: read_physio_puls.m,v 1.2 2008/02/14 01:20:27 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  [pathstr, filestr] = fileparts(filename);

  [fp, errstr] = fopen(filename, 'r', 'l');
  if ( fp == -1 ),
    error(errstr);
  end;

  physio.filename = filename;
%  physio.readtime = datestr(clock, 'yyyy-mmm-dd HH:MM:SS');


  PHYSIO_EOF = 6003;
  PHYSIO_EOL = 5003;
  PHYSIO_TRIGGER_ON  = 5000;  % special value for on
  PHYSIO_TRIGGER_OFF = 6000;  % special value for off


  %==--------------------------------------------------------------------==%

  fseek(fp, 0, 'bof');
  header = fscanf(fp, '%d', 4);
  fpos1 = ftell(fp);

%  PhysioMethod = sprintf('0x%02d', header(1));
%
%  switch PhysioMethod,
%   case '0x01',
%    physio.lPmuECGModePub = 'METHOD_NONE';
%   case '0x02',
%    physio.lPmuECGModePub = 'METHOD_TRIGGERING';
%   case '0x04',
%    physio.lPmuECGModePub = 'METHOD_GATING';
%   case '0x08',
%    physio.lPmuECGModePub = 'METHOD_RETROGATING';
%   case '0x10',
%    physio.lPmuECGModePub = 'METHOD_SOPE';
%   case '0x1E',
%    physio.lPmuECGModePub = 'METHOD_ALL';
%  end;


  ArrythmiaDetection = sprintf('0x%02d', header(1));

  switch ArrythmiaDetection,
   case '0x01',
    physio.iPmuADPub = 'AD_NONE';
   case '0x02',
    physio.iPmuADPub = 'AD_TIMEBASED';
   case '0x04',
    physio.iPmuADPub = 'AD_PATTERNBASED';
  end;


  PhysioSignal = sprintf('0x%02d', header(2));

  switch PhysioSignal,
   case '0x01',
    physio.iPmuHighPrioTriggerSignal = 'SIGNAL_NONE';
   case '0x02',
    physio.iPmuHighPrioTriggerSignal = 'SIGNAL_EKG';
   case '0x04',
    physio.iPmuHighPrioTriggerSignal = 'SIGNAL_PULSE';
   case '0x08',
    physio.iPmuHighPrioTriggerSignal = 'SIGNAL_EXT';
   case '0x0E',
    physio.iPmuHighPrioTriggerSignal = 'SIGNAL_CARDIAC';
   case '0x10',
    physio.iPmuHighPrioTriggerSignal = 'SIGNAL_RESPIRATION';
   case '0x1E',
    physio.iPmuHighPrioTriggerSignal = 'SIGNAL_ALL';
  end;


  %==--------------------------------------------------------------------==%

  fseek(fp, fpos1, 'bof');
  signal = fscanf(fp, '%d', inf);
  fpos2 = ftell(fp);

  signal(end) = [];


  puls = signal;

  signal_mask = ones(size(signal));

  trigger_1_ind = find(signal == PHYSIO_TRIGGER_ON);
  trigger_0_ind = find(signal == PHYSIO_TRIGGER_OFF);

  puls(trigger_1_ind) = [];
  puls(trigger_0_ind) = [];

  signal_mask(trigger_1_ind) = 0;
  signal_mask(trigger_0_ind) = 0;

  signal_index = cumsum(signal_mask);

  physio.trigger_on  = zeros(size(puls));
  physio.trigger_off = zeros(size(puls));

  physio.trigger_on( signal_index(trigger_1_ind)) = 1;
  physio.trigger_off(signal_index(trigger_0_ind)) = 1;


  % PULS appears to be sampled at 50 Hz
  t = (1:length(puls))' * 2.0e-2;


  %==--------------------------------------------------------------------==%

  fseek(fp, fpos2, 'bof');
  footer = sscanf(fgetl(fp), 'ECG  Freq Per: %d %d');
  physio.ECG.Freq = footer(1);
  physio.ECG.Per  = footer(2);

  footer = sscanf(fgetl(fp), 'PULS Freq Per: %d %d');
  physio.PULS.Freq = footer(1);
  physio.PULS.Per  = footer(2);

  footer = sscanf(fgetl(fp), 'RESP Freq Per: %d %d');
  physio.RESP.Freq = footer(1);
  physio.RESP.Per  = footer(2);

  footer = sscanf(fgetl(fp), 'EXT  Freq Per: %d %d');
  physio.EXT.Freq = footer(1);
  physio.EXT.Per  = footer(2);

  footer = sscanf(fgetl(fp), 'ECG  Min Max Avg StdDiff: %d %d %d %d');
  physio.ECG.Min     = footer(1);
  physio.ECG.Max     = footer(2);
  physio.ECG.Avg     = footer(3);
  physio.ECG.StdDiff = footer(4);

  footer = sscanf(fgetl(fp), 'PULS Min Max Avg StdDiff: %d %d %d %d');
  physio.PULS.Min     = footer(1);
  physio.PULS.Max     = footer(2);
  physio.PULS.Avg     = footer(3);
  physio.PULS.StdDiff = footer(4);

  footer = sscanf(fgetl(fp), 'RESP Min Max Avg StdDiff: %d %d %d %d');
  physio.RESP.Min     = footer(1);
  physio.RESP.Max     = footer(2);
  physio.RESP.Avg     = footer(3);
  physio.RESP.StdDiff = footer(4);

  footer = sscanf(fgetl(fp), 'EXT  Min Max Avg StdDiff: %d %d %d %d');
  physio.EXT.Min     = footer(1);
  physio.EXT.Max     = footer(2);
  physio.EXT.Avg     = footer(3);
  physio.EXT.StdDiff = footer(4);


  footer = sscanf(fgetl(fp), 'NrTrig NrMP NrArr AcqWin: %d %d %d %d');
  physio.NrTrig = footer(1);
  physio.NrMP   = footer(2);
  physio.NrArr  = footer(3);
  physio.AcqWin = footer(4);

  physio.LogStartMDHTime  = sscanf(fgetl(fp), 'LogStartMDHTime:  %d');
  physio.LogStopMDHTime   = sscanf(fgetl(fp), 'LogStopMDHTime:   %d');
  physio.LogStartMPCUTime = sscanf(fgetl(fp), 'LogStartMPCUTime: %d');
  physio.LogStopMPCUTime  = sscanf(fgetl(fp), 'LogStopMPCUTime:  %d');


  %==--------------------------------------------------------------------==%

  eof = fscanf(fp, '%d', 1);

  if ( eof ~= PHYSIO_EOF ),
    error('wrong EOF character encountered');
  end;


  fseek(fp, fpos2, 'bof');
  fclose(fp);


  if ( nargout == 0 ),

    figure('name', mfilename); axis; hold on; box on;
    stem(t, physio.trigger_on*max(puls),  'g', 'Marker', 'none');
    stem(t, physio.trigger_off*max(puls), 'r', 'Marker', 'none');
    plot(t, puls, 'k.-');
    xlabel('time (s)');
    puls_interval = mean(diff(t(find(physio.trigger_on))));
    title(sprintf('PF: %d / min, PP: %d ms', ...
                  round(puls_interval*60), round(puls_interval*1000)));



  end;



  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_physio_puls.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:



%%%
%%% Some comments on ECG data files (subject to change)
%%% The stored data is a string of numbers that looks like:
%%% 8 1 2 40 280 2292 2300 2306 ....... followed by a tail of average values etc..
%%% The 1st NUMBER encodes the method
%%% 1. lPmuECGModePub:
%%%   enum PhysioMethod
%%%   {
%%%     METHOD_NONE = 0x01,
%%%     METHOD_TRIGGERING = 0x02,
%%%     METHOD_GATING = 0x04,
%%%     METHOD_RETROGATING = 0x08,
%%%     METHOD_SOPE = 0x10,
%%%     METHOD_ALL = 0x1E
%%%   };
%%% The 2nd NUMBER encodes the ArrhythmiaDetection
%%% 2. iPmuADPub
%%%   enum ArrhythmiaDetection
%%%   {
%%%     AD_NONE = 0x01,
%%%     AD_TIMEBASED = 0x02,
%%%     AD_PATTERNBASED = 0x04
%%%   };
%%% The 3rd NUMBER encodes the signal used
%%% 3. iPmuHighPrioTriggerSignal (source of beep)
%%%   enum PhysioSignal
%%%   {
%%%     SIGNAL_NONE = 0x01,
%%%     SIGNAL_EKG = 0x02,
%%%     SIGNAL_PULSE = 0x04,
%%%     SIGNAL_EXT = 0x08,
%%%     SIGNAL_CARDIAC = 0x0E, /* the sequence usually takes this */
%%%     SIGNAL_RESPIRATION = 0x10,
%%%     SIGNAL_ALL = 0x1E,
%%%   };
%%% The 4th & 5th NUMBER encode gate open and close times in tick-time unit.
%%% current tick time is 2.5 ms (see example values above)
%%% // 4. ulECGGateOnCountPub (i.e. gate opens 100 ms after R wave)
%%% // 5. ulECGGateOffCountPub (i.e. gate closes 700 ms after R wave)
%%% All following numbers are signal values as function of sampling interval. The sample rate for the
%%% ECG/Cardiac and external signal is 400 Hz. The special value 5000 is used to mark a trigger on signal. The
%%% value 6000 is a trigger off mark.
%%%
