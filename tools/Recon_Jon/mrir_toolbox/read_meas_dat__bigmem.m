function [sos, mdh] = read_meas_dat__bigmem(filename, options)
%READ_MEAS_DAT  read in Siemens format raw data, VB13A/VB15A-style "meas.dat"
%
% data = read_meas_dat(filename, <options>)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/mar/17
% $Id: read_meas_dat.m,v 1.10 2008/03/04 04:10:31 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.10 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %------------------------------------------------------------------------%

  [fp, errstr] = fopen(filename, 'r', 'l');
  if ( fp == -1 ),
    error(errstr);
  end;

  % find beginning of binary data
  data_start = fread(fp, 1, 'uint32');


  %------------------------------------------------------------------------%

  % read header into one string for parsing
  header = fscanf(fp, '%c', data_start-4);



  param_list = {'NColMeas', 'NLinMeas', 'NChaMeas', 'NSetMeas', 'NEcoMeas', ...
                'NPhsMeas', 'NRepMeas', 'NSegMeas', 'NParMeas', 'NSlcMeas', ...
                'NIdaMeas', 'NIdbMeas', 'NIdcMeas', 'NIddMeas', 'NIdeMeas', ...
                'NAveMeas'};

  dimensions = cell2struct(cell(length(param_list),1), param_list, 1);
  dim = [];

  % scan through header for each of the ICE dimension values
  for ind = 1:length(param_list),
    param = param_list{ind};

    % exploit MATLAB regexp machinery to pull out parameter/value pairs
    match = regexp(header, ['(?<param>' param, ').{0,5}\{\s*(?<value>\d*)\s*\}'], 'names');

    % check if no match is found
    if ( isempty(match) ),
      if ( IS__RELEASE_VERSION ),
        IS__RELEASE_VERSION = 0;
        warning('SIEMENS:IO:versioning', 'missing header data---check ICE version');
      end;
      continue;
    end;

    % consider only last match (there can be as many as three in Config_.evp)
    match = match(end);

    % empty means number of elements in this dimension = 1
    if ( isempty(match.value) ),
      match.value = '1';
    end;

    % save out struct and numerical array
    dim(ind) = str2double(match.value);
    dimensions.(param_list{ind}) = dim(ind);

  end;

  % jump past ascii header to binary data
  fseek(fp, data_start, 'bof');

  % read in first MDH
  mdh = read_meas_dat__mdh(fp);


  fclose(fp);


  %------------------------------------------------------------------------%

  dimensions.NColMeas = double(mdh.ushSamplesInScan)

  mdh_length_float32 = 32;

  %  dims = [dimensions.NChaMeas, (dimensions.NColMeas*2 + mdh_length_float32), dimensions.NParMeas, dimensions.NLinMeas];
  dims = [(dimensions.NColMeas*2 + mdh_length_float32), dimensions.NChaMeas, dimensions.NParMeas];


  mmap = memmapfile(filename, 'format', {'single' dims 'ice'}, 'repeat', 1, 'writable', false, 'offset', data_start);


  disp('generated memory map!');


  sos = zeros(dimensions.NColMeas/2, dimensions.NLinMeas, dimensions.NParMeas, 'single');

  for iCha = 19:dimensions.NChaMeas,

    disp(sprintf('channel %02d', iCha));

    chan = zeros(dims(1), dimensions.NLinMeas, dimensions.NParMeas, 'single');

    tic;
    for iLin = 1:dimensions.NLinMeas,

      fprintf('%3d...', iLin);

      mmap.offset = data_start + prod(dims)*4*(iLin-1);

      chan(:,iLin,:) = mmap.Data.ice(:, iCha, :);

    end;
    toc;

    chan_cplx = complex(chan(33:2:end,:,:), chan(34:2:end,:,:));


    tic; img = mrir_iDFT(chan_cplx, 1); toc;
    tic; img = mrir_iDFT(img, 2); toc;
    tic; img = mrir_image_crop(mrir_iDFT(img, 3)); toc;

    sos = sos + abs(img).^2;


    disp('splitting up partitions...');
    tic;
    uncomb1(:,:,:) = img(:,:,  1:128);
    uncomb2(:,:,:) = img(:,:,129:256);
    uncomb3(:,:,:) = img(:,:,257:384);
    uncomb4(:,:,:) = img(:,:,385:512);
    uncomb5(:,:,:) = img(:,:,513:640);
    uncomb6(:,:,:) = img(:,:,641:768);
    uncomb7(:,:,:) = img(:,:,769:896);
    uncomb8(:,:,:) = img(:,:,897:end);
    toc;


    disp('saving...');
    tic;
    eval(sprintf('save uncomb_%02d uncomb* iLin iCha', iCha));
    toc;

  end;

  jnotify;


  return;



  ACQEND_length_float32 = 256/4;



  %************************************************************************%
  %%% $Source: /space/repo/1/dev/dev/matlab/read_meas_dat.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
