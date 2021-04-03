function [rss, mdh] = read_meas_test(filename, options)
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
  
  mmap = memmapfile(filename, 'format', 'single', 'writable', false, 'offset', data_start);
  disp('generated memory map!');

  mdh_length_float32 = 32;

  float32_per_lin = (dimensions.NColMeas*2 + mdh_length_float32) * dimensions.NChaMeas * dimensions.NParMeas;

  fprintf('allocating memory...');
  tic; data = zeros(float32_per_lin, dimensions.NLinMeas, 'single');
  fprintf('done! '); toc;

  fprintf('importing data......');
  tic;
  for iLin = 1:dimensions.NLinMeas,
    data(:, iLin) = mmap.data([float32_per_lin*(iLin-1)+1]:[float32_per_lin*(iLin-0)]);
  end;
  fprintf('done! '); toc;

  fprintf('repacking data......');
  data = reshape(data.', dimensions.NLinMeas, (dimensions.NColMeas*2 + mdh_length_float32), dimensions.NChaMeas, 1, 1, 1, 1, 1, dimensions.NParMeas);
  fprintf('done! '); toc;
  
  time data = permute(data, [2, 1, 3, 4, 5, 6, 7, 8, 9]);

  
  %                                         1 2 3 4 5 6 7 8 9
  time mdh_binary = data(1:mdh_length_float32,:,:,:,:,:,:,:,:);
  data(1:mdh_length_float32,:,:,:,:,:,:,:,:) = [];


  fprintf('calculating reconstruction...');
  tic;

  sos = single(0);

  for iCha = 1:dimensions.NChaMeas,

    raw = complex(data(1:2:end, :, iCha, 1, 1, 1, 1, 1, :), data(2:2:end, :, iCha, 1, 1, 1, 1, 1, :));
    img = mrir_conventional_3d(raw);

    sos = sos + abs(img).^2;

  end;
  fprintf('done! '); toc;


  rss = sqrt(sos);


  return;


  %************************************************************************%
  %%% $Source: /space/repo/1/dev/dev/matlab/read_meas_dat.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
