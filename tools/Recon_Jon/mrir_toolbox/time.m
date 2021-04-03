function time(varargin)
% TIME  display begin time, end time, and calculate duration of a command
%
% example:
%
%   time 'd = randn(10);'

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 10/13/2005
% $Id: time.m,v 1.1 2007/12/18 03:13:15 jonnyreb Exp $
%**************************************************************************%

  eval_str = '';

  for arg = 1:nargin,
    eval_str = cat(2, eval_str, varargin{arg});
    eval_str = cat(2, eval_str, ' ');
  end;

  eval_str = [deblank(eval_str), ';'];


  fprintf(1, '\nt0: %s\n\n', datestr(now));

  if ( isempty(varargin) ),
    return;
  end;

  fprintf(1, '  %s\n', eval_str);


  t0 = clock;
  evalin('caller', sprintf('%s', eval_str));
  t1 = clock;

  fprintf(1, '\nt1: %s\n', datestr(now));

  runtime_seconds = etime(t1, t0);

  runtime_str = sprintf('[[ %02dh %02dm %02ds ]]', fix(runtime_seconds/60/60), ...
                        rem(fix(runtime_seconds/60), 60), ...
                        rem(fix(runtime_seconds), 60));


  fprintf(1, '\n\nruntime: %s\n\n', runtime_str);

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/jonnyreb/MATLAB/time.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
