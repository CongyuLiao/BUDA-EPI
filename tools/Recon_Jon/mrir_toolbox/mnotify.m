function result = mnotify(email, subject, message)
% MNOTIFY  send email to user alerting that a job is completed
%
%  result = mnotify(email, subject, message)

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %------------------------------------------------------------------------%

  subject = sprintf('[MATLAB] %s', subject);
  message = sprintf('%s\\n\\n%s', message, datestr(now, 0));

  [status, result] = system(sprintf('echo -e "%s" | mail -v -s "%s" %s', message, ...
                                    subject, email));

  if ( status ~= 0 ),
    errstr = sprintf('email notification to "%s" failed', email);
    error(sprintf('==> [%s]: %s\n', mfilename, errstr));
  end;

  return;
