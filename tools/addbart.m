
% bart_path = '/homes/9/berkin/Downloads/bart-0.2.06/';
% bart_path = '/cluster/kawin/berkin/Matlab_Code_New/LIBRARY/bart-0.2.06/';
bart_path = [pwd, '/bart-0.2.06/'];

setenv('TOOLBOX_PATH', bart_path)
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));
setenv('PATH', strcat(getenv('TOOLBOX_PATH'), ':', getenv('PATH')));
setenv('LD_LIBRARY_PATH', '');
