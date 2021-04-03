


bart_path='/autofs/cluster/kawin/Congyu/CODE/gSlider_BUDA/tools/bart-0.2.08/';
setenv('TOOLBOX_PATH', bart_path)
if isempty(getenv('TOOLBOX_PATH'))
    error('Environment variable TOOLBOX_PATH is not set /autofs/cluster/kawin/Congyu/CODE/gSlider_BUDA/tools/bart-0.2.08/');
end
% set path to bart/matlab folder
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));
setenv('PATH', strcat(getenv('TOOLBOX_PATH'), ':', getenv('PATH')));
setenv('LD_LIBRARY_PATH', '');


