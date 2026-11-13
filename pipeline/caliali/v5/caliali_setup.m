function caliali_setup()
% CALIALI_SETUP  Put CaliAli (and, on old MATLAB without clim, the clim shim) on the path.
here = fileparts(mfilename('fullpath'));
if ~exist('CaliAli_parameters', 'file')
    addpath(genpath('/mnt/data1/CaliAli/CaliAli-main'));      % server install (patched in place)
end
if verLessThan('matlab', '9.12') && ~exist('clim', 'file') ...  % clim exists from R2022a (9.12)
        && isfolder(fullfile(here, 'compat_R2020b'))
    addpath(fullfile(here, 'compat_R2020b'));
end
end
