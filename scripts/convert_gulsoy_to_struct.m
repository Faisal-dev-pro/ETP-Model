function convert_gulsoy_to_struct(inMat, outMat, varargin)
% CONVERT_GULSOY_TO_STRUCT  One-shot converter from the Gulsoy
% TR_dataTable.mat (MATLAB table format with opaque objects, which
% scipy.io.loadmat cannot read) to a plain struct-array .mat that any
% downstream tool can consume.
%
% Also populates the .heater_power field (not present in source data)
% via configure_gulsoy_heaters with default settings.
%
% Run this ONCE after downloading the dataset:
%
%   convert_gulsoy_to_struct();                                     % defaults
%   convert_gulsoy_to_struct('data/TR_dataTable.mat', ...
%                            'data/gulsoy_parsed.mat')
%   convert_gulsoy_to_struct(in, out, 'wattage', 80)                % override
%
% Produces:
%   data/gulsoy_parsed.mat  containing a single variable 'data', a struct
%   array of length N (one entry per test). See load_gulsoy_data.m for
%   field documentation.
%
% Use -v7 (not -v7.3) so the output remains scipy-readable.
%
% Author: <your name>, 2026.  License: MIT.

    if nargin < 1 || isempty(inMat),  inMat  = fullfile('data','TR_dataTable.mat'); end
    if nargin < 2 || isempty(outMat), outMat = fullfile('data','gulsoy_parsed.mat'); end

    fprintf('convert_gulsoy_to_struct\n  in : %s\n  out: %s\n', inMat, outMat);

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    data = load_gulsoy_data(inMat);

    % Populate heater profiles (forwards any extra varargin)
    data = configure_gulsoy_heaters(data, varargin{:});

    save(outMat, 'data', '-v7');
    fprintf('\nSaved %d test record(s) to %s (%.1f MB)\n', ...
            numel(data), outMat, dir_size_mb(outMat));
end

function mb = dir_size_mb(p)
    d = dir(p);
    if isempty(d), mb = 0; else, mb = d(1).bytes/1024/1024; end
end
