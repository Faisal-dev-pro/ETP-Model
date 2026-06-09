function data = load_gulsoy_data(matPath)
% LOAD_GULSOY_DATA  Parse the Gulsoy et al. 2025 TR_dataTable.mat into
% a clean struct-array indexed by test cell.
%
% Source layout (verified 2026-05-26):
%   Wide table, 3 rows = 3 tests, 11 columns.
%   Two distinct sample rates per test:
%     - High-rate (~9.56M samples): ExpTime, IntPre, CellVoltage
%     - Low-rate (~9558 samples)  : ExpTimeTemp, MidIntTemp, MidSurfTemp,
%                                   NegSurfTemp, PosSurfTemp,
%                                   VentPos5mmAway, VentPos10mmAway
%   Plus TestID (categorical, 3x1) identifying each row.
%
% Returns:
%   data : (N x 1) struct array, one entry per test cell, with fields
%     .cell_id        char     test identifier (from TestID)
%
%     -- HIGH-RATE TIMESERIES (~10 kHz typical) --
%     .t_fast         [s]      time vector for pressure / voltage
%     .P_internal     [bar a]  internal gas pressure
%     .V_cell         [V]      cell voltage
%
%     -- LOW-RATE TIMESERIES (~10 Hz typical) --
%     .t              [s]      time vector for temperature signals
%     .T_internal     [degC]   internal axial-midpoint temperature
%     .T_surface_mid  [degC]   surface midpoint
%     .T_surface_pos  [degC]   surface near positive terminal
%     .T_surface_neg  [degC]   surface near negative terminal
%     .T_vent_5mm     [degC]   vent gas at 5 mm from vent
%     .T_vent_10mm    [degC]   vent gas at 10 mm from vent
%
%     -- METADATA --
%     .fs_fast        [Hz]     estimated high-rate sampling frequency
%     .fs_slow        [Hz]     estimated low-rate sampling frequency
%     .heater_power   [W]      EMPTY in source data; configure per test
%                              (Gulsoy paper Methods section)
%
% Note: the development plan asks for a single time vector 'data(i).t'.
% This loader assigns the temperature time grid to '.t' because the
% Phase-2 optimiser fits T_internal (slow-rate), and adds '.t_fast' for
% the pressure data Phase-3 will fit against.  Resampling between the
% two should be done by the consumer (Phase 3) so the loader stays
% lossless.
%
% Author: <your name>, 2026. License: MIT.

    if nargin < 1 || isempty(matPath)
        matPath = fullfile('data','TR_dataTable.mat');
    end
    assert(isfile(matPath), 'File not found: %s', matPath);

    S = load(matPath);

    % find the TR_dataTable
    tbl = []; tbl_name = '';
    fns = fieldnames(S);
    for k = 1:numel(fns)
        if istable(S.(fns{k}))
            tbl = S.(fns{k}); tbl_name = fns{k}; break;
        end
    end
    assert(~isempty(tbl), 'No table variable found in %s', matPath);
    fprintf('load_gulsoy_data: loaded table "%s" (%d rows x %d vars)\n', ...
            tbl_name, height(tbl), width(tbl));

    % Column mapping (verified against the actual Gulsoy 2025 export)
    %   canonical_field  <-  source column name
    COLMAP = { ...
        'cell_id',         'TestID'         ; ...
        % high-rate (pressure / voltage stream)
        't_fast',          'ExpTime'        ; ...
        'P_internal',      'IntPre'         ; ...
        'V_cell',          'CellVoltage'    ; ...
        % low-rate (temperature stream)
        't',               'ExpTimeTemp'    ; ...
        'T_internal',      'MidIntTemp'     ; ...
        'T_surface_mid',   'MidSurfTemp'    ; ...
        'T_surface_neg',   'NegSurfTemp'    ; ...
        'T_surface_pos',   'PosSurfTemp'    ; ...
        'T_vent_5mm',      'VentPos5mmAway' ; ...
        'T_vent_10mm',     'VentPos10mmAway'};

    varNames = tbl.Properties.VariableNames;
    N = height(tbl);
    data = repmat(empty_record(), N, 1);

    % Verify every expected column exists
    missing = {};
    for k = 1:size(COLMAP,1)
        if ~ismember(COLMAP{k,2}, varNames)
            missing{end+1} = COLMAP{k,2}; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        warning('Missing expected columns: %s', strjoin(missing, ', '));
    end

    % Extract per-test
    for i = 1:N
        for k = 1:size(COLMAP,1)
            fld = COLMAP{k,1};
            src = COLMAP{k,2};
            if ~ismember(src, varNames), continue; end
            val = tbl{i, src};
            if iscell(val), val = val{1}; end
            % TestID is categorical -> convert to char
            if iscategorical(val)
                data(i).(fld) = char(val);
            elseif isstring(val)
                data(i).(fld) = char(val);
            elseif isnumeric(val)
                data(i).(fld) = double(val(:));
            else
                data(i).(fld) = val;
            end
        end
        if isempty(data(i).cell_id)
            data(i).cell_id = sprintf('cell_%02d', i);
        end

        % derived sample-rate metadata
        data(i).fs_slow = estimate_fs(data(i).t);
        data(i).fs_fast = estimate_fs(data(i).t_fast);
        data(i).heater_power = [];  % not in source; populate per-test
    end

    % --- report ---------------------------------------------------------
    fprintf('\nload_gulsoy_data: parsed %d test(s):\n\n', N);
    for i = 1:N
        report_one(data(i), i);
    end
    fprintf(['\nNOTE: heater_power is empty (not in source data). The Gulsoy\n' ...
             '      methods section specifies the heater profile; set it per\n' ...
             '      test before running Phase 2, e.g.:\n' ...
             '        data(1).heater_power = constant_W * ones(size(data(1).t));\n']);
end

% =====================================================================
function fs = estimate_fs(t)
    if isempty(t) || numel(t) < 2
        fs = NaN; return;
    end
    dt = median(diff(t(1:min(end,1000))));
    if dt <= 0, fs = NaN; else, fs = 1/dt; end
end

function rec = empty_record()
    rec = struct( ...
        'cell_id', '', ...
        't', [], 'T_internal', [], ...
        'T_surface_mid', [], 'T_surface_pos', [], 'T_surface_neg', [], ...
        'T_vent_5mm', [], 'T_vent_10mm', [], ...
        't_fast', [], 'P_internal', [], 'V_cell', [], ...
        'fs_slow', NaN, 'fs_fast', NaN, ...
        'heater_power', []);
end

function report_one(d, i)
    fprintf('  [%d] cell_id = %s\n', i, d.cell_id);
    order = {'t','T_internal','T_surface_mid','T_surface_pos','T_surface_neg', ...
             'T_vent_5mm','T_vent_10mm','t_fast','P_internal','V_cell'};
    for k = 1:numel(order)
        v = d.(order{k});
        if isempty(v)
            fprintf('       %-15s : (empty)\n', order{k});
        else
            fprintf('       %-15s : %8d samples  [min=%9.3g  max=%9.3g]\n', ...
                order{k}, numel(v), min(v), max(v));
        end
    end
    if ~isnan(d.fs_slow)
        fprintf('       fs_slow         : %.2f Hz\n', d.fs_slow);
    end
    if ~isnan(d.fs_fast)
        fprintf('       fs_fast         : %.2f Hz\n', d.fs_fast);
    end
    fprintf('\n');
end
