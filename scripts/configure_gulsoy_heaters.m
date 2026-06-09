function data = configure_gulsoy_heaters(data, varargin)
% CONFIGURE_GULSOY_HEATERS  Populate the .heater_power field for each
% test in the parsed Gulsoy dataset.
%
% The TR_dataTable export does not include the heater wattage trace.
% The Gulsoy 2025 Methods section specifies the protocol as a constant
% heater output until thermal runaway is detected, then heater off.
%
% Default protocol (override as needed for your specific experiment):
%   - Constant heater output for the entire recording.
%     The Gulsoy 2024 methodology paper (J. Power Sources 617, 235147)
%     specifies the heater wattage; default here is 100 W as a placeholder.
%
% Usage:
%   data = configure_gulsoy_heaters(data);
%   data = configure_gulsoy_heaters(data, 'wattage', 80);
%   data = configure_gulsoy_heaters(data, 'wattage', [100 100 80]);  % per-test
%   data = configure_gulsoy_heaters(data, 'mode', 'cutoff_at_runaway');
%
% Options:
%   wattage         scalar or N-vector       (default 100 W)
%   mode            'always_on' (default) or 'cutoff_at_runaway'
%   off_delay       seconds after onset before cutoff (default 10)
%   skip_seconds    ignore the first N s when finding onset (default 30)
%   smooth_window   samples for dT/dt smoothing (default 11)
%
% After this, every data(i).heater_power is a vector aligned with data(i).t.
%
% Author: <your name>, 2026. License: MIT.

    opts.wattage       = 100;
    opts.mode          = 'always_on';
    opts.off_delay     = 10;
    opts.skip_seconds  = 30;
    opts.smooth_window = 11;
    opts.verbose       = true;
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end

    N = numel(data);
    if isscalar(opts.wattage), W = repmat(opts.wattage, N, 1);
    else, assert(numel(opts.wattage) == N, ...
            'wattage must be scalar or have %d entries', N);
         W = opts.wattage(:); end

    if opts.verbose
        fprintf('configure_gulsoy_heaters: mode=%s\n', opts.mode);
    end

    for i = 1:N
        t = data(i).t;
        if isempty(t)
            warning('Test %d has empty time vector; skipping', i);
            continue;
        end
        q = W(i) * ones(size(t));

        if strcmp(opts.mode, 'cutoff_at_runaway')
            % Robust onset detection: skip the first opts.skip_seconds
            % to avoid gradient() boundary artefacts, and smooth dT/dt.
            T = data(i).T_internal;
            if ~isempty(T) && numel(T) == numel(t)
                dT = gradient(T, t);
                if opts.smooth_window > 1
                    dT = movmean(dT, opts.smooth_window);
                end
                idx_start = find(t >= opts.skip_seconds, 1, 'first');
                if isempty(idx_start), idx_start = 1; end
                iOn_rel = find(dT(idx_start:end) > 1, 1, 'first');
                if ~isempty(iOn_rel)
                    iOn = iOn_rel + idx_start - 1;
                    t_off = t(iOn) + opts.off_delay;
                    q(t > t_off) = 0;
                    if opts.verbose
                        fprintf('  [%d] %s: %.0f W until t=%.1fs (onset detected at %.1fs)\n', ...
                                i, data(i).cell_id, W(i), t_off, t(iOn));
                    end
                else
                    if opts.verbose
                        fprintf('  [%d] %s: %.0f W (no onset detected, full duration)\n', ...
                                i, data(i).cell_id, W(i));
                    end
                end
            end
        else
            if opts.verbose
                fprintf('  [%d] %s: %.0f W for full duration (%.0f s)\n', ...
                        i, data(i).cell_id, W(i), t(end));
            end
        end

        data(i).heater_power = q;
    end
end
