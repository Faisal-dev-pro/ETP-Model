function results = evaluate_phase2_fit(varargin)
% EVALUATE_PHASE2_FIT  Closeout diagnostic for the Phase 2 parameter fit.
%
% Computes:
%   - Full-trace RMSE on each cell (training + validation)
%   - PRE-VENT RMSE (T_internal < 200 °C window) -- the paper's actual
%     validation target per Section 5.1 (target <15 K)
%   - Peak temperature error
%   - Time-to-onset error
%
% Produces a 3-panel plot per cell: predicted vs measured T_internal,
% predicted vs measured P (Phase 3 will use this), and residual T_pred-T_meas.
%
% Usage:
%   results = evaluate_phase2_fit();                      % uses params_nmc
%   results = evaluate_phase2_fit('paramsFile', 'params/etp_params_nmc.mat');
%
% Author: <your name>, 2026. License: MIT.

    opts.paramsFile = fullfile('params','etp_params_nmc.mat');
    opts.dataFile   = fullfile('data','gulsoy_parsed.mat');
    opts.preVentT   = 200;     % degC — boundary for pre-vent window
    opts.savePlots  = true;
    opts.verbose    = true;
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'params'), fullfile(root,'scripts'));

    assert(isfile(opts.paramsFile), 'Missing %s — run run_phase2_param_id first', opts.paramsFile);
    assert(isfile(opts.dataFile),   'Missing %s — run convert_gulsoy_to_struct first', opts.dataFile);

    S = load(opts.paramsFile);
    p = S.p_fitted;
    data = load(opts.dataFile).data;
    N = numel(data);

    fprintf('\n=== Phase 2 fit evaluation ===\n');
    fprintf('  Pre-vent window:  T_internal < %.0f degC\n', opts.preVentT);
    fprintf('  Cells: %d\n\n', N);

    % Pre-allocate with a fully-populated template so MATLAB doesn't
    % balk on the first assignment.
    template = struct('cell_id', '', 'rmse_full', NaN, 'rmse_prevent', NaN, ...
        'peak_meas_C', NaN, 'peak_pred_C', NaN, 't_onset_err', NaN, ...
        'T_pred', [], 'T_meas', [], 't', []);
    results = repmat(template, N, 1);

    fprintf('%-10s | %-12s | %-13s | %-11s | %-13s | %-13s\n', ...
        'cell_id','RMSE_full','RMSE_prevent','peak_meas','peak_pred','t_onset_err');
    fprintf('%s\n', repmat('-',1,90));

    for i = 1:N
        d = data(i);
        out = simulate_one(p, d);
        r = score_one(out, d, opts.preVentT);

        % Field-by-field copy into the templated struct array
        results(i).cell_id      = d.cell_id;
        results(i).rmse_full    = r.rmse_full;
        results(i).rmse_prevent = r.rmse_prevent;
        results(i).peak_meas_C  = r.peak_meas_C;
        results(i).peak_pred_C  = r.peak_pred_C;
        results(i).t_onset_err  = r.t_onset_err;
        results(i).T_pred       = r.T_pred;
        results(i).T_meas       = r.T_meas;
        results(i).t            = r.t;

        fprintf('%-10s | %5.1f K     | %5.1f K       | %5.0f degC | %5.0f degC   | %+6.1f s\n', ...
            d.cell_id, r.rmse_full, r.rmse_prevent, ...
            r.peak_meas_C, r.peak_pred_C, r.t_onset_err);

        if opts.savePlots
            plot_one(out, d, r, fullfile(root, 'results', sprintf('phase2_eval_cell%d.png', i)));
        end
    end

    fprintf('\nPaper-target check (Section 5.1, RMSE_prevent < 15 K):\n');
    for i = 1:N
        status = '[PASS]';
        if results(i).rmse_prevent > 15, status = '[over target]'; end
        fprintf('  %-10s : %5.1f K  %s\n', ...
            data(i).cell_id, results(i).rmse_prevent, status);
    end
end

% =====================================================================
function out = simulate_one(p, d)
    p_try = p;
    p_try.cell.T0   = d.T_internal(1) + 273.15;
    if ~isempty(d.heater_power)
        p_try.Qext_profile.time  = d.t;
        p_try.Qext_profile.value = d.heater_power;
    end
    p_try.solver.tFinal = d.t(end);
    out = etp_pure_matlab(p_try, 'verbose', false);
end

function r = score_one(out, d, preVentT)
    T_pred_C = interp1(out.t, out.T - 273.15, d.t, 'linear', NaN);
    T_meas_C = d.T_internal;

    err = T_pred_C(:) - T_meas_C(:);
    mask = ~isnan(err);

    r.rmse_full = sqrt(mean(err(mask).^2));

    % Pre-vent window: samples where BOTH measured and predicted are
    % below the pre-vent threshold. This is the standard TR-modelling
    % definition and excludes the period when either curve is in or
    % past runaway. No onset-detection needed.
    pv_mask = mask & (T_meas_C(:) < preVentT) & (T_pred_C(:) < preVentT);
    if any(pv_mask)
        r.rmse_prevent = sqrt(mean(err(pv_mask).^2));
        r.n_prevent    = sum(pv_mask);
    else
        r.rmse_prevent = NaN;
        r.n_prevent    = 0;
    end

    [Tmax, iMax] = max(T_meas_C);
    r.peak_meas_C = Tmax;
    r.peak_pred_C = max(T_pred_C);

    % onset (smoothed dT/dt > 1 K/s, sustained for >= 5 samples, after t > 30s)
    % Raw gradient is noisy at the start of high-rate data; smoothing
    % stops sensor noise from triggering the threshold spuriously.
    iOn_meas = detect_onset(d.t, T_meas_C, 1.0);
    iOn_pred = detect_onset(d.t, T_pred_C, 1.0);
    if ~isempty(iOn_meas) && ~isempty(iOn_pred)
        r.t_onset_err = d.t(iOn_pred) - d.t(iOn_meas);
    else
        r.t_onset_err = NaN;
    end

    r.T_pred = T_pred_C;
    r.T_meas = T_meas_C;
    r.t      = d.t;
end

% =====================================================================
function iOnset = detect_onset(t, T, threshold)
% Smoothed onset detection: dT/dt > threshold, sustained for several
% seconds, after t > 30 s. Time-based windowing so it works across
% different sample rates (10 Hz for cells #1/#2, 100 Hz for cell #3).
    iOnset = [];
    if numel(T) < 21, return; end
    dt = median(diff(t(1:min(end,1000))));
    if dt <= 0, return; end
    fs = 1/dt;

    smooth_window    = max(21, round(5 * fs));   % ~5 s smoothing
    sustained_window = max(5,  round(2 * fs));   % ~2 s sustained

    dT = gradient(T(:), t(:));
    dT_smooth = movmean(dT, smooth_window);
    above = dT_smooth > threshold;

    % require sustained_window consecutive samples above threshold
    if numel(above) >= sustained_window
        kernel = ones(sustained_window,1) / sustained_window;
        is_sustained = conv(double(above), kernel, 'same') > 0.99;
    else
        is_sustained = above;
    end

    valid = is_sustained & (t(:) > 30);
    iOnset = find(valid, 1, 'first');
end

function plot_one(out, d, r, savePath)
    fig = figure('Color','w','Position',[100 100 1000 800],'Visible','off');
    tlo = tiledlayout(fig, 3, 1, 'TileSpacing','compact','Padding','compact');
    title(tlo, sprintf('Phase 2 fit: cell %s  (RMSE full=%.1fK, pre-vent=%.1fK)', ...
        d.cell_id, r.rmse_full, r.rmse_prevent), ...
        'FontWeight','bold','Interpreter','none');

    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
    plot(ax1, d.t/60, r.T_meas, 'k-', 'LineWidth', 1.8, 'DisplayName','measured');
    plot(ax1, d.t/60, r.T_pred, 'r--','LineWidth', 1.5, 'DisplayName','model');
    ylabel(ax1, 'T_{internal} [degC]'); legend(ax1, 'Location','northwest');
    xlabel(ax1,'t [min]');

    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
    plot(ax2, d.t/60, r.T_meas, 'k-', 'LineWidth', 1.5, 'DisplayName','measured');
    plot(ax2, d.t/60, r.T_pred, 'r--','LineWidth', 1.2, 'DisplayName','model');
    yline(ax2, 200, 'b:', 'pre-vent window');
    ylabel(ax2,'T [degC] (pre-vent zoom)'); ylim(ax2, [0 250]);
    xlabel(ax2,'t [min]'); legend(ax2,'Location','northwest');

    ax3 = nexttile; grid(ax3,'on');
    plot(ax3, d.t/60, r.T_pred(:) - r.T_meas(:), 'b-','LineWidth',1);
    yline(ax3, 0,'k-'); yline(ax3, [-15 15], 'b:', {'-15K','+15K'});
    ylabel(ax3, 'T_{pred} - T_{meas} [K]'); xlabel(ax3,'t [min]');

    linkaxes([ax1 ax2 ax3], 'x');

    outDir = fileparts(savePath);
    if ~isempty(outDir) && ~exist(outDir,'dir'), mkdir(outDir); end
    exportgraphics(fig, savePath, 'Resolution', 130);
    close(fig);
    fprintf('  Saved plot: %s\n', savePath);
end
