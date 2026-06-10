function results = evaluate_phase3_fit(varargin)
% EVALUATE_PHASE3_FIT  Evaluate the Phase 3 pressure sub-model against
% all three Gulsoy test cells.
%
% Computes:
%   - Pre-vent pressure RMSE (target < 2 bar per plan section 4.1)
%   - Soft-vent activation time error
%   - Also re-reports temperature RMSE for completeness
%
% Usage:
%   results = evaluate_phase3_fit();
%   results = evaluate_phase3_fit('paramsFile','params/etp_params_nmc_p3.mat');
%
% Author: <your name>, 2026. License: MIT.

    opts.paramsFile   = fullfile('params','etp_params_nmc_p3.mat');
    opts.dataFile     = fullfile('data','gulsoy_parsed.mat');
    opts.P_burst_clip = 25;
    opts.savePlots    = true;
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'params'), fullfile(root,'scripts'));

    assert(isfile(opts.paramsFile), 'Missing %s - run Phase 3 fit first', opts.paramsFile);
    assert(isfile(opts.dataFile),   'Missing %s', opts.dataFile);

    p = load(opts.paramsFile).p_fitted;
    data = load(opts.dataFile).data;
    N = numel(data);

    fprintf('\n=== Phase 3 pressure evaluation ===\n');
    fprintf('  P_burst clip: %.0f bar\n', opts.P_burst_clip);
    fprintf('  Cells: %d\n\n', N);

    template = struct('cell_id','', 'rmse_P_prevent',NaN, ...
        'rmse_T_prevent',NaN, 'rmse_T_full',NaN, ...
        't_vent_err',NaN, 'P_peak_meas',NaN, 'P_peak_pred',NaN);
    results = repmat(template, N, 1);

    fprintf('%-10s | %-13s | %-13s | %-11s | %-11s | %-12s\n', ...
        'cell_id','RMSE_P_pvent','RMSE_T_pvent','P_pk_meas','P_pk_pred','t_vent_err');
    fprintf('%s\n', repmat('-',1,85));

    for i = 1:N
        d = data(i);
        r = score_cell(p, d, opts);
        results(i).cell_id        = d.cell_id;
        results(i).rmse_P_prevent = r.rmse_P;
        results(i).rmse_T_prevent = r.rmse_T;
        results(i).rmse_T_full    = r.rmse_T_full;
        results(i).t_vent_err     = r.t_vent_err;
        results(i).P_peak_meas    = r.P_peak_meas;
        results(i).P_peak_pred    = r.P_peak_pred;

        fprintf('%-10s | %5.3f bar    | %5.1f K       | %5.1f bar  | %5.1f bar  | %+6.1f s\n', ...
            d.cell_id, r.rmse_P, r.rmse_T, ...
            r.P_peak_meas, r.P_peak_pred, r.t_vent_err);

        if opts.savePlots
            plot_cell(r, d, fullfile(root,'results', ...
                sprintf('phase3_eval_cell%d.png', i)));
        end
    end

    fprintf('\nPaper-target check (Section 5.2, pre-vent P RMSE < 2 bar):\n');
    for i = 1:N
        if results(i).rmse_P_prevent <= 2
            status = '[PASS]';
        else
            status = '[over target]';
        end
        fprintf('  %-10s : %.3f bar  %s\n', ...
            data(i).cell_id, results(i).rmse_P_prevent, status);
    end
end

% =====================================================================
function r = score_cell(p, d, opts)
    p.pressure.enable = true;
    p.cell.T0 = d.T_internal(1) + 273.15;
    if ~isempty(d.heater_power)
        p.Qext_profile.time  = d.t;
        p.Qext_profile.value = d.heater_power;
    end
    p.solver.tFinal = d.t(end);
    out = etp_pure_matlab(p, 'verbose', false);

    % Predicted P in bar on the slow time grid
    P_pred_bar = interp1(out.t, out.P/1e5, d.t, 'linear', NaN);
    T_pred_C   = interp1(out.t, out.T-273.15, d.t, 'linear', NaN);

    % Measured P on the slow grid (interpolated from fast)
    P_meas_bar = interp1(d.t_fast, d.P_internal, d.t, 'linear', NaN);

    % Pre-vent mask: before either measured or predicted crosses clip
    iVent_meas = find(P_meas_bar > opts.P_burst_clip, 1, 'first');
    iVent_pred = find(P_pred_bar > opts.P_burst_clip, 1, 'first');
    if isempty(iVent_meas), [~, iVent_meas] = max(P_meas_bar); end  % fallback: peak P
    if isempty(iVent_pred), [~, iVent_pred] = max(P_pred_bar); end  % fallback: peak P
    iEnd = min(iVent_meas, iVent_pred);

    pv = (1:numel(P_meas_bar))' <= iEnd;
    valid = pv & ~isnan(P_meas_bar) & ~isnan(P_pred_bar);

    err_P = P_pred_bar - P_meas_bar;
    if any(valid)
        r.rmse_P = sqrt(mean(err_P(valid).^2));
    else
        r.rmse_P = NaN;
    end

    % Temperature pre-vent RMSE (same window as Phase 2 eval)
    err_T = T_pred_C - d.T_internal;
    T_pv = valid & (d.T_internal < 200) & (T_pred_C < 200);
    if any(T_pv)
        r.rmse_T = sqrt(mean(err_T(T_pv).^2));
    else
        r.rmse_T = NaN;
    end
    mask_all = ~isnan(err_T);
    r.rmse_T_full = sqrt(mean(err_T(mask_all).^2));

    % Peak pressure
    r.P_peak_meas = max(P_meas_bar);
    r.P_peak_pred = max(P_pred_bar);

    % Vent timing: first time P crosses 15 bar
    vent_thresh = 15;
    iV_meas = find(P_meas_bar > vent_thresh, 1, 'first');
    iV_pred = find(P_pred_bar > vent_thresh, 1, 'first');
    if ~isempty(iV_meas) && ~isempty(iV_pred)
        r.t_vent_err = d.t(iV_pred) - d.t(iV_meas);
    else
        r.t_vent_err = NaN;
    end

    % Stash for plotting
    r.P_pred = P_pred_bar;
    r.P_meas = P_meas_bar;
    r.T_pred = T_pred_C;
    r.t = d.t;
    r.iEnd = iEnd;
end

function plot_cell(r, d, savePath)
    fig = figure('Color','w','Position',[100 100 1000 700],'Visible','off');
    tlo = tiledlayout(fig, 3, 1, 'TileSpacing','compact','Padding','compact');
    title(tlo, sprintf('Phase 3 eval: cell %s  (P RMSE_{pvent}=%.3f bar)', ...
        d.cell_id, r.rmse_P), 'FontWeight','bold','Interpreter','none');

    % Panel 1: pressure
    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
    plot(ax1, r.t/60, r.P_meas, 'k-', 'LineWidth',1.5, 'DisplayName','measured');
    plot(ax1, r.t/60, r.P_pred, 'r--','LineWidth',1.5, 'DisplayName','model');
    ylim(ax1, [min(r.P_meas)-5, max(r.P_meas(1:r.iEnd))*1.5]);  % clip to pre-vent range
    if r.iEnd <= numel(r.t)
        xline(ax1, r.t(r.iEnd)/60, 'b:', 'pre-vent end');
    end
    ylabel(ax1,'P [bar]'); legend(ax1,'Location','northwest');
    title(ax1,'Internal gas pressure');

    % Panel 2: pressure residual (pre-vent zoom)
    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
    resid = r.P_pred - r.P_meas;
    plot(ax2, r.t/60, resid, 'b-', 'LineWidth',1);
    yline(ax2, 0,'k-'); yline(ax2, [-2 2], 'b:',{'-2 bar','+2 bar'});
    ylabel(ax2,'P_{pred}-P_{meas} [bar]'); title(ax2,'Pressure residual');
    if r.iEnd <= numel(r.t)
        xlim(ax2, [0 r.t(min(r.iEnd+100,numel(r.t)))/60]);
    end

    % Panel 3: temperature (confirmation that T fit still holds)
    ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on');
    plot(ax3, r.t/60, d.T_internal, 'k-', 'LineWidth',1.5, 'DisplayName','T meas');
    plot(ax3, r.t/60, r.T_pred,     'r--','LineWidth',1.2, 'DisplayName','T model');
    ylabel(ax3,'T [degC]'); xlabel(ax3,'t [min]');
    legend(ax3,'Location','northwest'); title(ax3,'Temperature (confirmation)');

    linkaxes([ax1 ax2 ax3], 'x');

    outDir = fileparts(savePath);
    if ~isempty(outDir) && ~exist(outDir,'dir'), mkdir(outDir); end
    exportgraphics(fig, savePath, 'Resolution', 130);
    close(fig);
    fprintf('  Saved plot: %s\n', savePath);
end
