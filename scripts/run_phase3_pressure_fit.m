function results = run_phase3_pressure_fit(varargin)
% RUN_PHASE3_PRESSURE_FIT  Phase 3 of the development plan:
%   Fit the gas-generation sub-model coefficients against the pre-vent
%   internal pressure trace of the training cell.
%
% Fits four parameters:
%   nu_vap  - mol gas per mol electrolyte vapourised
%   nu_e    - mol gas per mol electrolyte decomposed
%   nu_O2   - mol O2 per mol cathode active material
%   V_int   - cell internal free volume [m^3]
%
% Prerequisites:
%   - data/gulsoy_parsed.mat
%   - params/etp_params_nmc.mat (Phase 2 output)
%
% Usage:
%   results = run_phase3_pressure_fit();
%   results = run_phase3_pressure_fit('trainingCell', 2);
%
% Author: <your name>, 2026. License: MIT.

    opts.trainingCell = 2;
    opts.maxFunEvals  = 500;
    opts.verbose      = true;
    opts.P_burst_clip = 25;     % bar - ignore data above this (post-vent)
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'params'), fullfile(root,'scripts'));

    % --- load Phase 2 fitted params + data ----------------------------
    pFile = fullfile(root,'params','etp_params_nmc.mat');
    assert(isfile(pFile), 'Missing %s - run Phase 2 first', pFile);
    p = load(pFile).p_fitted;

    dFile = fullfile(root,'data','gulsoy_parsed.mat');
    assert(isfile(dFile), 'Missing %s - run convert_gulsoy_to_struct first', dFile);
    data = load(dFile).data;
    d = data(opts.trainingCell);

    fprintf('\n=== Phase 3 pressure sub-model fit ===\n');
    fprintf('  Training cell: #%d (%s)\n', opts.trainingCell, d.cell_id);

    % --- enable the pressure model ------------------------------------
    p.pressure.enable = true;

    % --- prepare measured pressure on the SLOW time grid --------------
    % P_internal is on the fast grid (10 kHz); we need it on the slow
    % grid (10 Hz) to match the model output from etp_pure_matlab.
    P_meas_slow = interp1(d.t_fast, d.P_internal, d.t, 'linear', NaN);

    % Identify the pre-vent window: P < P_burst_clip bar AND before
    % the first time P crosses that threshold.
    % Vent = first threshold crossing, OR peak pressure if threshold never reached
    iVent = find(P_meas_slow > opts.P_burst_clip, 1, 'first');
    if isempty(iVent)
        [~, iVent] = max(P_meas_slow);  % fallback: peak pressure = vent onset
    end
    pv_mask = (1:numel(P_meas_slow))' <= iVent;
    pv_mask = pv_mask & ~isnan(P_meas_slow);
    fprintf('  Pre-vent window: %d of %d samples (P < %.0f bar)\n', ...
            sum(pv_mask), numel(pv_mask), opts.P_burst_clip);

    % --- show the measured P range ------------------------------------
    P_pv = P_meas_slow(pv_mask);
    fprintf('  P_meas range in window: [%.2f, %.2f] bar\n', min(P_pv), max(P_pv));

    % --- set up fmincon -----------------------------------------------
    %   x = [nu_vap, nu_e, nu_O2, V_int]
    x0 = [0.0;  2.0;  1.0;  2.5e-6];
    lb = [0.0;  0.01;   0.01;   1.0e-6];  % Phase 3c lower bounds
    ub = [0.0;  10.0;  5.0;   5.0e-6];   % Phase 3c upper bounds

    fprintf('\n  Initial guesses:\n');
    labels = {'nu_vap','nu_e','nu_O2','V_int [mL]'};
    scales = [1, 1, 1, 1e6];
    for k = 1:4
        fprintf('    %-12s : [%.2f, %.2f],  init %.2f\n', ...
            labels{k}, lb(k)*scales(k), ub(k)*scales(k), x0(k)*scales(k));
    end

    objFn = @(x) objective_P(x, p, d, P_meas_slow, pv_mask);

    optsFmin = optimoptions('fmincon','Display','iter', ...
        'MaxFunctionEvaluations', opts.maxFunEvals, ...
        'OptimalityTolerance', 1e-6, ...
        'FiniteDifferenceStepSize', [0.05; 0.05; 0.05; 1e-7]);

    fprintf('\n  Running fmincon...\n\n');
    [x_fit, fval, exitflag] = fmincon(objFn, x0, [], [], [], [], ...
                                       lb, ub, [], optsFmin);
    rmse_P = sqrt(fval);

    % --- update params ------------------------------------------------
    p.pressure.nu_vap = x_fit(1);
    p.pressure.nu_e   = x_fit(2);
    p.pressure.nu_O2  = x_fit(3);
    p.pressure.V_int  = x_fit(4);

    % --- summary ------------------------------------------------------
    fprintf('\n--- Phase 3 results ---\n');
    fprintf('  nu_vap  : %.2f -> %.2f\n', 1.0, x_fit(1));
    fprintf('  nu_e    : %.2f -> %.2f\n', 1.5, x_fit(2));
    fprintf('  nu_O2   : %.2f -> %.2f\n', 0.5, x_fit(3));
    fprintf('  V_int   : %.2f -> %.2f mL\n', 2.0, x_fit(4)*1e6);
    fprintf('  RMSE(P_internal, pre-vent) on training cell: %.3f bar\n', rmse_P);

    if rmse_P > 2
        warning('Pre-vent pressure RMSE = %.2f bar > 2 bar target (plan section 4.1)', rmse_P);
    else
        fprintf('  -> PASS: under 2 bar target\n');
    end

    % --- run one final sim for plotting --------------------------------
    out = simulate_with_pressure(p, d);
    P_pred = out.P / 1e5;   % Pa -> bar
    P_pred_slow = interp1(out.t, P_pred, d.t, 'linear', NaN);

    % --- plot ----------------------------------------------------------
    fig = figure('Color','w','Position',[100 100 1000 600]);
    tlo = tiledlayout(fig, 2, 1, 'TileSpacing','compact','Padding','compact');
    title(tlo, sprintf('Phase 3 fit: cell %s  (P RMSE=%.3f bar)', ...
        d.cell_id, rmse_P), 'FontWeight','bold','Interpreter','none');

    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
    plot(ax1, d.t/60, P_meas_slow, 'k-', 'LineWidth',1.5, 'DisplayName','measured');
    plot(ax1, d.t/60, P_pred_slow, 'r--', 'LineWidth',1.5, 'DisplayName','model');
    if iVent <= numel(d.t)
        xline(ax1, d.t(iVent)/60, 'b:', 'pre-vent end');
    end
    ylabel(ax1, 'P_{internal} [bar]'); legend(ax1,'Location','northwest');
    title(ax1, 'Internal gas pressure');

    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
    resid = P_pred_slow - P_meas_slow;
    plot(ax2, d.t/60, resid, 'b-', 'LineWidth',1);
    yline(ax2, 0, 'k-'); yline(ax2, [-2 2], 'b:', {'-2 bar','+2 bar'});
    ylabel(ax2, 'P_{pred} - P_{meas} [bar]'); xlabel(ax2, 't [min]');
    title(ax2, 'Pressure residual');

    linkaxes([ax1 ax2], 'x');

    outDir = fullfile(root,'results');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    savePath = fullfile(outDir, 'phase3_pressure_fit.png');
    exportgraphics(fig, savePath, 'Resolution', 130);
    fprintf('  Saved plot: %s\n', savePath);

    % --- save ----------------------------------------------------------
    outPath = fullfile(root,'params','etp_params_nmc_p3.mat');
    p_fitted = p;
    save(outPath, 'p_fitted');
    fprintf('  Saved fitted params: %s\n', outPath);

    results.x_fit    = x_fit;
    results.rmse_P   = rmse_P;
    results.p_fitted = p;
    results.exitflag = exitflag;
end

% =====================================================================
function L = objective_P(x, p, d, P_meas_slow, pv_mask)
    p.pressure.nu_vap = x(1);
    p.pressure.nu_e   = x(2);
    p.pressure.nu_O2  = x(3);
    p.pressure.V_int  = x(4);

    try
        out = simulate_with_pressure(p, d);
        P_pred_bar = out.P / 1e5;   % Pa -> bar
        P_pred_slow = interp1(out.t, P_pred_bar, d.t, 'linear', NaN);
        err = P_pred_slow(pv_mask) - P_meas_slow(pv_mask);
        mask2 = ~isnan(err);
        if sum(mask2) < 0.3 * sum(pv_mask)
            L = 1e9;
        else
            L = mean(err(mask2).^2);
        end
    catch
        L = 1e9;
    end
end

function out = simulate_with_pressure(p, d)
    p.pressure.enable = true;
    p.cell.T0 = d.T_internal(1) + 273.15;
    if ~isempty(d.heater_power)
        p.Qext_profile.time  = d.t;
        p.Qext_profile.value = d.heater_power;
    end
    p.solver.tFinal = d.t(end);
    out = etp_pure_matlab(p, 'verbose', false);
end
