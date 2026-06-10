function results = run_phase4_sensitivity(varargin)
% RUN_PHASE4_SENSITIVITY  Phase 4 of the development plan:
%   One-at-a-time sensitivity sweeps on three parameters:
%     1. State of Charge (SoC): 50%, 60%, 70%, 80%, 90%, 100%
%     2. Heater power: +/-25% about nominal
%     3. Effective heat capacity: +/-15% about nominal
%
%   Output metrics:
%     - Time-to-onset (dT/dt > 1 K/s)
%     - Pre-vent pressure rise rate (dP/dt averaged over last 60s before onset)
%
%   Produces:
%     - Normalised sensitivity indices S_i = (dy/dx_i)*(x_i/y)
%     - Bar charts for paper Figures 4 and 5
%     - Results table for paper Section 5.3
%
% Prerequisites:
%   - params/etp_params_nmc_p3.mat (Phase 3 output with pressure enabled)
%
% Usage:
%   results = run_phase4_sensitivity();
%
% Author: <your name>, 2026. License: MIT.

    opts.paramsFile = fullfile('params','etp_params_nmc_p3.mat');
    opts.savePlots  = true;
    opts.verbose    = true;
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'params'), fullfile(root,'scripts'));

    assert(isfile(opts.paramsFile), 'Missing %s - run Phase 3 first', opts.paramsFile);
    p_base = load(opts.paramsFile).p_fitted;
    p_base.pressure.enable = true;

    fprintf('\n=== Phase 4 sensitivity analysis ===\n\n');

    % ---- nominal run -------------------------------------------------
    fprintf('Running nominal case...\n');
    [t_onset_nom, dPdt_nom] = run_one(p_base, opts);
    fprintf('  Nominal: t_onset = %.1f s, dP/dt = %.4f bar/s\n\n', ...
            t_onset_nom, dPdt_nom);

    % ==================================================================
    %  Sweep 1: State of Charge
    % ==================================================================
    SoC_values = [0.50, 0.60, 0.70, 0.80, 0.90, 1.00];
    n_soc = numel(SoC_values);
    t_onset_soc = zeros(n_soc, 1);
    dPdt_soc    = zeros(n_soc, 1);

    fprintf('--- Sweep 1: State of Charge ---\n');
    for k = 1:n_soc
        p_try = p_base;
        soc = SoC_values(k);

        % SoC affects initial reaction extents:
        %   c_ne (anode lithium): scales linearly with SoC
        %   alpha (cathode delithiation): higher SoC = more delithiated
        %   nu_O2 effect: more oxygen at higher SoC (already in the model
        %     via alpha initial condition)
        % SoC affects:
        %   c_ne: more intercalated Li at higher SoC -> more AnE fuel
        %   alpha: always 0.04 (reaction progress, not delithiation)
        %   nu_O2 effective: more O2 release at higher SoC (scale linearly)
        p_try.rxn.x0(2) = 0.75 * soc;                   % c_ne scales with SoC
        p_try.rxn.x0(3) = 0.04;                          % alpha always starts here
        p_try.pressure.nu_O2 = p_base.pressure.nu_O2 * (0.5 + 0.5*soc); % more O2 at high SoC

        [t_onset_soc(k), dPdt_soc(k)] = run_one(p_try, opts);
        fprintf('  SoC=%3.0f%%: t_onset=%7.1f s, dP/dt=%.4f bar/s\n', ...
                soc*100, t_onset_soc(k), dPdt_soc(k));
    end

    % ==================================================================
    %  Sweep 2: Heater power (+/-25%)
    % ==================================================================
    Q_nom = p_base.Qext_profile.value(1);
    Q_factors = [0.75, 0.85, 0.95, 1.00, 1.05, 1.15, 1.25];
    n_q = numel(Q_factors);
    t_onset_q = zeros(n_q, 1);
    dPdt_q    = zeros(n_q, 1);

    fprintf('\n--- Sweep 2: Heater power (nominal = %.1f W) ---\n', Q_nom);
    for k = 1:n_q
        p_try = p_base;
        p_try.Qext_profile.value = Q_nom * Q_factors(k) * ...
            ones(size(p_try.Qext_profile.value));

        [t_onset_q(k), dPdt_q(k)] = run_one(p_try, opts);
        fprintf('  Q=%5.1f W (%+3.0f%%): t_onset=%7.1f s, dP/dt=%.4f bar/s\n', ...
                Q_nom*Q_factors(k), (Q_factors(k)-1)*100, ...
                t_onset_q(k), dPdt_q(k));
    end

    % ==================================================================
    %  Sweep 3: Heat capacity (+/-15%)
    % ==================================================================
    cp_nom = p_base.cell.cp;
    cp_factors = [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15];
    n_cp = numel(cp_factors);
    t_onset_cp = zeros(n_cp, 1);
    dPdt_cp    = zeros(n_cp, 1);

    fprintf('\n--- Sweep 3: Heat capacity (nominal = %.1f J/(kg K)) ---\n', cp_nom);
    for k = 1:n_cp
        p_try = p_base;
        p_try.cell.cp = cp_nom * cp_factors(k);

        [t_onset_cp(k), dPdt_cp(k)] = run_one(p_try, opts);
        fprintf('  cp=%6.1f (%+3.0f%%): t_onset=%7.1f s, dP/dt=%.4f bar/s\n', ...
                cp_nom*cp_factors(k), (cp_factors(k)-1)*100, ...
                t_onset_cp(k), dPdt_cp(k));
    end

    % ==================================================================
    %  Compute normalised sensitivity indices at the nominal point
    % ==================================================================
    fprintf('\n--- Normalised sensitivity indices ---\n');
    fprintf('  S_i = (dy/dx_i) * (x_i / y)  evaluated at nominal\n\n');

    % For SoC: use the 80% and 100% points to estimate derivative
    i80  = find(SoC_values == 0.80);
    i100 = find(SoC_values == 1.00);
    if ~isempty(i80) && ~isempty(i100)
        dx_soc = SoC_values(i100) - SoC_values(i80);
        S_soc_tonset = ((t_onset_soc(i100) - t_onset_soc(i80)) / dx_soc) * ...
                       (0.90 / t_onset_nom);
        S_soc_dPdt   = ((dPdt_soc(i100) - dPdt_soc(i80)) / dx_soc) * ...
                       (0.90 / dPdt_nom);
    else
        S_soc_tonset = NaN; S_soc_dPdt = NaN;
    end

    % For heater power: use +/-25% points
    i_lo = find(Q_factors == 0.75);
    i_hi = find(Q_factors == 1.25);
    i_nm = find(Q_factors == 1.00);
    if ~isempty(i_lo) && ~isempty(i_hi)
        dQ = Q_nom * (Q_factors(i_hi) - Q_factors(i_lo));
        S_Q_tonset = ((t_onset_q(i_hi) - t_onset_q(i_lo)) / dQ) * ...
                     (Q_nom / t_onset_nom);
        S_Q_dPdt   = ((dPdt_q(i_hi) - dPdt_q(i_lo)) / dQ) * ...
                     (Q_nom / dPdt_nom);
    else
        S_Q_tonset = NaN; S_Q_dPdt = NaN;
    end

    % For heat capacity: use +/-15% points
    i_lo = find(cp_factors == 0.85);
    i_hi = find(cp_factors == 1.15);
    if ~isempty(i_lo) && ~isempty(i_hi)
        dcp = cp_nom * (cp_factors(i_hi) - cp_factors(i_lo));
        S_cp_tonset = ((t_onset_cp(i_hi) - t_onset_cp(i_lo)) / dcp) * ...
                      (cp_nom / t_onset_nom);
        S_cp_dPdt   = ((dPdt_cp(i_hi) - dPdt_cp(i_lo)) / dcp) * ...
                      (cp_nom / dPdt_nom);
    else
        S_cp_tonset = NaN; S_cp_dPdt = NaN;
    end

    fprintf('  %-20s | S(t_onset) | S(dP/dt)\n', 'Parameter');
    fprintf('  %s\n', repmat('-', 1, 55));
    fprintf('  %-20s | %+8.3f   | %+8.3f\n', 'SoC', S_soc_tonset, S_soc_dPdt);
    fprintf('  %-20s | %+8.3f   | %+8.3f\n', 'Heater power', S_Q_tonset, S_Q_dPdt);
    fprintf('  %-20s | %+8.3f   | %+8.3f\n', 'Heat capacity', S_cp_tonset, S_cp_dPdt);

    % ==================================================================
    %  Store results
    % ==================================================================
    results.nominal.t_onset = t_onset_nom;
    results.nominal.dPdt    = dPdt_nom;

    results.soc.values      = SoC_values;
    results.soc.t_onset     = t_onset_soc;
    results.soc.dPdt        = dPdt_soc;

    results.heater.factors  = Q_factors;
    results.heater.Q_nom    = Q_nom;
    results.heater.t_onset  = t_onset_q;
    results.heater.dPdt     = dPdt_q;

    results.cp.factors      = cp_factors;
    results.cp.cp_nom       = cp_nom;
    results.cp.t_onset      = t_onset_cp;
    results.cp.dPdt         = dPdt_cp;

    results.sensitivity.S_soc_tonset = S_soc_tonset;
    results.sensitivity.S_soc_dPdt   = S_soc_dPdt;
    results.sensitivity.S_Q_tonset   = S_Q_tonset;
    results.sensitivity.S_Q_dPdt     = S_Q_dPdt;
    results.sensitivity.S_cp_tonset  = S_cp_tonset;
    results.sensitivity.S_cp_dPdt    = S_cp_dPdt;

    % ==================================================================
    %  Plots
    % ==================================================================
    if opts.savePlots
        outDir = fullfile(root, 'results');
        if ~exist(outDir,'dir'), mkdir(outDir); end
        plot_soc_sweep(results, outDir);
        plot_sensitivity_bars(results, outDir);
    end

    % Save results
    outPath = fullfile(root, 'results', 'phase4_sensitivity.mat');
    save(outPath, 'results');
    fprintf('\nSaved results: %s\n', outPath);
    fprintf('\n=== Phase 4 complete ===\n');
end

% =====================================================================
%  Run a single simulation and extract metrics
% =====================================================================
function [t_onset, dPdt_avg] = run_one(p, opts)
    p.pressure.enable = true;
    try
        out = etp_pure_matlab(p, 'verbose', false);
    catch ME
        if opts.verbose
            fprintf('    SIM FAILED: %s\n', ME.message);
        end
        t_onset = NaN; dPdt_avg = NaN;
        return;
    end

    T_C = out.T - 273.15;
    P_bar = out.P / 1e5;

    % Onset: smoothed dT/dt > 1 K/s, sustained, after t > 30 s
    dT = movmean(gradient(T_C, out.t), 21);
    above = dT > 1;
    if numel(above) >= 5
        sustained = conv(double(above), ones(5,1)/5, 'same') > 0.99;
    else
        sustained = above;
    end
    valid = sustained & (out.t > 30);
    iOn = find(valid, 1, 'first');

    if isempty(iOn)
        t_onset = NaN;
        dPdt_avg = NaN;
        return;
    end
    t_onset = out.t(iOn);

    % Pre-vent pressure rise rate: average dP/dt over the 60 s before onset
    t_window_start = max(0, t_onset - 60);
    mask = out.t >= t_window_start & out.t <= t_onset;
    if sum(mask) > 5
        dP = gradient(P_bar, out.t);
        dPdt_avg = mean(dP(mask));
    else
        dPdt_avg = NaN;
    end
end

% =====================================================================
%  Plot: SoC sweep (Figure 4 in the paper)
% =====================================================================
function plot_soc_sweep(r, outDir)
    fig = figure('Color','w','Position',[100 100 900 500],'Visible','off');
    tlo = tiledlayout(fig, 1, 2, 'TileSpacing','compact','Padding','compact');
    title(tlo, 'Sensitivity to State of Charge', 'FontWeight','bold');

    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
    bar(ax1, r.soc.values*100, r.soc.t_onset, 'FaceColor', [0.2 0.4 0.8]);
    xlabel(ax1, 'SoC [%]'); ylabel(ax1, 'Time to onset [s]');
    title(ax1, 'Time to onset vs SoC');

    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
    bar(ax2, r.soc.values*100, r.soc.dPdt, 'FaceColor', [0.8 0.3 0.2]);
    xlabel(ax2, 'SoC [%]'); ylabel(ax2, 'Pre-vent dP/dt [bar/s]');
    title(ax2, 'Pressure rise rate vs SoC');

    savePath = fullfile(outDir, 'phase4_soc_sweep.png');
    exportgraphics(fig, savePath, 'Resolution', 150);
    close(fig);
    fprintf('  Saved: %s\n', savePath);
end

% =====================================================================
%  Plot: Normalised sensitivity bar chart (Figure 5 in the paper)
% =====================================================================
function plot_sensitivity_bars(r, outDir)
    S = r.sensitivity;
    params = {'SoC', 'Heater power', 'Heat capacity'};
    S_tonset = [S.S_soc_tonset, S.S_Q_tonset, S.S_cp_tonset];
    S_dPdt   = [S.S_soc_dPdt,   S.S_Q_dPdt,   S.S_cp_dPdt];

    fig = figure('Color','w','Position',[100 100 800 450],'Visible','off');
    tlo = tiledlayout(fig, 1, 2, 'TileSpacing','compact','Padding','compact');
    title(tlo, 'Normalised Sensitivity Indices', 'FontWeight','bold');

    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
    b1 = bar(ax1, categorical(params, params), abs(S_tonset), ...
        'FaceColor', [0.2 0.4 0.8]);
    ylabel(ax1, '|S_i| (time to onset)');
    title(ax1, 'S_i for t_{onset}');
    % Add value labels
    for k = 1:3
        text(ax1, k, abs(S_tonset(k))+0.02, sprintf('%+.2f', S_tonset(k)), ...
            'HorizontalAlignment','center','FontSize',10);
    end

    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
    b2 = bar(ax2, categorical(params, params), abs(S_dPdt), ...
        'FaceColor', [0.8 0.3 0.2]);
    ylabel(ax2, '|S_i| (pressure rise rate)');
    title(ax2, 'S_i for dP/dt');
    for k = 1:3
        text(ax2, k, abs(S_dPdt(k))+0.02, sprintf('%+.2f', S_dPdt(k)), ...
            'HorizontalAlignment','center','FontSize',10);
    end

    savePath = fullfile(outDir, 'phase4_sensitivity_bars.png');
    exportgraphics(fig, savePath, 'Resolution', 150);
    close(fig);
    fprintf('  Saved: %s\n', savePath);
end
