function results = run_phase2_param_id(varargin)
% RUN_PHASE2_PARAM_ID  Phase 2 of the development plan:
%   2.2  Re-fit reaction kinetics + enthalpies to NMC (training cell)
%   2.3  Identify thermal loss coefficient h and emissivity epsilon
%
% Usage:
%   results = run_phase2_param_id()
%   results = run_phase2_param_id('trainingCell', 2)
%   results = run_phase2_param_id('boundsMode', 'wide')      % wider dH_cat
%   results = run_phase2_param_id('boundsMode', 'wide', 'fitAnE', true)
%
% Options:
%   trainingCell    1, 2, or 3                       (default 1)
%   boundsMode      'tight' | 'wide'                 (default 'tight')
%                   tight: per dev plan (Ea ±15%, dH ×0.7..×1.4)
%                   wide:  NMC-adjusted (Ea ±20%, dH_cat ×0.7..×2.0)
%   fitAnE          include AnE enthalpy in fit     (default false)
%                   AnE bounds are always wide (×0.7..×1.5)
%                   because AnE is the dominant heat producer
%   skipKineticsFit only do thermal-loss ID (2.3)    (default false)
%   maxFunEvals     fmincon evaluation cap           (default 500)
%
% Author: <your name>, 2026. License: MIT.

    opts.skipKineticsFit = false;
    opts.trainingCell    = 1;
    opts.maxFunEvals     = 500;
    opts.boundsMode      = 'tight';
    opts.fitAnE          = false;
    opts.verbose         = true;
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'params'), fullfile(root,'scripts'));

    % --- load data -----------------------------------------------------
    parsedPath = fullfile(root,'data','gulsoy_parsed.mat');
    if ~isfile(parsedPath)
        rawPath = fullfile(root,'data','TR_dataTable.mat');
        if isfile(rawPath)
            fprintf('Phase 2: parsed dataset not found, converting from raw...\n');
            convert_gulsoy_to_struct(rawPath, parsedPath);
        else
            error(['Phase 2: dataset not found.\n', ...
                   '  Place TR_dataTable.mat under data/ and re-run.\n', ...
                   '  Or run convert_gulsoy_to_struct() with explicit paths.']);
        end
    end
    S = load(parsedPath);
    data = S.data;
    n_cells = numel(data);
    assert(n_cells >= 1, 'No test records found in %s', parsedPath);
    assert(opts.trainingCell <= n_cells, ...
        'Training cell index %d > available cells %d', opts.trainingCell, n_cells);

    fprintf('\n=== Phase 2 parameter identification ===\n');
    fprintf('  Dataset: %d test(s), training on cell #%d\n\n', ...
            n_cells, opts.trainingCell);

    p = etp_params();
    d_train = data(opts.trainingCell);

    % ---- 2.3: identify h and epsilon (pre-onset) ---------------------
    [h_fit, eps_fit, info_thermal] = identify_thermal_losses(d_train, p, opts);
    p.cell.h_conv = h_fit;
    p.cell.eps    = eps_fit;

    % ---- 2.2: re-fit cathode kinetics to NMC -------------------------
    if opts.skipKineticsFit
        fprintf('\n  Skipping kinetics fit (skipKineticsFit=true).\n');
        params_fit = [p.rxn.A(3); p.rxn.Ea(3); p.rxn.dH(3)];
        if opts.fitAnE, params_fit(4) = p.rxn.dH(2); end
        rmse_train = compute_rmse(p, d_train);
    else
        [params_fit, rmse_train, info_kin] = fit_cathode_kinetics(p, d_train, opts);
        p.rxn.A(3)  = params_fit(1);
        p.rxn.Ea(3) = params_fit(2);
        p.rxn.dH(3) = params_fit(3);
        if opts.fitAnE
            p.rxn.dH(2) = params_fit(4);
        end
    end

    % ---- summary ------------------------------------------------------
    fprintf('\n--- Phase 2 results ---\n');
    fprintf('  bounds  : %s\n', opts.boundsMode);
    fprintf('  fit AnE : %s\n', tf(opts.fitAnE));
    fprintf('  h_conv  : %.3f -> %.3f  W/(m^2*K)\n', 15.0, h_fit);
    fprintf('  eps     : %.3f -> %.3f\n', 0.80, eps_fit);
    fprintf('  A_cat   : %.3e -> %.3e  1/s\n', 6.667e13, params_fit(1));
    fprintf('  Ea_cat  : %.4e -> %.4e J/mol (Δ = %+.1f%%)\n', ...
            1.396e5, params_fit(2), 100*(params_fit(2)/1.396e5 - 1));
    fprintf('  dH_cat  : %.3e -> %.3e J/kg  (×%.2f literature)\n', ...
            3.14e5, params_fit(3), params_fit(3)/3.14e5);
    if opts.fitAnE
        fprintf('  dH_AnE  : %.3e -> %.3e J/kg  (×%.2f literature)\n', ...
                1.714e6, params_fit(4), params_fit(4)/1.714e6);
    end
    fprintf('  RMSE(T_internal) on training cell: %.2f K\n', rmse_train);

    if rmse_train > 30
        warning(['Training-cell RMSE = %.1f K > 30 K target.  ', ...
                 'Check parameter bounds and/or column mapping in load_gulsoy_data.'], ...
                rmse_train);
    end

    % ---- save -------------------------------------------------------
    outPath = fullfile(root,'params','etp_params_nmc.mat');
    p_fitted = p;
    save(outPath, 'p_fitted');
    fprintf('\nSaved fitted params: %s\n', outPath);

    results.params_fit = params_fit;
    results.h_fit      = h_fit;
    results.eps_fit    = eps_fit;
    results.rmse_train = rmse_train;
    results.p_fitted   = p;
    results.opts       = opts;
end

function s = tf(b), if b, s='true'; else, s='false'; end, end

% =====================================================================
function [h_fit, eps_fit, info] = identify_thermal_losses(d, p, opts)
% Use the pre-onset portion (T < 100 degC) where chemistry is negligible.
% Fit (h, eps) by minimising T error on the slow-rise segment.

    % Find pre-onset window: from t=0 until T_internal crosses 100 degC
    T_C = d.T_internal;
    t   = d.t;
    iPre = find(T_C < 100, 1, 'last');
    if isempty(iPre) || iPre < 50
        warning('Pre-onset window too short for thermal-loss ID; using defaults');
        h_fit = p.cell.h_conv; eps_fit = p.cell.eps;
        info = struct('window_size', 0);
        return;
    end
    t_pre = t(1:iPre);
    T_pre_C = T_C(1:iPre);
    fprintf('  [2.3] Thermal-loss ID: pre-onset window 1..%d (%.0f s)\n', ...
            iPre, t_pre(end));

    if isempty(d.heater_power)
        warning('No heater_power in dataset; using p.Qext_profile default');
        Qext_fn = @(tq) interp1(p.Qext_profile.time, p.Qext_profile.value, ...
                                tq, 'previous', p.Qext_profile.value(end));
    else
        Qext_fn = @(tq) interp1(d.t, d.heater_power, tq, 'linear', d.heater_power(end));
    end

    % Solve simplified pre-onset balance: m*cp*dT/dt = Qext - Qloss(T; h, eps)
    obj = @(x) loss_thermal(x(1), x(2), t_pre, T_pre_C+273.15, p, Qext_fn);

    x0 = [15; 0.8];
    lb = [5;  0.5];
    ub = [40; 0.95];
    optsFmin = optimoptions('fmincon','Display','none', ...
                            'MaxFunctionEvaluations', 200);
    try
        x_fit = fmincon(obj, x0, [], [], [], [], lb, ub, [], optsFmin);
    catch ME
        warning('fmincon for thermal losses failed: %s\n  Using defaults.', ME.message);
        x_fit = x0;
    end
    h_fit = x_fit(1);
    eps_fit = x_fit(2);
    info.window_size = iPre;
    info.rmse = sqrt(obj(x_fit));
end

function L = loss_thermal(h, eps, t_grid, T_meas_K, p, Qext_fn)
    p.cell.h_conv = h;
    p.cell.eps    = eps;
    % integrate pre-onset segment with reaction kinetics frozen at IC
    odeopts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1);
    y0 = [p.rxn.x0(1); p.rxn.x0(2); p.rxn.x0(3); p.rxn.x0(4); ...
          T_meas_K(1); p.rxn.tSEI_init];
    rhsfn = @(tt,yy) thermal_only_rhs(tt,yy,p,Qext_fn);
    try
        [~, ysol] = ode15s(rhsfn, [t_grid(1) t_grid(end)], y0, odeopts);
        T_pred = ysol(:,5);
        % interpolate predicted onto t_grid
        sol = ode15s(rhsfn, [t_grid(1) t_grid(end)], y0, odeopts);
        T_pred_grid = deval(sol, t_grid, 5);
        L = mean((T_pred_grid(:) - T_meas_K(:)).^2);
    catch
        L = 1e9;
    end
end

function dy = thermal_only_rhs(t, y, p, Qext_fn)
% Same RHS as etp_pure_matlab but with reaction rates near-zero
% (we're in the pre-onset window so they don't matter).
    T = y(5);
    Qext = Qext_fn(t);
    Qconv = p.cell.h_conv * p.cell.As * (T - p.T_amb);
    Qrad  = p.cell.eps * p.sigma * p.cell.As * (T^4 - p.T_amb^4);
    dy = zeros(6,1);
    dy(5) = (Qext - Qconv - Qrad) / (p.cell.m * p.cell.cp);
end

% =====================================================================
function [params_fit, rmse_train, info] = fit_cathode_kinetics(p, d, opts)
% Fit cathode A, Ea, dH (and optionally AnE dH) with mode-dependent bounds.
    fprintf('\n  [2.2] Reaction kinetics fit (cell #%d, bounds=%s, fitAnE=%s)\n', ...
            opts.trainingCell, opts.boundsMode, mat2str(opts.fitAnE));

    x0 = [p.rxn.A(3); p.rxn.Ea(3); p.rxn.dH(3)];

    switch opts.boundsMode
        case 'tight'    % per development plan
            lb_mul = [0.30; 0.85; 0.70];
            ub_mul = [3.00; 1.15; 1.40];
        case 'wide'     % NMC-adjusted (dH_cat up to 2x literature)
            lb_mul = [0.30; 0.80; 0.70];
            ub_mul = [3.00; 1.20; 2.00];
        otherwise
            error('Unknown boundsMode "%s"', opts.boundsMode);
    end

    if opts.fitAnE
        x0(4)     = p.rxn.dH(2);
        lb_mul(4) = 0.70;
        ub_mul(4) = 1.50;          % AnE always wide; it's the dominant heat source
    end

    lb = x0 .* lb_mul;
    ub = x0 .* ub_mul;

    fprintf('     A_cat:  [%.2e, %.2e],  init %.2e\n', lb(1), ub(1), x0(1));
    fprintf('     Ea_cat: [%.3e, %.3e],  init %.3e\n', lb(2), ub(2), x0(2));
    fprintf('     dH_cat: [%.2e, %.2e],  init %.2e\n', lb(3), ub(3), x0(3));
    if opts.fitAnE
        fprintf('     dH_AnE: [%.2e, %.2e],  init %.2e\n', lb(4), ub(4), x0(4));
    end

    objFn = @(x) objective_T_internal(x, p, d, opts.fitAnE);
    optsFmin = optimoptions('fmincon','Display','iter', ...
        'MaxFunctionEvaluations', opts.maxFunEvals, ...
        'OptimalityTolerance', 1e-4);

    [params_fit, fval] = fmincon(objFn, x0, [], [], [], [], lb, ub, [], optsFmin);
    rmse_train = sqrt(fval);
    info.fval = fval;
end

function L = objective_T_internal(x, p, d, fitAnE)
% L = mean-square error between predicted and measured T_internal
    p_try = p;
    p_try.rxn.A(3)  = x(1);
    p_try.rxn.Ea(3) = x(2);
    p_try.rxn.dH(3) = x(3);
    if fitAnE
        p_try.rxn.dH(2) = x(4);
    end
    p_try.cell.T0   = d.T_internal(1) + 273.15;

    if ~isempty(d.heater_power)
        p_try.Qext_profile.time  = d.t;
        p_try.Qext_profile.value = d.heater_power;
    end
    p_try.solver.tFinal = d.t(end);

    try
        out = etp_pure_matlab(p_try, 'verbose', false);
        T_pred_C = interp1(out.t, out.T, d.t) - 273.15;
        err = T_pred_C(:) - d.T_internal(:);
        mask = ~isnan(err);
        if sum(mask) < 0.5 * numel(err)
            L = 1e9;
        else
            L = mean(err(mask).^2);
        end
    catch
        L = 1e9;
    end
end

function rmse = compute_rmse(p, d)
    p_try = p;
    p_try.cell.T0   = d.T_internal(1) + 273.15;
    if ~isempty(d.heater_power)
        p_try.Qext_profile.time  = d.t;
        p_try.Qext_profile.value = d.heater_power;
    end
    p_try.solver.tFinal = d.t(end);
    out = etp_pure_matlab(p_try, 'verbose', false);
    T_pred_C = interp1(out.t, out.T, d.t) - 273.15;
    err = T_pred_C(:) - d.T_internal(:);
    rmse = sqrt(mean(err(~isnan(err)).^2));
end
