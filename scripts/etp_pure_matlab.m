function out = etp_pure_matlab(p, varargin)
% ETP_PURE_MATLAB  Pure-MATLAB ode15s integration of the coupled
% Electro-Thermal-Pressure model.  Functionally identical to the Simulink
% build but ~20x faster to iterate on and easier to debug.
%
% Usage:
%   p   = etp_params();
%   out = etp_pure_matlab(p);                  % default 30 W, 1800 s
%   out = etp_pure_matlab(p, 'tFinal', 3600);  % override sim length
%   plot(out.t/60, out.T - 273.15); xlabel('t [min]'); ylabel('T [degC]');
%
% State vector y(t) (six variables):
%   y(1) = x_sei     SEI fraction remaining             [-]
%   y(2) = c_ne      anode reactive content fraction    [-]
%   y(3) = alpha     cathode conversion progress        [-]
%   y(4) = c_e       electrolyte fraction remaining     [-]
%   y(5) = T         cell temperature                   [K]
%   y(6) = tSEI      SEI thickness (drives g_AnE)       [-]
%
% Optional Phase-3 pressure state appended when p.pressure.enable == true:
%   y(7) = n_g       moles of gas in V_int              [mol]
%
% Returns:
%   out.t            (Nx1) time vector [s]
%   out.y            (NxNs) full state
%   out.T            (Nx1) cell T [K]
%   out.x            (Nx4) reaction extents [SEI, AnE, Cat, Elec]
%   out.Q_rxn        (Nx1) total reaction heat [W]
%   out.Q_ext        (Nx1) external heat input [W]
%   out.Q_loss       (Nx1) heat loss to ambient [W]
%   out.P            (Nx1) internal pressure [Pa]  (if enabled)
%   out.events       struct of detected events (onset, vent, peak)
%
% Author: <your name>, 2026. License: MIT.

    % ---- parse options -------------------------------------------------
    opts.tFinal     = p.solver.tFinal;
    opts.verbose    = true;
    opts            = parse_opts(opts, varargin);

    % ---- initial state -------------------------------------------------
    y0 = [p.rxn.x0(1);     % x_sei
          p.rxn.x0(2);     % c_ne
          p.rxn.x0(3);     % alpha
          p.rxn.x0(4);     % c_e
          p.cell.T0;       % T  [K]
          p.rxn.tSEI_init];% tSEI

    if p.pressure.enable
        % derive initial n_g from P0, V_int, T0 (subtract so P(0) = P0)
        n_g0 = p.pressure.P0 * p.pressure.V_int / (p.Ru * p.cell.T0);
        y0(7) = n_g0;
        p.pressure.n_g0 = n_g0;  % stash for diagnostic
    end

    % ---- solver --------------------------------------------------------
    odeopts = odeset('RelTol', p.solver.relTol, ...
                     'AbsTol', p.solver.absTol, ...
                     'MaxStep', p.solver.maxStep, ...
                     'NonNegative', 1:4);    % reaction extents >= 0

    if opts.verbose
        fprintf('etp_pure_matlab: integrating with ode15s, tFinal=%.0f s\n', opts.tFinal);
    end
    tic;
    [t, y] = ode15s(@(t,y) rhs(t,y,p), [0 opts.tFinal], y0, odeopts);
    wallTime = toc;
    if opts.verbose
        fprintf('  done. %d steps, wall=%.2f s, RTF=%.4f\n', ...
                numel(t), wallTime, wallTime/opts.tFinal);
    end

    % ---- post-process diagnostics --------------------------------------
    N = numel(t);
    Qrxn  = zeros(N,1);
    Qext  = zeros(N,1);
    Qloss = zeros(N,1);
    for k = 1:N
        [Qrxn(k), Qext(k), Qloss(k)] = diag_q(t(k), y(k,:).', p);
    end

    out.t        = t;
    out.y        = y;
    out.T        = y(:,5);
    out.tSEI     = y(:,6);
    out.x        = y(:,1:4);
    out.Q_rxn    = Qrxn;
    out.Q_ext    = Qext;
    out.Q_loss   = Qloss;

    if p.pressure.enable
        out.n_g  = y(:,7);
        out.P    = y(:,7) .* p.Ru .* y(:,5) ./ p.pressure.V_int;   % Pa
    end

    % ---- event detection ----------------------------------------------
    out.events = detect_events(out, p);
end

% =====================================================================
%  RHS — coupled state derivatives
% =====================================================================
function dy = rhs(t, y, p)
    x_sei = y(1);
    c_ne  = y(2);
    alpha = y(3);
    c_e   = y(4);
    T     = y(5);
    tSEI  = y(6);

    % reaction rates per Kim eq (4):  R_i = A * exp(-Ea/RT) * f(x) * g(t)
    R = zeros(4,1);
    R(1) = p.rxn.A(1) * exp(-p.rxn.Ea(1)/(p.Ru*T)) * max(x_sei,0);
    g_AnE = exp(-tSEI / p.rxn.tSEI0);
    R(2) = p.rxn.A(2) * exp(-p.rxn.Ea(2)/(p.Ru*T)) * max(c_ne,0) * g_AnE;
    R(3) = p.rxn.A(3) * exp(-p.rxn.Ea(3)/(p.Ru*T)) * max(alpha,0) * max(1-alpha,0);
    R(4) = p.rxn.A(4) * exp(-p.rxn.Ea(4)/(p.Ru*T)) * max(c_e,0);

    % reaction-extent ODEs:  dx/dt = c_i * R_i
    dx_sei = p.rxn.c(1) * R(1);
    dc_ne  = p.rxn.c(2) * R(2);
    dalpha = p.rxn.c(3) * R(3);
    dc_e   = p.rxn.c(4) * R(4);
    dtSEI  = R(2);     % SEI grows at rate of AnE reaction (Kim eq 5')

    % volumetric heat generation: Q_i [W] = w_i * V_cell * dH_i * R_i
    Qrxn = sum(p.rxn.w .* p.cell.V .* p.rxn.dH .* R);

    % external heat input (piecewise-linear interp on profile)
    Qext = interp1(p.Qext_profile.time, p.Qext_profile.value, t, ...
                   'previous', p.Qext_profile.value(end));

    % heat loss: convection + radiation to ambient
    Qconv = p.cell.h_conv * p.cell.As * (T - p.T_amb);
    Qrad  = p.cell.eps * p.sigma * p.cell.As * (T^4 - p.T_amb^4);
    Qloss = Qconv + Qrad;

    % lumped energy balance:  m*cp*dT/dt = Qrxn + Qext - Qloss
    dT = (Qrxn + Qext - Qloss) / (p.cell.m * p.cell.cp);

    dy = [dx_sei; dc_ne; dalpha; dc_e; dT; dtSEI];

    % ----- optional pressure state ---------------------------------
    if p.pressure.enable
        % electrolyte vapourisation
        if T < p.pressure.T_vap_off
            k_vap = p.pressure.k_vap_pre * ...
                    exp((T - p.pressure.T_vap_on)/p.pressure.T_vap_scale);
        else
            k_vap = 0;
        end
        % remaining liquid electrolyte (rough: fraction = c_e)
        m_e = p.pressure.m_e0 * max(c_e,0);

        dn_g = p.pressure.nu_vap * k_vap * (m_e / p.pressure.MW_e) ...
             + p.pressure.nu_e   * abs(dc_e)   * p.pressure.n_e0 ...
             + p.pressure.nu_O2  * max(dalpha,0) * p.pressure.n_c0;
        dy(7,1) = dn_g;
    end
end

% =====================================================================
%  Diagnostic Q breakdown (computed post-hoc for plotting)
% =====================================================================
function [Qrxn, Qext, Qloss] = diag_q(t, y, p)
    x_sei = y(1); c_ne = y(2); alpha = y(3); c_e = y(4); T = y(5); tSEI = y(6);
    R(1) = p.rxn.A(1)*exp(-p.rxn.Ea(1)/(p.Ru*T))*max(x_sei,0);
    R(2) = p.rxn.A(2)*exp(-p.rxn.Ea(2)/(p.Ru*T))*max(c_ne,0)*exp(-tSEI/p.rxn.tSEI0);
    R(3) = p.rxn.A(3)*exp(-p.rxn.Ea(3)/(p.Ru*T))*max(alpha,0)*max(1-alpha,0);
    R(4) = p.rxn.A(4)*exp(-p.rxn.Ea(4)/(p.Ru*T))*max(c_e,0);
    Qrxn = sum(p.rxn.w .* p.cell.V .* p.rxn.dH .* R(:));
    Qext = interp1(p.Qext_profile.time, p.Qext_profile.value, t, ...
                   'previous', p.Qext_profile.value(end));
    Qloss = p.cell.h_conv*p.cell.As*(T-p.T_amb) ...
          + p.cell.eps*p.sigma*p.cell.As*(T^4-p.T_amb^4);
end

% =====================================================================
%  Event detection — onset, peak, vent
% =====================================================================
function ev = detect_events(out, p)
    % onset: dT/dt first crosses 1 K/s
    dT = gradient(out.T, out.t);
    iOnset = find(dT > 1, 1, 'first');
    if ~isempty(iOnset)
        ev.t_onset = out.t(iOnset);
        ev.T_onset = out.T(iOnset);
    else
        ev.t_onset = NaN;  ev.T_onset = NaN;
    end
    [ev.T_peak, iPk] = max(out.T);
    ev.t_peak = out.t(iPk);

    if p.pressure.enable
        iVent = find(out.P > p.pressure.P_burst, 1, 'first');
        if ~isempty(iVent)
            ev.t_vent = out.t(iVent);
            ev.P_vent = out.P(iVent);
        else
            ev.t_vent = NaN;  ev.P_vent = NaN;
        end
    end
end

% =====================================================================
function opts = parse_opts(opts, args)
    for k = 1:2:numel(args)
        opts.(args{k}) = args{k+1};
    end
end
