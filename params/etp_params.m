function p = etp_params()
% ETP_PARAMS  Returns a struct of all parameters for the coupled
% Electro-Thermal-Pressure (ETP) model of 21700 NMC thermal runaway.
%
% These are the Phase-1 / Phase-2 starting values per the development plan.
% They mix:
%   - Hatchard-Kim reaction kinetics (Kim et al. 2007, Table 2) — LCO-fitted,
%     re-fit to NMC in Phase 2.2
%   - 21700 NMC physical constants (Gulsoy et al. 2024 cell, refined in Phase 2)
%   - Pressure sub-model coefficients (Phase 3 — stubbed with chemistry-based
%     initial guesses; fit in Phase 3.3)
%
% Usage:
%   p = etp_params();           % returns the param struct
%   save('params/etp_params.mat','p');   % cache for fast load
%
% Author: <your name>, 2026.  License: MIT.

% ------------------------------------------------------------------------
% Universal constants
% ------------------------------------------------------------------------
p.Ru     = 8.314;        % J/(mol*K)  Universal gas constant
p.sigma  = 5.670374e-8;  % W/(m^2*K^4) Stefan-Boltzmann
p.T_amb  = 298.15;       % K (= 25 degC) — ambient / heater off-state

% ------------------------------------------------------------------------
% Cell geometry & thermophysical (21700 NMC, Gulsoy et al. 2024 cell)
% ------------------------------------------------------------------------
p.cell.R       = 0.0105;  % m   — 21700 radius
p.cell.L       = 0.0700;  % m   — 21700 length
p.cell.m       = 0.069;   % kg  — typical 21700 NMC mass (refine from datasheet)
p.cell.cp      = 1100;    % J/(kg*K) — effective heat capacity
p.cell.As      = 2*pi*p.cell.R*p.cell.L + 2*pi*p.cell.R^2;  % m^2 — cylinder + ends
p.cell.V       = pi*p.cell.R^2*p.cell.L;                    % m^3 — cell volume
p.cell.rho     = p.cell.m / p.cell.V;                       % kg/m^3 — derived
p.cell.eps     = 0.80;    % surface emissivity (will fit in Phase 2.3)
p.cell.h_conv  = 15;      % W/(m^2*K) natural convection (will fit in Phase 2.3)

% Initial cell temperature for Phase-1 sanity check (overridden by data load)
p.cell.T0      = 298.15;  % K (25 degC)

% ------------------------------------------------------------------------
% Reaction kinetics (Hatchard 2001 / Kim 2007 — LCO baseline)
% Order: 1=SEI, 2=AnE (anode-electrolyte), 3=Cat (cathode), 4=Elec (electrolyte)
% ------------------------------------------------------------------------
% Pre-exponential frequency factors A_i [1/s]
p.rxn.A     = [1.667e15;  2.5e13;   6.667e13; 5.14e25];
% Activation energies Ea_i [J/mol]
p.rxn.Ea    = [1.3508e5;  1.3508e5; 1.396e5;  2.74e5];
% Specific enthalpies of reaction Δh_i [J/kg]
p.rxn.dH    = [2.57e5;    1.714e6;  3.14e5;   1.55e5];
% Specific mass of reactants w_i [kg/m^3]
p.rxn.w     = [610.4;     610.4;    1221.0;   406.9];
% Reaction-order exponents n1_i, n2_i (see Kim eq. 4)
p.rxn.n1    = [1; 1; 1; 1];
p.rxn.n2    = [0; 0; 1; 0];
% Sign on dx/dt: -1 (consumed) for SEI/AnE/Elec; +1 (progress) for cathode
p.rxn.c     = [-1; -1; +1; -1];
% Initial reaction-extent values x_i(0) per Kim Table 2
p.rxn.x0    = [0.15;      0.75;     0.04;     1.0];
% SEI thickness scale for AnE modifier g_AnE(t) = exp(-t_SEI/t_SEI,0)
p.rxn.tSEI0 = 0.033;
p.rxn.tSEI_init = 0.033;  % initial SEI thickness

% Convenience labels
p.rxn.labels = {'SEI','AnE','Cat','Elec'};

% ------------------------------------------------------------------------
% External heat input (Phase-1 sanity-check default = 30 W step)
% Override by setting p.Qext_profile.time / .value before running.
% ------------------------------------------------------------------------
p.Qext_profile.time  = [0; 1e9];    % s — step that lasts the whole sim
p.Qext_profile.value = [30; 30];    % W

% ------------------------------------------------------------------------
% Solver settings (per development plan §1.3)
% ------------------------------------------------------------------------
p.solver.name    = 'ode15s';
p.solver.relTol  = 1e-6;
p.solver.absTol  = 1e-9;
p.solver.maxStep = 0.5;     % s
p.solver.tFinal  = 1800;    % s — 30 minutes per §1.4 spec

% ------------------------------------------------------------------------
% Pressure sub-model (Phase 3 — stubbed; off until you set p.pressure.enable=true)
% ------------------------------------------------------------------------
p.pressure.enable = false;        % Phase 1 = OFF.  Set true in Phase 3.
p.pressure.V_int  = 2.0e-6;       % m^3 (2.0 mL) — fit in Phase 3.3
p.pressure.P0     = 1.013e5;      % Pa (1 atm) — initial cell internal pressure
p.pressure.P_burst= 30e5;         % Pa (30 bar) — burst-disk threshold per plan §3.4

% Stoichiometric coefficients (Phase 3 initial guesses per plan §3.3)
p.pressure.nu_vap = 1.0;   % mol gas per mol electrolyte vapourised
p.pressure.nu_e   = 1.5;   % mol gas per mol electrolyte decomposed
p.pressure.nu_O2  = 0.5;   % mol O2 per mol cathode active material delithiated

% Initial moles available for gas-producing reactions
% (derived from cell mass and chemistry; refine in Phase 3 with cell teardown data)
% Crude estimates — electrolyte mass ~15% of cell, MW ~100 g/mol
p.pressure.m_e0   = 0.15 * p.cell.m;   % kg liquid electrolyte initially
p.pressure.MW_e   = 0.100;             % kg/mol average electrolyte MW
p.pressure.n_e0   = p.pressure.m_e0 / p.pressure.MW_e;   % mol electrolyte
% Cathode active material ~30% of cell mass, MW ~96 g/mol for NMC
p.pressure.m_c0   = 0.30 * p.cell.m;
p.pressure.MW_c   = 0.096;
p.pressure.n_c0   = p.pressure.m_c0 / p.pressure.MW_c;

% Electrolyte vapourisation rate constant (simplified)
% k_vap(T) = exp((T-T_vap_on)/T_vap_scale) for T < T_vap_off, else 0
p.pressure.T_vap_on    = 353.15;  % K (80 degC) — vapourisation turn-on
p.pressure.T_vap_off   = 473.15;  % K (200 degC) — liquid exhausted by here
p.pressure.T_vap_scale = 20;      % K — scale of turn-on
p.pressure.k_vap_pre   = 0;      % disabled: see Phase 3 notes    % 1/s — vapourisation rate pre-factor

end
