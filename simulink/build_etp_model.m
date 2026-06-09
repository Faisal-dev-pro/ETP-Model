function build_etp_model(modelName)
% BUILD_ETP_MODEL  Programmatically construct the coupled ETP Simulink
% model.  Running this script (re)creates simulink/etp_model.slx from
% scratch, so the model is fully version-controlled via this .m file
% rather than as a binary .slx blob.
%
%   build_etp_model()              % default: builds 'etp_model'
%   build_etp_model('my_test')     % custom name
%
% Design choice
% -------------
% Each reaction is implemented as a MATLAB Function block that takes
% (T, x, tSEI) and returns (dx/dt, Q_i).  This avoids the wiring trap
% of trying to express the Arrhenius+f(x) product with primitive blocks
% (which fails for AnE because the Fcn block does not support
% multi-input expressions like u(1)*exp(-u(2)/...) ).
%
% Phase scope:
%   Phase 1  -- IMPLEMENTED.  4 reactions, lumped thermal balance, Q_ext
%               from base workspace, all signals logged.
%   Phase 3  -- Pressure subsystem is a stub.  Replace the Constant block
%               with a MATLAB Function that reads alpha, c_e, T and
%               integrates dn_g/dt.
%
% Inputs from base workspace (auto-populated if absent):
%   p        -- params struct from etp_params()
%   Qext_ts  -- timeseries (auto-built in model InitFcn)
%
% Outputs (logged to base workspace 'simOut.yout' in this order):
%   1: T [K]    2: x_SEI    3: c_ne    4: alpha    5: c_e
%   6: tSEI     7: Q_rxn    8: P [Pa]
%
% Author: <your name>, 2026. License: MIT.

    if nargin < 1, modelName = 'etp_model'; end

    % --- prepare workspace --------------------------------------------
    if bdIsLoaded(modelName), close_system(modelName, 0); end
    new_system(modelName);
    open_system(modelName);

    % --- configure solver per plan section 1.3 ------------------------
    cs = getActiveConfigSet(modelName);
    set_param(cs, 'SolverType',    'Variable-step');
    set_param(cs, 'Solver',        'ode15s');
    set_param(cs, 'RelTol',        '1e-6');
    set_param(cs, 'AbsTol',        '1e-9');
    set_param(cs, 'MaxStep',       '0.5');
    set_param(cs, 'StartTime',     '0');
    set_param(cs, 'StopTime',      '1800');
    set_param(cs, 'SaveOutput',    'on');
    set_param(cs, 'SaveFormat',    'Dataset');
    set_param(cs, 'ReturnWorkspaceOutputs', 'on');

    % Auto-create Qext_ts from p.Qext_profile if user has not supplied one
    set_param(modelName, 'InitFcn', [...
        'if ~evalin(''base'',''exist(''''Qext_ts'''',''''var'''')''),' newline ...
        '  p_ = evalin(''base'',''p'');' newline ...
        '  Qext_ts = timeseries(p_.Qext_profile.value, p_.Qext_profile.time);' newline ...
        '  assignin(''base'',''Qext_ts'',Qext_ts);' newline ...
        'end']);

    % =================================================================
    %  Layout (approximate pixel coords):
    %     x ~  20.. 140  inputs
    %     x ~ 200.. 460  reaction MATLAB Function blocks (4 stacked)
    %     x ~ 500.. 600  Q_rxn sum + tSEI integrator
    %     x ~ 640.. 840  energy balance, T integrator
    %     x ~ 900..1000  outports
    % =================================================================

    add_block('built-in/Note', [modelName '/Title'], ...
        'Position', [20 5 400 25], ...
        'Text', 'Coupled ETP Model - Hatchard-Kim 0D (Phase 1)', ...
        'FontSize', '14', 'FontWeight', 'bold');

    % ----------------------------------------------------------------
    %  Q_ext source (From Workspace)
    % ----------------------------------------------------------------
    add_block('simulink/Sources/From Workspace', [modelName '/Qext_src'], ...
        'Position', [40 380 110 420], ...
        'VariableName', 'Qext_ts', ...
        'OutputAfterFinalValue', 'Holding final value');

    % ----------------------------------------------------------------
    %  Heat loss (single Fcn block - only depends on T)
    %  Constants embedded literally so the block doesn't need to
    %  re-evaluate the params struct at every step.
    % ----------------------------------------------------------------
    p_now = evalin('base','p');
    hAs  = p_now.cell.h_conv * p_now.cell.As;
    epsSigmaAs = p_now.cell.eps * p_now.sigma * p_now.cell.As;
    Tamb = p_now.T_amb;
    add_block('built-in/From', [modelName '/From_T_loss'], ...
        'Position', [40 460 80 480], 'GotoTag', 'Tcell');
    add_block('simulink/User-Defined Functions/Fcn', [modelName '/Q_loss_fcn'], ...
        'Position', [110 450 320 490], ...
        'Expr', sprintf('%.6e*(u-%.6e) + %.6e*(u^4-%.6e)', ...
                        hAs, Tamb, epsSigmaAs, Tamb^4));
    add_line(modelName, 'From_T_loss/1', 'Q_loss_fcn/1', 'autorouting','on');

    % ----------------------------------------------------------------
    %  Reaction subsystems
    % ----------------------------------------------------------------
    rxn_names = {'SEI','AnE','Cat','Elec'};
    rxn_y0    = [60, 160, 260, 360];

    for i = 1:4
        build_reaction_block(modelName, i, rxn_names{i}, rxn_y0(i));
    end

    % ----------------------------------------------------------------
    %  Sum the four Q_i contributions
    % ----------------------------------------------------------------
    add_block('built-in/Sum', [modelName '/Q_rxn_sum'], ...
        'Position', [500 100 530 240], ...
        'Inputs', '++++', 'IconShape', 'rectangular');
    for i = 1:4
        add_line(modelName, sprintf('Rxn_%s/2', rxn_names{i}), ...
                            sprintf('Q_rxn_sum/%d', i), 'autorouting','on');
    end

    add_block('built-in/Goto', [modelName '/Q_rxn_goto'], ...
        'Position', [555 160 615 180], ...
        'GotoTag', 'Q_rxn', 'TagVisibility', 'global');
    add_line(modelName, 'Q_rxn_sum/1', 'Q_rxn_goto/1');

    % ----------------------------------------------------------------
    %  Energy balance:  dT/dt = (Q_rxn + Q_ext - Q_loss) / (m*cp)
    % ----------------------------------------------------------------
    add_block('built-in/From', [modelName '/From_Q_rxn'], ...
        'Position', [620 290 660 310], 'GotoTag', 'Q_rxn');
    add_block('built-in/Sum', [modelName '/E_balance'], ...
        'Position', [680 280 710 340], ...
        'Inputs', '++-', 'IconShape', 'round');
    add_line(modelName, 'From_Q_rxn/1', 'E_balance/1', 'autorouting','on');
    add_line(modelName, 'Qext_src/1',   'E_balance/2', 'autorouting','on');
    add_line(modelName, 'Q_loss_fcn/1', 'E_balance/3', 'autorouting','on');

    add_block('built-in/Gain', [modelName '/inv_mcp'], ...
        'Position', [730 295 780 325], ...
        'Gain', sprintf('%.8e', 1/(p_now.cell.m*p_now.cell.cp)));
    add_line(modelName, 'E_balance/1', 'inv_mcp/1', 'autorouting','on');

    add_block('built-in/Integrator', [modelName '/T_integrator'], ...
        'Position', [800 295 840 325], ...
        'InitialCondition', sprintf('%.6f', p_now.cell.T0));
    add_line(modelName, 'inv_mcp/1', 'T_integrator/1', 'autorouting','on');

    add_block('built-in/Goto', [modelName '/Tcell_goto'], ...
        'Position', [860 300 920 320], ...
        'GotoTag', 'Tcell', 'TagVisibility', 'global');
    add_line(modelName, 'T_integrator/1', 'Tcell_goto/1', 'autorouting','on');

    % ----------------------------------------------------------------
    %  tSEI integrator (driven by AnE block, port 3 = dtSEI/dt)
    % ----------------------------------------------------------------
    add_block('built-in/Integrator', [modelName '/tSEI_integrator'], ...
        'Position', [500 200 540 230], ...
        'InitialCondition', sprintf('%.8e', p_now.rxn.tSEI_init));
    add_line(modelName, 'Rxn_AnE/3', 'tSEI_integrator/1', 'autorouting','on');
    add_block('built-in/Goto', [modelName '/tSEI_goto'], ...
        'Position', [560 205 620 225], ...
        'GotoTag', 'tSEI', 'TagVisibility', 'global');
    add_line(modelName, 'tSEI_integrator/1', 'tSEI_goto/1', 'autorouting','on');

    % ----------------------------------------------------------------
    %  Outports
    % ----------------------------------------------------------------
    add_outport_from(modelName, 'T_out',     1, 'Tcell',    [950 295 980 315]);
    add_outport_from(modelName, 'xSEI_out',  2, 'x_SEI',    [950  60 980  80]);
    add_outport_from(modelName, 'cne_out',   3, 'x_AnE',    [950 160 980 180]);
    add_outport_from(modelName, 'alpha_out', 4, 'x_Cat',    [950 260 980 280]);
    add_outport_from(modelName, 'ce_out',    5, 'x_Elec',   [950 360 980 380]);
    add_outport_from(modelName, 'tSEI_out',  6, 'tSEI',     [950 205 980 225]);
    add_outport_from(modelName, 'Qrxn_out',  7, 'Q_rxn',    [950 400 980 420]);

    % Pressure subsystem (Phase-3 stub)
    add_block('built-in/Subsystem', [modelName '/PressureModel'], ...
        'Position', [200 480 540 560]);
    build_pressure_subsystem([modelName '/PressureModel']);
    add_block('simulink/Sinks/Out1', [modelName '/P_out'], ...
        'Position', [950 510 980 530], 'Port', '8');
    add_line(modelName, 'PressureModel/1', 'P_out/1', 'autorouting','on');

    % --- save ----------------------------------------------------------
    saveDir = fileparts(mfilename('fullpath'));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    outPath = fullfile(saveDir, [modelName '.slx']);
    save_system(modelName, outPath);
    fprintf('build_etp_model: saved %s\n', outPath);
end

% =====================================================================
%   add_outport_from  -- helper to wire a Goto-tag into a numbered Outport
% =====================================================================
function add_outport_from(modelName, outName, port, gotoTag, pos)
    fromName = ['From_' outName];
    add_block('built-in/From', [modelName '/' fromName], ...
        'Position', [pos(1)-50 pos(2) pos(1)-10 pos(2)+20], ...
        'GotoTag', gotoTag);
    add_block('simulink/Sinks/Out1', [modelName '/' outName], ...
        'Position', pos, 'Port', sprintf('%d', port));
    add_line(modelName, [fromName '/1'], [outName '/1'], 'autorouting','on');
end

% =====================================================================
%   build_reaction_block  -- MATLAB Function block per reaction
% =====================================================================
function build_reaction_block(modelName, idx, label, yTop)
% Creates a subsystem 'Rxn_<label>' with:
%   Inputs : T (port 1), tSEI (port 2)
%   Outputs: 1=x_i  2=Q_i   [3=dtSEI/dt -- only for AnE]
%
% Inside: a MATLAB Function block computes the rate, and an Integrator
% holds the reaction extent.

    sysPath = [modelName '/Rxn_' label];
    add_block('built-in/Subsystem', sysPath, ...
        'Position', [200 yTop 460 yTop+80]);

    % --- clear the default In1->Out1 wiring ---------------------------
    inner_lines = find_system(sysPath, 'SearchDepth', 1, 'FindAll','on', 'type','line');
    for k = 1:numel(inner_lines)
        try, delete_line(inner_lines(k)); catch, end
    end
    inner = find_system(sysPath, 'SearchDepth', 1, 'LookUnderMasks','all');
    for k = 1:numel(inner)
        if ~strcmp(inner{k}, sysPath)
            try, delete_block(inner{k}); catch, end
        end
    end

    % --- Inports ------------------------------------------------------
    add_block('simulink/Sources/In1', [sysPath '/T_in'], ...
        'Position', [20 30 50 50], 'Port', '1');
    add_block('simulink/Sources/In1', [sysPath '/tSEI_in'], ...
        'Position', [20 70 50 90], 'Port', '2');

    % --- MATLAB Function block ----------------------------------------
    add_block('simulink/User-Defined Functions/MATLAB Function', ...
        [sysPath '/rate_fcn'], 'Position', [110 20 220 100]);
    set_rate_function_code(sysPath, label);
    % Give Stateflow a moment to re-parse the chart so the I/O port
    % count matches the function signature we just wrote.  Without this
    % the subsequent add_line('rate_fcn/3', ...) for AnE can fail
    % because the chart still has its default single-input signature.
    pause(0.2);

    % --- Reaction-extent integrator -----------------------------------
    p_now = evalin('base','p');
    add_block('built-in/Integrator', [sysPath '/x_int'], ...
        'Position', [260 30 300 60], ...
        'InitialCondition', sprintf('%.8e', p_now.rxn.x0(idx)));

    % --- Wiring -------------------------------------------------------
    add_line(sysPath, 'T_in/1',     'rate_fcn/1', 'autorouting','on');
    add_line(sysPath, 'x_int/1',    'rate_fcn/2', 'autorouting','on');
    add_line(sysPath, 'tSEI_in/1',  'rate_fcn/3', 'autorouting','on');

    add_line(sysPath, 'rate_fcn/1', 'x_int/1',    'autorouting','on');

    % --- Outports -----------------------------------------------------
    add_block('simulink/Sinks/Out1', [sysPath '/x_out'], ...
        'Position', [330 30 360 50], 'Port', '1');
    add_line(sysPath, 'x_int/1', 'x_out/1', 'autorouting','on');

    add_block('simulink/Sinks/Out1', [sysPath '/Q_out'], ...
        'Position', [330 60 360 80], 'Port', '2');
    add_line(sysPath, 'rate_fcn/2', 'Q_out/1', 'autorouting','on');

    if strcmp(label,'AnE')
        add_block('simulink/Sinks/Out1', [sysPath '/dtSEI_out'], ...
            'Position', [330 90 360 110], 'Port', '3');
        add_line(sysPath, 'rate_fcn/3', 'dtSEI_out/1', 'autorouting','on');
    end

    % --- outside the subsystem: feed T and tSEI in --------------------
    add_block('built-in/From', [modelName '/From_T_' label], ...
        'Position', [160 yTop+20 200 yTop+40], 'GotoTag', 'Tcell');
    add_line(modelName, ['From_T_' label '/1'], ['Rxn_' label '/1'], ...
        'autorouting','on');

    if strcmp(label,'AnE')
        add_block('built-in/From', [modelName '/From_tSEI_AnE'], ...
            'Position', [160 yTop+50 200 yTop+70], 'GotoTag', 'tSEI');
        add_line(modelName, 'From_tSEI_AnE/1', 'Rxn_AnE/2', 'autorouting','on');
    else
        % Other reactions ignore tSEI -- ground the port
        add_block('simulink/Sources/Ground', [modelName '/Gnd_tSEI_' label], ...
            'Position', [160 yTop+50 200 yTop+70]);
        add_line(modelName, ['Gnd_tSEI_' label '/1'], ['Rxn_' label '/2'], ...
            'autorouting','on');
    end

    % --- Publish the reaction extent for downstream use ---------------
    add_block('built-in/Goto', [modelName '/x_' label '_goto'], ...
        'Position', [475 yTop+10 545 yTop+30], ...
        'GotoTag', ['x_' label], 'TagVisibility', 'global');
    add_line(modelName, ['Rxn_' label '/1'], ['x_' label '_goto/1'], ...
        'autorouting','on');
end

% =====================================================================
%   set_rate_function_code  -- write the MATLAB code for each reaction
%
%   Parameter values are embedded as numeric literals at build time.
%   This avoids evalin('base',...) inside the chart (slow + not codegen
%   compatible) and means the model rebuild has to be re-run whenever
%   params change.  That's a feature, not a bug: it makes parameter
%   provenance auditable from the .slx alone.
% =====================================================================
function set_rate_function_code(sysPath, label)
    fcnPath = [sysPath '/rate_fcn'];
    sfRoot = sfroot;
    chart = sfRoot.find('-isa','Stateflow.EMChart','Path',fcnPath);
    if isempty(chart)
        pause(0.2);
        chart = sfRoot.find('-isa','Stateflow.EMChart','Path',fcnPath);
    end
    assert(~isempty(chart), 'Could not locate MATLAB Function chart at %s', fcnPath);

    % Pull current parameter values from base workspace
    p = evalin('base','p');
    Ru = p.Ru;
    Vc = p.cell.V;

    switch label
        case 'SEI'
            i = 1;
            code = sprintf([ ...
                'function [dxdt, Q] = rate_fcn(T, x, tSEI_unused)\n' ...
                '%%#codegen\n' ...
                'A    = %.8e;\n' ...
                'Ea   = %.8e;\n' ...
                'dH   = %.8e;\n' ...
                'w    = %.8e;\n' ...
                'V    = %.8e;\n' ...
                'Ru   = %.8e;\n' ...
                'c    = %d;\n' ...
                'R = A * exp(-Ea/(Ru*T)) * max(x,0);\n' ...
                'dxdt = c * R;\n' ...
                'Q    = w * V * dH * R;\n' ...
                'end\n'], ...
                p.rxn.A(i), p.rxn.Ea(i), p.rxn.dH(i), p.rxn.w(i), Vc, Ru, p.rxn.c(i));
        case 'AnE'
            i = 2;
            code = sprintf([ ...
                'function [dxdt, Q, dtSEI] = rate_fcn(T, x, tSEI)\n' ...
                '%%#codegen\n' ...
                'A     = %.8e;\n' ...
                'Ea    = %.8e;\n' ...
                'dH    = %.8e;\n' ...
                'w     = %.8e;\n' ...
                'V     = %.8e;\n' ...
                'Ru    = %.8e;\n' ...
                'c     = %d;\n' ...
                'tSEI0 = %.8e;\n' ...
                'g     = exp(-tSEI/tSEI0);\n' ...
                'R     = A * exp(-Ea/(Ru*T)) * max(x,0) * g;\n' ...
                'dxdt  = c * R;\n' ...
                'Q     = w * V * dH * R;\n' ...
                'dtSEI = R;\n' ...
                'end\n'], ...
                p.rxn.A(i), p.rxn.Ea(i), p.rxn.dH(i), p.rxn.w(i), Vc, Ru, ...
                p.rxn.c(i), p.rxn.tSEI0);
        case 'Cat'
            i = 3;
            code = sprintf([ ...
                'function [dxdt, Q] = rate_fcn(T, x, tSEI_unused)\n' ...
                '%%#codegen\n' ...
                'A    = %.8e;\n' ...
                'Ea   = %.8e;\n' ...
                'dH   = %.8e;\n' ...
                'w    = %.8e;\n' ...
                'V    = %.8e;\n' ...
                'Ru   = %.8e;\n' ...
                'c    = %d;\n' ...
                'R = A * exp(-Ea/(Ru*T)) * max(x,0) * max(1-x,0);\n' ...
                'dxdt = c * R;\n' ...
                'Q    = w * V * dH * R;\n' ...
                'end\n'], ...
                p.rxn.A(i), p.rxn.Ea(i), p.rxn.dH(i), p.rxn.w(i), Vc, Ru, p.rxn.c(i));
        case 'Elec'
            i = 4;
            code = sprintf([ ...
                'function [dxdt, Q] = rate_fcn(T, x, tSEI_unused)\n' ...
                '%%#codegen\n' ...
                'A    = %.8e;\n' ...
                'Ea   = %.8e;\n' ...
                'dH   = %.8e;\n' ...
                'w    = %.8e;\n' ...
                'V    = %.8e;\n' ...
                'Ru   = %.8e;\n' ...
                'c    = %d;\n' ...
                'R = A * exp(-Ea/(Ru*T)) * max(x,0);\n' ...
                'dxdt = c * R;\n' ...
                'Q    = w * V * dH * R;\n' ...
                'end\n'], ...
                p.rxn.A(i), p.rxn.Ea(i), p.rxn.dH(i), p.rxn.w(i), Vc, Ru, p.rxn.c(i));
    end
    chart.Script = code;
end

% =====================================================================
%   PRESSURE SUBSYSTEM (Phase 3 stub)
% =====================================================================
function build_pressure_subsystem(sysPath)
    % Clear default blocks
    inner_lines = find_system(sysPath, 'SearchDepth', 1, 'FindAll','on', 'type','line');
    for k = 1:numel(inner_lines)
        try, delete_line(inner_lines(k)); catch, end
    end
    inner = find_system(sysPath, 'SearchDepth', 1, 'LookUnderMasks','all');
    for k = 1:numel(inner)
        if ~strcmp(inner{k}, sysPath)
            try, delete_block(inner{k}); catch, end
        end
    end

    p_now = evalin('base','p');
    add_block('simulink/Sources/Constant', [sysPath '/P_stub'], ...
        'Position', [40 30 90 60], 'Value', sprintf('%.6e', p_now.pressure.P0));
    add_block('simulink/Sinks/Out1', [sysPath '/P_out'], ...
        'Position', [120 30 150 60], 'Port', '1');
    add_line(sysPath, 'P_stub/1', 'P_out/1', 'autorouting','on');

    add_block('built-in/Note', [sysPath '/PhaseNote'], ...
        'Position', [40 80 280 100], ...
        'Text', 'PRESSURE - Phase 3 stub.  See build_etp_model.m header.', ...
        'FontSize', '10', 'FontAngle', 'italic');
end
