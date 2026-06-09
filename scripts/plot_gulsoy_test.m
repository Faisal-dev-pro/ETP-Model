function fig = plot_gulsoy_test(d, varargin)
% PLOT_GULSOY_TEST  Visualise one parsed Gulsoy test for sanity-checking.
%
%   fig = plot_gulsoy_test(data(1));
%   fig = plot_gulsoy_test(data(1), 'savePath', 'results/test1.png');
%
% Five panels:
%   1. Temperatures (internal, surface mid/pos/neg, vent 5/10 mm)
%   2. Internal pressure
%   3. Cell voltage
%   4. Heater power profile (if populated)
%   5. dT_internal / dt — useful for spotting onset
%
% Author: <your name>, 2026. License: MIT.

    opts.savePath = '';
    opts.title    = '';
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end
    if isempty(opts.title)
        opts.title = sprintf('Gulsoy test: %s', d.cell_id);
    end

    fig = figure('Color','w','Position',[100 80 1000 900]);
    tlo = tiledlayout(fig, 5, 1, 'TileSpacing','compact','Padding','compact');
    title(tlo, opts.title, 'FontWeight', 'bold', 'Interpreter', 'none');

    % --- 1. temperatures (low-rate stream) ---
    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
    plotIfNotEmpty(ax1, d.t, d.T_internal,    'k-',  1.8, 'T_{internal}');
    plotIfNotEmpty(ax1, d.t, d.T_surface_mid, 'b-',  1.2, 'T_{surf mid}');
    plotIfNotEmpty(ax1, d.t, d.T_surface_pos, 'r-',  1.0, 'T_{surf +}');
    plotIfNotEmpty(ax1, d.t, d.T_surface_neg, 'g-',  1.0, 'T_{surf -}');
    plotIfNotEmpty(ax1, d.t, d.T_vent_5mm,    'm--', 1.0, 'T_{vent 5mm}');
    plotIfNotEmpty(ax1, d.t, d.T_vent_10mm,   'c--', 1.0, 'T_{vent 10mm}');
    ylabel(ax1, 'T [degC]'); legend(ax1, 'Location','northwest','NumColumns',2);
    title(ax1, 'Temperatures');

    % --- 2. internal pressure (high-rate stream) ---
    ax2 = nexttile; grid(ax2,'on');
    if ~isempty(d.t_fast) && ~isempty(d.P_internal)
        % decimate to keep the plot snappy
        step = max(1, floor(numel(d.t_fast)/50000));
        plot(ax2, d.t_fast(1:step:end), d.P_internal(1:step:end), 'b-', 'LineWidth', 1.2);
    end
    ylabel(ax2, 'P_{internal} [bar]'); title(ax2, 'Internal gas pressure');

    % --- 3. cell voltage (high-rate stream) ---
    ax3 = nexttile; grid(ax3,'on');
    if ~isempty(d.t_fast) && ~isempty(d.V_cell)
        step = max(1, floor(numel(d.t_fast)/50000));
        plot(ax3, d.t_fast(1:step:end), d.V_cell(1:step:end), 'r-', 'LineWidth', 1.2);
    end
    ylabel(ax3, 'V_{cell} [V]'); title(ax3, 'Cell voltage');

    % --- 4. heater profile ---
    ax4 = nexttile; grid(ax4,'on');
    if ~isempty(d.heater_power)
        plot(ax4, d.t, d.heater_power, 'k-', 'LineWidth', 1.5);
        ylabel(ax4, 'Q_{ext} [W]');
    else
        text(ax4, 0.5, 0.5, '(heater_power not configured)', ...
            'HorizontalAlignment','center', 'Units','normalized', ...
            'Interpreter','none');
    end
    title(ax4, 'Heater profile');

    % --- 5. dT/dt ---
    ax5 = nexttile; grid(ax5,'on');
    if ~isempty(d.t) && ~isempty(d.T_internal)
        dT = gradient(d.T_internal, d.t);
        plot(ax5, d.t, dT, 'k-', 'LineWidth', 1.0); hold(ax5,'on');
        yline(ax5, 1, 'r:', '1 K/s onset threshold');
        ylabel(ax5, 'dT_{int}/dt [K/s]');
        set(ax5, 'YScale','log');
        ylim(ax5, [1e-2 1e3]);
    end
    title(ax5, 'Internal temperature rate (log)');
    xlabel(ax5, 't [s]');

    linkaxes([ax1 ax2 ax3 ax4 ax5], 'x');

    if ~isempty(opts.savePath)
        outDir = fileparts(opts.savePath);
        if ~isempty(outDir) && ~exist(outDir,'dir'), mkdir(outDir); end
        exportgraphics(fig, opts.savePath, 'Resolution', 130);
        fprintf('Saved plot: %s\n', opts.savePath);
    end
end

function plotIfNotEmpty(ax, t, y, style, lw, lbl)
    if isempty(t) || isempty(y), return; end
    plot(ax, t, y, style, 'LineWidth', lw, 'DisplayName', lbl);
end
