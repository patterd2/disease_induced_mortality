%% run_cycle_periods.m

%  Time unit: weeks (periods are converted to years for plotting, /52).

clear; close all; clc;

%% ---- shared / baseline parameters (manuscript Table 1 demography) ----
base = struct();
base.Lam   = 1/2600;     % recruitment (week^-1)
base.mu    = 1/2600;     % natural mortality (week^-1)  => S0 = Lam/mu = 1
base.b     = 0.02;       % baseline per-infection mortality f(0)
base.m     = 100;        % fold reduction at high infection: f(inf) = b/m
base.n     = 2;          % Hill coefficient (n=2 gives the richer dynamics)
base.delta = 0;          % (set per model below)

% common options: tighten/relax horizons & tolerances here if needed
opts          = struct();
opts.verbose  = true;
opts.relamp_min = 1e-3;  % below this relative amplitude = "no cycle" (NaN)

% =====================================================================
sis = base; sis.a = -100; sis.gamma = 0.5;
o1  = opts; o1.y0 = [0.985 0.012]; o1.Tbase = 6e4;
R0_sis = linspace(1.2, 1.55, 100);
[~, sisR0] = cycle_period_sweep('sis','R0', R0_sis, sis, o1);

%% (2) SIS:  period vs gamma
sis2 = base; sis2.a = -100; sis2.R0 = 1.15;
o2   = opts; o2.y0 = [0.985 0.012]; o2.Tbase = 6e4;
g_sis = linspace(0.111, 0.35, 50);
[~, sisG] = cycle_period_sweep('sis','gamma', g_sis, sis2, o2);

%% =====================================================================
%  (3) SIRS: period vs R0      (a = -100, gamma = 0.25, delta = 0.25)
%  Cycle region is bracketed by two Hopf points (R0 ~ 1.30 and ~1.55):
%  the period traces a U-shape. Seed mid-region.
% =====================================================================
sir = base; sir.a = -100; sir.gamma = 0.25; sir.delta = 0.25;
o3  = opts; o3.y0 = [0.70 0.02 0.10]; o3.Tbase = 4e4;
R0_sir = linspace(1.3, 1.45, 100);
[~, sirR0] = cycle_period_sweep('sirs','R0', R0_sir, sir, o3);

%% (4) SIRS: period vs delta   (a = -100, gamma = 0.25, R0 = 1.42)
sir2 = base; sir2.a = -100; sir2.gamma = 0.25; sir2.R0 = 1.2;
o4   = opts; o4.y0 = [0.70 0.02 0.10]; o4.Tbase = 4e4;
d_sir = linspace(0.15, 1, 100);
[~, sirD] = cycle_period_sweep('sirs','delta', d_sir, sir2, o4);

%% ----------------------------- plot -----------------------------
% ggplot2-style aesthetic shared with the rest of the repo (see
% set_ggplot_style / single_trajectory_SIR). Rate parameters are shown on the
% x-axis as their reciprocals (1/gamma = avg. infectious period, 1/delta =
% avg. immune period, in weeks), matching the convention used elsewhere.
set_ggplot_style();
gg = evalin('base', 'gg_colours');
LW = 1.4;  FN = 'Helvetica';
col_sis = gg(1,:);   % pink/salmon  — SIS panels
col_sir = gg(4,:);   % purple       — SIR(S) panels
ylab = 'limit-cycle period (years)';

fig = figure('Units','inches','Position',[1 1 9 7],'Color','w');

% (A) SIS: period vs R0
ax1 = subplot(2,2,1); hold(ax1,'on'); ax1.TickDir = 'out';
plot(ax1, sisR0.param_vals, sisR0.periods_yr, 'o-', 'Color', col_sis, ...
     'MarkerFaceColor', col_sis, 'MarkerSize', 4, 'LineWidth', LW);
xlabel(ax1, '\itR\rm_0', 'Interpreter','tex');
ylabel(ax1, ylab, 'Interpreter','tex');
title(ax1, 'SIS (\gamma = 0.5,  a = -100)', 'Interpreter','tex', 'FontWeight','normal');
panel_label(ax1, 'A', FN);

% (B) SIS: period vs 1/gamma (avg. infectious period)
ax2 = subplot(2,2,2); hold(ax2,'on'); ax2.TickDir = 'out';
plot(ax2, 1./sisG.param_vals, sisG.periods_yr, 'o-', 'Color', col_sis, ...
     'MarkerFaceColor', col_sis, 'MarkerSize', 4, 'LineWidth', LW);
xlabel(ax2, 'avg. infectious period (1/\gamma, weeks)', 'Interpreter','tex');
ylabel(ax2, ylab, 'Interpreter','tex');
title(ax2, 'SIS (\itR\rm_0 = 1.15,  a = -100)', 'Interpreter','tex', 'FontWeight','normal');
panel_label(ax2, 'B', FN);

% (C) SIR(S): period vs R0
ax3 = subplot(2,2,3); hold(ax3,'on'); ax3.TickDir = 'out';
plot(ax3, sirR0.param_vals, sirR0.periods_yr, 'o-', 'Color', col_sir, ...
     'MarkerFaceColor', col_sir, 'MarkerSize', 4, 'LineWidth', LW);
xlabel(ax3, '\itR\rm_0', 'Interpreter','tex');
ylabel(ax3, ylab, 'Interpreter','tex');
title(ax3, 'SIR(S) (\gamma = 0.25, \delta = 0.25, a = -100)', 'Interpreter','tex', 'FontWeight','normal');
panel_label(ax3, 'C', FN);

% (D) SIR(S): period vs 1/delta (avg. immune period)
ax4 = subplot(2,2,4); hold(ax4,'on'); ax4.TickDir = 'out';
plot(ax4, 1./sirD.param_vals, sirD.periods_yr, 'o-', 'Color', col_sir, ...
     'MarkerFaceColor', col_sir, 'MarkerSize', 4, 'LineWidth', LW);
xlabel(ax4, 'avg. immune period (1/\delta, weeks)', 'Interpreter','tex');
ylabel(ax4, ylab, 'Interpreter','tex');
title(ax4, 'SIR(S) (\itR\rm_0 = 1.2, \gamma = 0.25, a = -100)', 'Interpreter','tex', 'FontWeight','normal');
panel_label(ax4, 'D', FN);

% save vector + raster copies for the manuscript / SI
set(fig,'PaperPositionMode','auto');
print(fig,'cycle_periods','-dpng','-r200');
try, print(fig,'cycle_periods','-dpdf','-bestfit'); catch, end

fprintf('\nSaved cycle_periods.png (and .pdf if supported).\n');

% ---- local helper: bold panel letter in the top-left margin of an axes ----
function panel_label(ax, s, FN)
xl = ax.XLim; yl = ax.YLim;
text(ax, xl(1)-0.10*range(xl), yl(2), s, 'FontName', FN, 'FontSize', 12, ...
     'FontWeight','bold', 'VerticalAlignment','bottom', ...
     'HorizontalAlignment','left', 'Clipping','off');
end
