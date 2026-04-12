%% SIR model with disease-induced mortality: single-trajectory analysis
%
%   ODE system
%   ----------
%   d    = -a*m / ((m-1)*b)
%   f(I) =  b + a*I^2 / (1 + d*I^2)          % I-dependent mortality rate
%
%   S'  = Lambda - beta*S*I - mu*S + delta*R
%   I'  = beta*S*I - (gamma_r+mu)*I - f(I)*I
%   R'  = gamma_r*I - (mu+delta)*R
%
%   where  beta = R0*(mu/Lambda)*(gamma_r+mu+b)
%
%   At the disease-free equilibrium (DFE): (S*, I*, R*) = (1, 0, 0)
%   since Lambda = mu.
%
%   Outputs: plots/single_trajectory.pdf

%% ── Style ────────────────────────────────────────────────────────────────

%addpath('/Users/denispatterson/Documents/MatCont7p4/');
set_ggplot_style();
gg = evalin('base', 'gg_colours');

%% ── Parameters ───────────────────────────────────────────────────────────

Lambda  = 0.0003846154;
R0_val  = 1.005;
mu      = 0.0003846154;
delta   = 0.07692308;
gamma_r = 1;  % 'gamma_r' avoids conflicts with the built-in gamma()
a       = -10000;
b       = 0.02;
m       = 100;

%% ── Derived quantities ───────────────────────────────────────────────────

% Mortality shape coefficient (always positive when a < 0)
d    = -a * m / ((m - 1) * b);

% Effective transmission coefficient
beta = R0_val * (mu / Lambda) * (gamma_r + mu + b);

% Disease-induced mortality rate as a function of I
%   f(0) = b  (baseline),  f(inf) -> b/m  (saturating, lower at high I)
f    = @(I) b + a .* I.^2 ./ (1 + d .* I.^2);

fprintf('\n── Derived parameters ───────────────────────────────────────\n');
fprintf('  d          = %12.4f\n', d);
fprintf('  beta       = %12.8f\n', beta);
fprintf('  f(I=0)     = %12.6f  (= b)\n', f(0));
fprintf('  f(I->inf)  = %12.6f  (= b/m)\n', b/m);
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% ── ODE right-hand side ──────────────────────────────────────────────────

odefun = @(t, u) [ ...
    Lambda - beta*u(1)*u(2) - mu*u(1)             + delta*u(3); ...
    beta*u(1)*u(2) - (gamma_r + mu)*u(2) - f(u(2))*u(2);       ...
    gamma_r*u(2)   - (mu + delta)*u(3)                          ];

%% ── Initial conditions and integration ───────────────────────────────────

u0   = [0.99; 0.01; 0];           % S(0), I(0), R(0) — near DFE with small I
tspan = [0, 1000];

odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', [1 2 3]);

[t, u] = ode45(odefun, tspan, u0, odeOpts);

S = u(:, 1);
I = u(:, 2);
R = u(:, 3);
N = S + I + R;   % total population (not conserved due to disease mortality)

fprintf('── Endpoint values at t = %g ─────────────────────────────────\n', tspan(end));
fprintf('  S = %.6f,  I = %.6f,  R = %.6f,  N = %.6f\n', S(end), I(end), R(end), N(end));
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% ── Colour assignments ───────────────────────────────────────────────────

col_SI  = gg(1, :);   % pink/salmon  — phase portrait & I(t)
col_S   = gg(2, :);   % olive-green  — S(t)
col_R   = gg(3, :);   % teal-cyan    — R(t)
mk_col  = [0.25 0.25 0.25];   % dark grey for IC / end-state markers

%% ── Figure ───────────────────────────────────────────────────────────────
%
%   Layout (2 × 3 grid):
%     Row 1:  A — 2D phase portrait (S,I)  |  B — 3D phase portrait (S,I,R) [double-wide]
%     Row 2:  C — S(t)  |  D — I(t)  |  E — R(t)

fig = figure('Units', 'inches', 'Position', [1 1 13 7], 'Color', 'w');

LW = 1.2;   % line width (set_ggplot_style default is 1.2)
FS = 10;    % font size  (set_ggplot_style default is 10)
FN = 'Helvetica';

col_3d = gg(4, :);   % purple — 3D phase portrait

%% Panel A — 2D phase portrait in (S, I) plane ────────────────────────────

ax1 = subplot(2, 3, 1);
hold(ax1, 'on');
ax1.TickDir = 'out';

plot(ax1, S, I, '-', 'Color', col_SI, 'LineWidth', LW);

% IC marker (filled circle)
plot(ax1, S(1), I(1), 'o', ...
    'Color', mk_col, 'MarkerFaceColor', col_SI, ...
    'MarkerSize', 6, 'LineWidth', 0.8);

% End-state marker (filled square)
plot(ax1, S(end), I(end), 's', ...
    'Color', mk_col, 'MarkerFaceColor', col_SI, ...
    'MarkerSize', 6, 'LineWidth', 0.8);

xlabel(ax1, 'Susceptible,  {\itS}', 'Interpreter', 'tex');
ylabel(ax1, 'Infected,  {\itI}',   'Interpreter', 'tex');

yt = ax1.YTick;
if numel(yt) > 1,  ax1.YTick = yt(2:end);  end

xl1 = ax1.XLim;  yl1 = ax1.YLim;
text(ax1, xl1(1), yl1(2), 'A', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Panel B — 3D phase portrait in (S, I, R) space ────────────────────────

ax2 = subplot(2, 3, [2, 3]);
hold(ax2, 'on');
ax2.TickDir = 'out';

plot3(ax2, S, I, R, '-', 'Color', col_3d, 'LineWidth', LW);

% IC marker (filled circle)
plot3(ax2, S(1), I(1), R(1), 'o', ...
    'Color', mk_col, 'MarkerFaceColor', col_3d, ...
    'MarkerSize', 6, 'LineWidth', 0.8);

% End-state marker (filled square)
plot3(ax2, S(end), I(end), R(end), 's', ...
    'Color', mk_col, 'MarkerFaceColor', col_3d, ...
    'MarkerSize', 6, 'LineWidth', 0.8);

xlabel(ax2, 'Susceptible,  {\itS}', 'Interpreter', 'tex');
ylabel(ax2, 'Infected,  {\itI}',   'Interpreter', 'tex');
zlabel(ax2, 'Recovered,  {\itR}',  'Interpreter', 'tex');

view(ax2, [-35, 22]);
ax2.Box = 'off';
grid(ax2, 'on');

yt = ax2.YTick;
if numel(yt) > 1,  ax2.YTick = yt(2:end);  end

xl2 = ax2.XLim;  yl2 = ax2.YLim;  zl2 = ax2.ZLim;
text(ax2, xl2(1), yl2(1), zl2(2), 'B', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Panel C — S(t) ─────────────────────────────────────────────────────────

ax3 = subplot(2, 3, 4);
hold(ax3, 'on');
ax3.TickDir = 'out';

plot(ax3, t, S, '-', 'Color', col_S, 'LineWidth', LW);

% ggplot2-style x-expansion: axis starts just below t=0, first tick at 0
t_pad = 0.05 * tspan(end);
xlim(ax3, [-t_pad, tspan(end) + t_pad]);
xticks(ax3, linspace(0, tspan(end), 5));

xlabel(ax3, 'Time,  {\itt}  (days)',         'Interpreter', 'tex');
ylabel(ax3, 'Susceptible,  {\itS}({\itt})', 'Interpreter', 'tex');

yt = ax3.YTick;
if numel(yt) > 1,  ax3.YTick = yt(2:end);  end

xl3 = ax3.XLim;  yl3 = ax3.YLim;
text(ax3, xl3(1), yl3(2), 'C', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Panel D — I(t) ─────────────────────────────────────────────────────────

ax4 = subplot(2, 3, 5);
hold(ax4, 'on');
ax4.TickDir = 'out';

plot(ax4, t, I, '-', 'Color', col_SI, 'LineWidth', LW);

xlim(ax4, [-t_pad, tspan(end) + t_pad]);
xticks(ax4, linspace(0, tspan(end), 5));

xlabel(ax4, 'Time,  {\itt}  (days)',       'Interpreter', 'tex');
ylabel(ax4, 'Infected,  {\itI}({\itt})', 'Interpreter', 'tex');

yt = ax4.YTick;
if numel(yt) > 1,  ax4.YTick = yt(2:end);  end

xl4 = ax4.XLim;  yl4 = ax4.YLim;
text(ax4, xl4(1), yl4(2), 'D', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Panel E — R(t) ─────────────────────────────────────────────────────────

ax5 = subplot(2, 3, 6);
hold(ax5, 'on');
ax5.TickDir = 'out';

plot(ax5, t, R, '-', 'Color', col_R, 'LineWidth', LW);

xlim(ax5, [-t_pad, tspan(end) + t_pad]);
xticks(ax5, linspace(0, tspan(end), 5));

xlabel(ax5, 'Time,  {\itt}  (days)',         'Interpreter', 'tex');
ylabel(ax5, 'Recovered,  {\itR}({\itt})', 'Interpreter', 'tex');

yt = ax5.YTick;
if numel(yt) > 1,  ax5.YTick = yt(2:end);  end

xl5 = ax5.XLim;  yl5 = ax5.YLim;
text(ax5, xl5(1), yl5(2), 'E', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% ── Save ─────────────────────────────────────────────────────────────────
%
% set(fig, 'PaperPositionMode', 'auto');
% 
% if ~exist('plots', 'dir')
%     mkdir('plots');
% end
% 
% exportgraphics(fig, 'plots/single_trajectory.pdf', 'ContentType', 'vector');
% fprintf('Saved → plots/single_trajectory.pdf\n');
