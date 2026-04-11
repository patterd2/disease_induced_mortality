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

addpath('/Users/denispatterson/Documents/MatCont7p4/');
set_ggplot_style();
gg = evalin('base', 'gg_colours');

%% ── Parameters ───────────────────────────────────────────────────────────

Lambda  = 0.0003846154;
R0_val  = 0.9999999999901998;
mu      = 0.0003846154;
delta   = 0.07692308;
gamma_r = 1;           % 'gamma_r' avoids shadowing MATLAB's built-in gamma()
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
tspan = [0, 100];

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

fig = figure('Units', 'inches', 'Position', [1 1 9 7], 'Color', 'w');

LW = 1.2;   % line width (set_ggplot_style default is 1.2)
FS = 10;    % font size  (set_ggplot_style default is 10)
FN = 'Helvetica';

%% Panel A — Phase portrait in (S, I) plane ───────────────────────────────

ax1 = subplot(2, 2, 1);
hold(ax1, 'on');

plot(ax1, S, I, '-', 'Color', col_SI, 'LineWidth', LW);

% IC marker (filled circle)
plot(ax1, S(1), I(1), 'o', ...
    'Color', mk_col, 'MarkerFaceColor', col_SI, ...
    'MarkerSize', 6, 'LineWidth', 0.8);

% End-state marker (filled square)
plot(ax1, S(end), I(end), 's', ...
    'Color', mk_col, 'MarkerFaceColor', col_SI, ...
    'MarkerSize', 6, 'LineWidth', 0.8);

xlabel(ax1, 'Susceptible,  $S$', 'Interpreter', 'latex');
ylabel(ax1, 'Infected,  $I$',    'Interpreter', 'latex');

% Panel letter A — placed after axes are rendered
xl1 = ax1.XLim;  yl1 = ax1.YLim;
text(ax1, xl1(1), yl1(2), 'A', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Panel B — S(t) ─────────────────────────────────────────────────────────

ax2 = subplot(2, 2, 2);
hold(ax2, 'on');

plot(ax2, t, S, '-', 'Color', col_S, 'LineWidth', LW);

xlabel(ax2, 'Time,  $t$',       'Interpreter', 'latex');
ylabel(ax2, 'Susceptible,  $S(t)$', 'Interpreter', 'latex');

xl2 = ax2.XLim;  yl2 = ax2.YLim;
text(ax2, xl2(1), yl2(2), 'B', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Panel C — I(t) ─────────────────────────────────────────────────────────

ax3 = subplot(2, 2, 3);
hold(ax3, 'on');

plot(ax3, t, I, '-', 'Color', col_SI, 'LineWidth', LW);

xlabel(ax3, 'Time,  $t$',     'Interpreter', 'latex');
ylabel(ax3, 'Infected,  $I(t)$', 'Interpreter', 'latex');

xl3 = ax3.XLim;  yl3 = ax3.YLim;
text(ax3, xl3(1), yl3(2), 'C', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Panel D — R(t) ─────────────────────────────────────────────────────────

ax4 = subplot(2, 2, 4);
hold(ax4, 'on');

plot(ax4, t, R, '-', 'Color', col_R, 'LineWidth', LW);

xlabel(ax4, 'Time,  $t$',        'Interpreter', 'latex');
ylabel(ax4, 'Recovered,  $R(t)$', 'Interpreter', 'latex');

xl4 = ax4.XLim;  yl4 = ax4.YLim;
text(ax4, xl4(1), yl4(2), 'D', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% ── Save ─────────────────────────────────────────────────────────────────

set(fig, 'PaperPositionMode', 'auto');

if ~exist('plots', 'dir')
    mkdir('plots');
end

exportgraphics(fig, 'plots/single_trajectory.pdf', 'ContentType', 'vector');
fprintf('Saved → plots/single_trajectory.pdf\n');
