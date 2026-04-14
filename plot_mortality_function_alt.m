%% Decreasing disease-induced mortality function f(I) — shape across values of a
%
%   f(I)  = f0 - (f0 - f1)/(1 + exp(-(I-theta)/s))
%
%     f0 = 0.02    (baseline at zero prevalence)
%     f1 = 0.0002  (floor at high prevalence)
%     s = 0.01

%% ── Style ────────────────────────────────────────────────────────────────

%addpath('/Users/denispatterson/Documents/MatCont7p4/');
set_ggplot_style();
gg = evalin('base', 'gg_colours');

LW = 1.2;
FN = 'Helvetica';
FS = 10;

%% ── Fixed parameters ─────────────────────────────────────────────────────

f0 = 0.02;    % f(I) at I = 0   (baseline mortality)
f1 = 0.0002;  % f(I) as I → ∞  (floor mortality)
s  = 0.0005;     % logistic scale  (steepness of transition)

%% ── Values of theta to compare ───────────────────────────────────────────
%
%   theta is the inflection point of the logistic; five values spanning the
%   I ∈ [0, 0.5] plotting range, one per gg colour.

theta_vals = [0.001, 0.005, 0.010, 0.015];

%% ── I grid ───────────────────────────────────────────────────────────────
%
%   [0, 0.015] spans the range of theta values so every curve shows its
%   inflection point and a visible portion of both plateau regions.

I_vec = linspace(0, 0.015, 500);

%% ── Colours ──────────────────────────────────────────────────────────────
%
%   Use gg_colours palette (one colour per theta value):
%     theta = 0.001  →  pink/salmon  gg(1,:)
%     theta = 0.005  →  olive-green  gg(2,:)
%     theta = 0.010  →  teal-cyan    gg(3,:)
%     theta = 0.015  →  purple       gg(4,:)

n_curves = length(theta_vals);
cmap     = gg(1:n_curves, :);

%% ── Figure ───────────────────────────────────────────────────────────────

fig = figure('Units', 'inches', 'Position', [1 1 5 4], 'Color', 'w');
ax  = axes(fig);
hold(ax, 'on');
ax.TickDir = 'out';   % explicit override — belt-and-braces

for k = 1:n_curves
    thetak = theta_vals(k);
    fk     = f0 - (f0 - f1) ./ (1 + exp(-(I_vec - thetak) / s));

    plot(ax, I_vec, fk, '-', ...
        'Color',       cmap(k, :), ...
        'LineWidth',   LW, ...
        'DisplayName', sprintf('\\theta = %.3f', thetak));
end

%% ── Axis limits: ggplot2-style 5% left-expansion ─────────────────────────
%
%   x-data lives on [0, 0.25]; pad left by 5% so axis starts just below 0
%   while the first tick label remains at I = 0.

I_max   = I_vec(end);
x_pad   = 0.05 * I_max;
xlim(ax, [-x_pad/2, I_max + x_pad]);
xticks(ax, 0 : 0.005 : I_max);   % ticks at 0, 0.005, 0.010, 0.015
ax.XAxis.Exponent = 0;           % force decimal display (suppress ×10⁻³)

%% Reference lines ─────────────────────────────────────────────────────────

% yline(ax, b, '--', ...
%     'Color',             [0.5 0.5 0.5], ...
%     'LineWidth',         0.8, ...
%     'Label',             '$f(0) = b$', ...
%     'Interpreter',       'latex', ...
%     'LabelHorizontalAlignment', 'left', ...
%     'HandleVisibility',  'off');
% 
% yline(ax, b/m, ':', ...
%     'Color',             [0.5 0.5 0.5], ...
%     'LineWidth',         0.8, ...
%     'Label',             '$f(\infty) = b/m$', ...
%     'Interpreter',       'latex', ...
%     'LabelHorizontalAlignment', 'left', ...
%     'HandleVisibility',  'off');

%% Labels and legend ───────────────────────────────────────────────────────

xlabel(ax, 'Infection levels (I)',                  'Interpreter', 'tex');
ylabel(ax, 'Per-infection mortality rate (f(I))', 'Interpreter', 'tex');

legend(ax, ...
    'Interpreter', 'tex', ...
    'Location',    'east', ...
    'FontSize',    FS, ...
    'Box',         'off');

ylim(ax, [0, f0 * 1.1]);   % small headroom above the baseline

% Suppress first y-tick mark and its label (ggplot2 style)
yt = ax.YTick;
if numel(yt) > 1,  ax.YTick = yt(2:end);  end

%% ── Save ─────────────────────────────────────────────────────────────────

% set(fig, 'PaperPositionMode', 'auto');
% 
% if ~exist('plots', 'dir')
%     mkdir('plots');
% end
% 
% exportgraphics(fig, 'plots/mortality_function.pdf', 'ContentType', 'vector');
% fprintf('Saved → plots/mortality_function.pdf\n');
