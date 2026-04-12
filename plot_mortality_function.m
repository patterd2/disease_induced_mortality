%% Decreasing disease-induced mortality function f(I) — shape across values of a
%
%   f(I)  = b + a*I^2 / (1 + d*I^2)
%   d(a)  = -a * m / ((m-1) * b)
%
%   For a < 0:
%     f(0)   = b      = 0.02    (baseline at zero prevalence)
%     f(inf) = b/m    = 0.0002  (floor at high prevalence)
%
%   Larger |a|  →  sharper drop at lower I  (inflection at I* = 1/sqrt(d))
%
%   Mirrors Figure 2A of the companion paper but for the decreasing case.
%   Output: plots/mortality_function.pdf

%% ── Style ────────────────────────────────────────────────────────────────

%addpath('/Users/denispatterson/Documents/MatCont7p4/');
set_ggplot_style();

LW = 1.2;
FN = 'Helvetica';
FS = 10;

%% ── Fixed parameters ─────────────────────────────────────────────────────

b = 0.02;    % baseline disease-induced mortality  (= f at I = 0)
m = 100;     % saturation parameter                (floor = b/m at I → ∞)

%% ── Values of a to compare ───────────────────────────────────────────────
%
%   5 values log-spaced from -10 to -10 000 (3 intermediate values)

a_vals = [-10, -100, -1000, -5000, -10000];

%   Precompute d and inflection I* for reference
fprintf('\n── Derived d and inflection I* ──────────────────────────────\n');
for k = 1:length(a_vals)
    dk = -a_vals(k) * m / ((m-1) * b);
    fprintf('  a = %7d  →  d = %9.1f,  I* = %.4f\n', ...
            a_vals(k), dk, 1/sqrt(dk));
end
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% ── I grid ───────────────────────────────────────────────────────────────
%
%   [0, 0.5] captures the full descent for all a values while keeping the
%   a = -10 curve visually separated from the others.

I_vec = linspace(0, 0.5, 2000);

%% ── Colours ──────────────────────────────────────────────────────────────
%
%   Sequential parula gradient: lightest = least negative a (slowest drop),
%   darkest = most negative a (sharpest drop near I = 0).

n_curves = length(a_vals);
cmap     = parula(n_curves);

%% ── Figure ───────────────────────────────────────────────────────────────

fig = figure('Units', 'inches', 'Position', [1 1 5 4], 'Color', 'w');
ax  = axes(fig);
hold(ax, 'on');

for k = 1:n_curves
    ak = a_vals(k);
    dk = -ak * m / ((m - 1) * b);
    fk = b + ak .* I_vec.^2 ./ (1 + dk .* I_vec.^2);

    plot(ax, I_vec, fk, '-', ...
        'Color',       cmap(k, :), ...
        'LineWidth',   LW, ...
        'DisplayName', sprintf('$a = %d$', ak));
end

%% Reference lines ─────────────────────────────────────────────────────────

yline(ax, b, '--', ...
    'Color',             [0.5 0.5 0.5], ...
    'LineWidth',         0.8, ...
    'Label',             '$f(0) = b$', ...
    'Interpreter',       'latex', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility',  'off');

yline(ax, b/m, ':', ...
    'Color',             [0.5 0.5 0.5], ...
    'LineWidth',         0.8, ...
    'Label',             '$f(\infty) = b/m$', ...
    'Interpreter',       'latex', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility',  'off');

%% Labels and legend ───────────────────────────────────────────────────────

xlabel(ax, 'Infected,  $I$',           'Interpreter', 'latex');
ylabel(ax, 'Mortality rate,  $f(I)$',  'Interpreter', 'latex');

legend(ax, ...
    'Interpreter', 'latex', ...
    'Location',    'northeast', ...
    'FontSize',    FS, ...
    'Box',         'off');

ylim(ax, [0, b * 1.15]);   % small headroom above the baseline

%% ── Save ─────────────────────────────────────────────────────────────────

set(fig, 'PaperPositionMode', 'auto');

if ~exist('plots', 'dir')
    mkdir('plots');
end

exportgraphics(fig, 'plots/mortality_function.pdf', 'ContentType', 'vector');
fprintf('Saved → plots/mortality_function.pdf\n');
