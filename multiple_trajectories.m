%% SIR model with disease-induced mortality: max/min infected vs R0
%
%   Sweeps R0 across a range, integrates each trajectory to large time,
%   and records the maximum and minimum of I(t) in the second half of
%   the time series (after transients have decayed).
%
%   Outputs: plot of I_max and I_min vs R0

%% ── Style ────────────────────────────────────────────────────────────────

set_ggplot_style();
gg = evalin('base', 'gg_colours');

%% ── Fixed parameters ─────────────────────────────────────────────────────

Lambda  = 0.0003846154;
mu      = Lambda;
delta   = 0.045;
gamma_r = 1;
a       = -100000;
b       = 0.02;
m       = 100;

d = -a * m / ((m - 1) * b);

f    = @(I, beta_) b + (a .* I.^2 ./ (1 + d .* I.^2));
dfdI = @(I)        2*a.*I ./ (1 + d.*I.^2).^2;

%% ── R0 sweep ─────────────────────────────────────────────────────────────

R0_vec = linspace(1.008, 1.0125, 500);
nR0    = numel(R0_vec);

I_max = nan(1, nR0);
I_min = nan(1, nR0);

u0    = [0.978; 0.00025; 0.00877];
tspan = [0, 150000];

odeOpts = odeset('RelTol',      1e-12,   ...
                 'AbsTol',      1e-12,   ...
                 'NonNegative', [1 2 3], ...
                 'MaxOrder',    5);

fprintf('Sweeping R0 from %.3f to %.3f (%d values)...\n', R0_vec(1), R0_vec(end), nR0);

for k = 1:nR0

    beta = R0_vec(k) * (gamma_r + mu + b);

    odefun = @(t, u) [ ...
        Lambda - beta*u(1)*u(2) - mu*u(1) + delta*u(3); ...
        beta*u(1)*u(2) - (gamma_r + mu)*u(2) - f(u(2))*u(2); ...
        gamma_r*u(2)   - (mu + delta)*u(3) ];

    jac = @(t, u) [ ...
        -(beta*u(2) + mu),  -beta*u(1),      delta; ...
         beta*u(2),  beta*u(1) - (gamma_r+mu) - f(u(2)) - dfdI(u(2))*u(2),  0; ...
         0,          gamma_r,  -(mu+delta) ];

    opts_k = odeset(odeOpts, 'Jacobian', jac);

    [t, u] = ode15s(odefun, tspan, u0, opts_k);

    I   = u(:, 2);
    mid = find(t >= tspan(end)/2, 1);

    I_max(k) = max(I(mid:end));
    I_min(k) = min(I(mid:end));

    if mod(k, 10) == 0
        fprintf('  R0 = %.3f  →  I_max = %.6f,  I_min = %.6f\n', ...
            R0_vec(k), I_max(k), I_min(k));
    end
end

fprintf('Done.\n\n');

%% ── Figure ───────────────────────────────────────────────────────────────

LW = 1.2;
FS = 10;
FN = 'Helvetica';

fig = figure('Units', 'inches', 'Position', [1 1 7 5], 'Color', 'w');
ax  = axes(fig);
hold(ax, 'on');
ax.TickDir = 'out';

plot(ax, R0_vec, I_max, '-', 'Color', gg(1,:), 'LineWidth', LW);
plot(ax, R0_vec, I_min, '-', 'Color', gg(2,:), 'LineWidth', LW);

xlabel(ax, 'Basic reproduction number, \itR\rm_0', 'Interpreter', 'tex', ...
    'FontName', FN, 'FontSize', FS);
ylabel(ax, 'Infected, \itI',                       'Interpreter', 'tex', ...
    'FontName', FN, 'FontSize', FS);

legend(ax, {'\itI\rm_{max}', '\itI\rm_{min}'}, ...
    'Location', 'northwest', 'Box', 'off', ...
    'FontName', FN, 'FontSize', FS);

xl = ax.XLim;  yl = ax.YLim;
text(ax, xl(1), yl(2), 'A', ...
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
% exportgraphics(fig, 'plots/multiple_trajectories.pdf', 'ContentType', 'vector');
% fprintf('Saved → plots/multiple_trajectories.pdf\n');
