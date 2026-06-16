%% SIS model with disease-induced mortality: single-trajectory analysis
%
%   ODE system
%   ----------
%   d    = -a*m / ((m-1)*b)
%   f(I) =  b + a*I^2 / (1 + d*I^2)          % I-dependent mortality rate
%
%   S'  = Lambda - beta*S*I - mu*S + gamma_r*I
%   I'  = beta*S*I - (gamma_r+mu)*I - f(I)*I
%
%   where  beta = R0*(mu/Lambda)*(gamma_r+mu+b)
%
%   At the disease-free equilibrium (DFE): (S*, I*) = (1, 0)
%   since Lambda = mu.
%
%   Outputs: plots/single_trajectory_SIS.pdf

%% ── Style ────────────────────────────────────────────────────────────────

%addpath('/Users/denispatterson/Documents/MatCont7p4/');
set_ggplot_style();
gg = evalin('base', 'gg_colours');

%% ── Parameters ───────────────────────────────────────────────────────────

Lambda  = 0.0003846154;
R0_val  = 1.1313;
mu      = Lambda;
gamma_r = 0.11111;  % 'gamma_r' avoids conflicts with the built-in gamma()
a       = -100;
b       = 0.02;
m       = 100;

%% ── Derived quantities ───────────────────────────────────────────────────

% Mortality shape coefficient (always positive when a < 0)
d    = -a * m / ((m - 1) * b);

% Effective transmission coefficient
beta = R0_val * (mu / Lambda) * (gamma_r + mu + b);

% Disease-induced mortality rate as a function of I
%   f(0) = b  (baseline),  f(inf) -> b/m  (saturating, lower at high I)
f    = @(I) b + (a .* I.^2 ./ (1 + d .* I.^2));

% Derivative of f w.r.t. I  (needed for analytical Jacobian)
%   f'(I) = 2aI / (1 + dI²)²
dfdI = @(I) 2*a.*I ./ (1 + d.*I.^2).^2;

% Analytical Jacobian of the ODE right-hand side
%
%   J = [ -(β I + µ),   -β S + γ                              ]
%       [   β I,         β S-(γ+µ)-f(I)-f'(I)·I               ]
%
jac  = @(t, u) [ ...
    -(beta*u(2) + mu),   -beta*u(1) + gamma_r; ...
     beta*u(2),           beta*u(1) - (gamma_r+mu) - f(u(2)) - dfdI(u(2))*u(2) ];

fprintf('\n── Derived parameters ──\n');
fprintf('  d          = %12.4f\n', d);
fprintf('  beta       = %12.8f\n', beta);
fprintf('  f(I=0)     = %12.6f  (= b)\n', f(0));
fprintf('  f(I->inf)  = %12.6f  (= b/m)\n', b/m);
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% ── ODE right-hand side ──────────────────────────────────────────────────

odefun = @(t, u) [ ...
    Lambda - beta*u(1)*u(2) - mu*u(1) + gamma_r*u(2); ...
    beta*u(1)*u(2) - (gamma_r + mu)*u(2) - f(u(2))*u(2) ];

%% ── Fixed points and their stability ─────────────────────────────────────
%
%   DFE:      I = 0,  S = Lambda/mu.
%   Endemic:  beta*S = gamma_r+mu+f(I)  =>  S = (gamma_r+mu+f(I))/beta,
%             and substituting into S'=0 gives the scalar root condition
%               g(I) = beta*(Lambda - (mu+f(I))*I) - mu*(gamma_r+mu+f(I)) = 0.
%   The n=2 mortality makes g(I) cubic-like, so several endemic roots may
%   coexist. Each fixed point is classified via the Jacobian eigenvalues:
%   stable (both Re < 0) -> red dot, unstable -> black dot.

eqS = Lambda/mu;                 % DFE susceptible level
eq_pts = [eqS, 0];               % rows: [S, I]

gI    = @(I) beta*(Lambda - (mu + f(I)).*I) - mu*(gamma_r + mu + f(I));
Iscan = linspace(1e-9, 0.20, 200000);
gv    = gI(Iscan);
sc    = find(diff(sign(gv)) ~= 0);   % bracket each sign change
for k = 1:numel(sc)
    Iroot = fzero(gI, [Iscan(sc(k)), Iscan(sc(k)+1)]);
    Sroot = (gamma_r + mu + f(Iroot)) / beta;
    eq_pts(end+1, :) = [Sroot, Iroot];   %#ok<SAGROW>
end

% classify each equilibrium: stable if both eigenvalues have Re < 0
eq_stable = false(size(eq_pts, 1), 1);
for k = 1:size(eq_pts, 1)
    ev = eig(jac(0, eq_pts(k, :).'));
    eq_stable(k) = all(real(ev) < 0);
end

fprintf('── Fixed points (S, I)  [stability] ─────────────────────────\n');
for k = 1:size(eq_pts, 1)
    lbl = 'unstable'; if eq_stable(k), lbl = 'stable'; end
    fprintf('  S = %.6f,  I = %.6f   [%s]\n', eq_pts(k,1), eq_pts(k,2), lbl);
end
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% ── Initial conditions and integration ───────────────────────────────────

u0   = [0.85; 0.01];   % S(0), I(0)
tspan = [0, 300000];

odeOpts = odeset('RelTol',      1e-12,   ...
                 'AbsTol',      1e-12,   ...
                 'NonNegative', [1 2],   ...
                 'Jacobian',    jac,     ...
                 'MaxOrder',    5);

[t, u] = ode15s(odefun, tspan, u0, odeOpts);

S = u(:, 1);
I = u(:, 2);

mid = find(t >= tspan(end)/2, 1);   % first index in the second half of tspan

fprintf('── Endpoint values at t = %g ─────────────────────────────────\n', tspan(end));
fprintf('  S = %.6f,  I = %.6f\n', S(end), I(end));
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% ── Colour assignments ───────────────────────────────────────────────────

col_SI  = gg(1, :);   % pink/salmon  — phase portrait (panel A)
col_S   = gg(2, :);   % olive-green  — (unused in time-series panels)

col_S_B = gg(4, :);   % purple       — S(t) in panel B
col_I_B = gg(5, :);   % warm orange  — I(t) in panel B

%% ── Figure ───────────────────────────────────────────────────────────────
%
%   Layout (1 × 2 grid):
%     A — 2D phase portrait (S,I)  |  B — S(t) and I(t) on dual y-axes

fig = figure('Units', 'inches', 'Position', [1 1 10 4], 'Color', 'w');

LW = 2;   % line width
FS = 12;    % font size
FN = 'Helvetica';

%% Panel A — 2D phase portrait in (S, I) plane ────────────────────────────

ax1 = subplot(1, 2, 1);
hold(ax1, 'on');
ax1.TickDir = 'out';

plot(ax1, S(mid:end), I(mid:end), '-', 'Color', col_SI, 'LineWidth', LW);

% fixed points: red = stable, black = unstable
col_stable   = [0.85 0.15 0.15];   % red
col_unstable = [0 0 0];            % black
for k = 1:size(eq_pts, 1)
    if eq_stable(k)
        cmk = col_stable;
    else
        cmk = col_unstable;
    end
    plot(ax1, eq_pts(k,1), eq_pts(k,2), 'o', ...
        'MarkerFaceColor', cmk, 'MarkerEdgeColor', cmk, 'MarkerSize', 9);
end

xlabel(ax1, 'Susceptible,  (S)', 'Interpreter', 'tex');
ylabel(ax1, 'Infected, (I)',   'Interpreter', 'tex');

yt = ax1.YTick;
if numel(yt) > 1,  ax1.YTick = yt(2:end);  end

xl1 = ax1.XLim;  yl1 = ax1.YLim;
text(ax1, xl1(1), yl1(2), 'A', ...
    'FontName', FN, 'FontSize', 12, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Panel B — S(t) and I(t) on dual y-axes ────────────────────────────────

ax2 = subplot(1, 2, 2);
ax2.TickDir = 'out';
ax2.Box     = 'off';

t_pad = 0.05 * tspan(end);

yyaxis(ax2, 'left');
plot(ax2, t, S, '-', 'Color', col_S_B, 'LineWidth', LW);
ylabel(ax2, 'Susceptibles (S)', 'Interpreter', 'tex');
ax2.YAxis(1).Color = col_S_B;
yticks(ax2, linspace(ax2.YLim(1), ax2.YLim(2), 5));

yyaxis(ax2, 'right');
plot(ax2, t, I, '-', 'Color', col_I_B, 'LineWidth', LW);
ylabel(ax2, 'Infected (I)', 'Interpreter', 'tex');
ax2.YAxis(2).Color = col_I_B;
yticks(ax2, linspace(ax2.YLim(1), ax2.YLim(2), 5));

xlim(ax2, [-t_pad, tspan(end) + t_pad]);
xticks(ax2, linspace(0, tspan(end), 5));
xlabel(ax2, 'Time (weeks)', 'Interpreter', 'tex');

yyaxis(ax2, 'left');
xl2 = ax2.XLim;  yl2 = ax2.YLim;
text(ax2, xl2(1), yl2(2), 'B', ...
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
% exportgraphics(fig, 'plots/single_trajectory_SIS.pdf', 'ContentType', 'vector');
% fprintf('Saved → plots/single_trajectory_SIS.pdf\n');
