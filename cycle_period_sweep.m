function [periods_wk, info] = cycle_period_sweep(model, param_name, param_vals, base, opts)
%CYCLE_PERIOD_SWEEP  Period of the stable limit cycle vs a system parameter.
%
%   For the SIS / SIRS models with non-constant disease-induced mortality
%   (DIM), this function sweeps one parameter (R0, gamma or delta), integrates
%   the ODEs onto the attractor, and returns the period of the stable limit
%   cycle (NaN where the attractor is an equilibrium, i.e. no cycle).
%
%   Models (decreasing DIM, f(I) = a I^n/(1 + d I^n) + b, a < 0):
%     SIS    S' = Lam - beta S I - mu S + gamma I
%            I' = beta S I - (gamma+mu) I - f(I) I
%     SIRS   S' = Lam - beta S I - mu S + delta R
%            I' = beta S I - (gamma+mu) I - f(I) I
%            R' = gamma I - (mu+delta) R           (SIR when delta = 0)
%
%   USAGE
%     [T, info] = cycle_period_sweep('sirs','R0', 1.32:0.01:1.54, base, opts)
%
%   INPUTS
%     model       'sis' or 'sirs'
%     param_name  'R0', 'gamma' or 'delta' (the swept parameter)
%     param_vals  vector of values for the swept parameter
%     base        struct of fixed parameters with fields:
%                   Lam, mu, gamma, delta, R0, a, b, m, n
%                 (delta ignored for 'sis'; the swept field is overwritten)
%     opts        (optional) struct of solver/detector options; see defaults
%                 below. Notable fields:
%                   y0          initial [S I] (sis) or [S I R] (sirs) used to
%                               SEED the first sweep point. Choose it inside
%                               the cycle region so the warm start lands on the
%                               cycle (see run_cycle_periods.m for examples).
%                   reverse     true to sweep param_vals back-to-front
%                   Tbase, Tmax integration horizon (weeks); auto-extended to
%                               capture long-period (near-homoclinic) cycles
%                   Nout        number of uniform output samples
%                   tail        fraction of the record used for detection
%                   relamp_min  min relative amplitude to count as a cycle
%
%   OUTPUTS
%     periods_wk  vector of cycle periods in WEEKS (NaN = no stable cycle)
%     info        struct: .param_vals, .periods_yr, .rel_amp, .ncross,
%                          .Iend (endpoint I), .model, .param_name
%
%   NOTES
%     * beta is set from R0 via  beta = R0 (gamma + mu + f(0)) / (Lam/mu),
%       with f(0) = b  (Eq. 3 of the manuscript; S0 = Lam/mu).
%     * d is fixed by (a,m,b):  d = -a m / (b (m-1)),  so f(0)=b, f(inf)=b/m.
%       (The Table 1 value d ~ 5.05e5 is inconsistent with a=-100 under this
%       relation, which gives d ~ 5050; this routine recomputes d from a,m,b.)
%     * WARM START: each sweep point is integrated from the previous point's
%       endpoint. This continues ONE attractor branch. In bistable regions
%       (stable EE coexisting with a cycle) the branch you trace depends on the
%       seed y0 and the sweep direction; sweep both ways to map hysteresis.
%     * Uses only base MATLAB (ode15s). No toolboxes required. Octave-friendly.
%
%   Companion: run_cycle_periods.m (driver that makes the figures).

% ----------------------- defaults -----------------------
if nargin < 5, opts = struct(); end
def = struct('Tbase',6e4,'Tmax',4e5,'Nout',1.5e5,'tail',0.6, ...
             'relamp_min',1e-3,'reltol',1e-9,'abstol',1e-12, ...
             'reverse',false,'min_cross',3,'max_extend',3,'verbose',true);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end
% the swept parameter is supplied per-iteration (p.(param_name)=pv(i)) so it
% need not be present in base; require only the remaining fixed parameters.
need = setdiff({'Lam','mu','gamma','delta','R0','a','b','m','n'}, ...
               {param_name}, 'stable');
for k = 1:numel(need)
    if ~isfield(base,need{k}), error('base.%s is required',need{k}); end
end
model = lower(model);
if ~ismember(model,{'sis','sirs'}), error('model must be ''sis'' or ''sirs'''); end

pv = param_vals(:).';
if opts.reverse, pv = fliplr(pv); end
np = numel(pv);

% default seed
if ~isfield(opts,'y0') || isempty(opts.y0)
    if strcmp(model,'sis'), opts.y0 = [0.985, 0.012];
    else,                   opts.y0 = [0.70, 0.02, 0.10]; end
end

periods_wk = nan(1,np);
rel_amp    = nan(1,np);
ncross     = zeros(1,np);
Iend       = nan(1,np);

y = opts.y0(:).';
odeopt = odeset('RelTol',opts.reltol,'AbsTol',opts.abstol);

for i = 1:np
    p = base;                       % copy fixed params
    p.(param_name) = pv(i);         % overwrite swept param

    [d,~] = dim_d(p.a,p.b,p.m); 
    beta  = p.R0*(p.gamma + p.mu + p.b)/(p.Lam/p.mu);
    f     = @(I) dim_f(I,p.a,p.b,p.m,p.n);

    if strcmp(model,'sis')
        rhs = @(t,Y) [ p.Lam - beta*Y(1)*Y(2) - p.mu*Y(1) + p.gamma*Y(2); ...
                       beta*Y(1)*Y(2) - (p.gamma+p.mu)*Y(2) - f(Y(2))*Y(2) ];
        Iidx = 2;
    else
        rhs = @(t,Y) [ p.Lam - beta*Y(1)*Y(2) - p.mu*Y(1) + p.delta*Y(3); ...
                       beta*Y(1)*Y(2) - (p.gamma+p.mu)*Y(2) - f(Y(2))*Y(2); ...
                       p.gamma*Y(2) - (p.mu+p.delta)*Y(3) ];
        Iidx = 2;
    end

    % adaptive integration: extend horizon until >~12 periods are captured
    T = opts.Tbase; per = NaN; amp = 0; nc = 0; Yend = y;
    for attempt = 1:opts.max_extend
        tg = linspace(0,T,opts.Nout);
        sol = ode15s(rhs,[0 T],y,odeopt);
        Y   = deval(sol,tg);
        Yend = Y(:,end).';
        [per,nc,amp] = detect_period(tg, Y(Iidx,:), opts.tail, ...
                                     opts.relamp_min, opts.min_cross);
        if ~isnan(per) && T < 12*per
            T = min(15*per, opts.Tmax); continue
        end
        break
    end

    y = Yend;                       % warm start for next parameter value
    periods_wk(i) = per;
    rel_amp(i)    = amp;
    ncross(i)     = nc;
    Iend(i)       = Yend(Iidx);

    if opts.verbose
        fprintf('%-4s %s=%.4f : T=%8.2f yr  relamp=%.2e  ncross=%d\n', ...
            model, param_name, pv(i), per/52, amp, nc);
    end
end

% undo reversal so outputs align with the input order of param_vals
if opts.reverse
    pv         = fliplr(pv);
    periods_wk = fliplr(periods_wk);
    rel_amp    = fliplr(rel_amp);
    ncross     = fliplr(ncross);
    Iend       = fliplr(Iend);
end

info = struct('param_vals',pv,'periods_yr',periods_wk/52,'rel_amp',rel_amp, ...
              'ncross',ncross,'Iend',Iend,'model',model,'param_name',param_name);
end

% ===================== local functions =====================

function v = dim_f(I,a,b,m,n)
% Hill-type DIM function f(I) = a I^n/(1 + d I^n) + b  (decreasing for a<0).
d = dim_d(a,b,m);
v = a.*I.^n ./ (1 + d.*I.^n) + b;
end

function [d,fwd] = dim_d(a,b,m)
% Saturation constant fixed by f(0)=b and f(inf)=b/m (a<0): d = -a m/(b(m-1)).
d   = -a*m/(b*(m-1));
fwd = b/m;          % limiting value f(inf)
end

function [per, nc, relamp] = detect_period(t, y, tail, relamp_min, min_cross)
%DETECT_PERIOD  Period from a (possibly anharmonic) oscillatory time series.
%   Uses upward crossings of the signal mean (linearly interpolated). Robust
%   for long, spiky, near-homoclinic cycles. Returns NaN if the record looks
%   like a fixed point (relative amplitude below relamp_min) or too few cycles.
n  = numel(t);
i0 = max(1, floor(n*(1-tail)));
tt = t(i0:end);  yy = y(i0:end);
rng = max(yy) - min(yy);
mn  = mean(yy);
relamp = rng / max(abs(mn), eps);
if relamp < relamp_min
    per = NaN; nc = 0; return         % settled to an equilibrium
end
c = yy - mn;
cr = [];                              % upward zero-crossing times
for k = 1:numel(c)-1
    if c(k) < 0 && c(k+1) >= 0
        cr(end+1) = tt(k) + (tt(k+1)-tt(k)) * (-c(k))/(c(k+1)-c(k)); %#ok<AGROW>
    end
end
nc = numel(cr);
if nc < min_cross
    per = NaN; return
end
dintervals = diff(cr(2:end));         % drop first interval (often transient)
per = median(dintervals);
end
