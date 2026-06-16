function analytic_bifurcations_SIS
%
%     dS/dt = Lambda - beta*S*I - mu*S + gamma*I
%     dI/dt = beta*S*I - (gamma+mu)*I - f(I)*I
%
%     f(I) = a*I^2/(1+d*I^2) + b ,   a < 0 (decreasing),   b > 0
%     d    = -m*a/((m-1)*b)         so that f(0)=b and f(Inf)=b/m
%
%
%   transcritical : R0 = 1                              (vertical line)
%   saddle-node   : det J = 0  ->  beta_SN(Ihat)
%   Hopf          : tr  J = 0  ->  beta_H (Ihat)   (with det J > 0)
%   BT (codim-2)  : beta_H(Ihat) = beta_SN(Ihat)
%   GH (codim-2)  : first Lyapunov coefficient l1(Ihat) = 0 on Hopf branch

% The GLOBAL curves in Fig. 3B (fold of limit cycles, heteroclinic/
% homoclinic) are computed numerically (via AUTO) and then plotted
% manually, so this data is only valid for the fixed parameter set below.

% ====================================================================
b      = 0.02;        % baseline per-infection DIM, f(0)=b   [week^-1]
m      = 100;         % fold reduction at high infection, f(Inf)=b/m
a      = -100;        % curvature (a<0 => decreasing DIM)
mu = 1/2600;      % recruitment                          [week^-1]
Lambda = mu;      % natural mortality                    [week^-1]

% ---- derived DIM constants and handles --------------------------------
d   = -m*a/((m-1)*b);                       % saturation parameter
f   = @(I) a.*I.^2./(1+d.*I.^2) + b;        % DIM function
fp  = @(I) 2*a.*I./(1+d.*I.^2).^2;          % f'(I)
fpp = @(I) 2*a.*(1-3*d.*I.^2)./(1+d.*I.^2).^3;        % f''(I)
fppp= @(I) 24*a*d.*I.*(d.*I.^2-1)./(1+d.*I.^2).^4;    % f'''(I)

gamma_of = @(I,beta) (beta.*(Lambda - mu.*I - I.*f(I)) - mu.^2 - mu.*f(I))./mu;
R0_of    = @(beta,gam) beta.*Lambda./(mu.*(gam + mu + b));

% ---- closed-form bifurcation branches (functions of Ihat) -------------
betaH  = @(I) -mu./I - fp(I);                          % Hopf:  tr J = 0
betaSN = @(I) -mu.*fp(I)./(mu + f(I) + I.*fp(I));      % SN:    det J = 0
detH   = @(I,beta) I.*(beta.*(mu + f(I) + I.*fp(I)) + mu.*fp(I));

% ====================================================================
Ig = linspace(1e-6, 0.20, 1200000).';

% --- Hopf branch ---
bH = betaH(Ig);  gH = gamma_of(Ig,bH);  dH = detH(Ig,bH);
okH = (bH>0) & (gH>0) & (dH>0);          % physical & genuine Hopf (det>0)
R0H = R0_of(bH,gH);

% --- Saddle-node branch ---
denSN = (mu + f(Ig) + Ig.*fp(Ig));
bSN = betaSN(Ig);  gSN = gamma_of(Ig,bSN);
okSN = (bSN>0) & (gSN>0) & (denSN>0);
R0SN = R0_of(bSN,gSN);

% --- Saddle-node TYPE via the trace of J at the fold ------------------------
% At an endemic equilibrium beta*S = gamma+mu+f(I), so J(2,2) = -I*f'(I) and
%   tr J = -(beta*I + mu + I*f'(I)).
% At a fold det J = 0, hence the eigenvalues are exactly {0, tr J}. The two
% colliding equilibria are a node (det J>0 nearby; both eigenvalues = sign of
% tr J) and a saddle. Therefore:
%   tr J < 0  ->  stable node + saddle      (ordinary fold)
%   tr J > 0  ->  unstable node + saddle    (degenerate: two unstable points)
% The switch occurs where tr J = 0 (the Hopf condition), i.e. exactly at BT.
trSN  = -(bSN.*Ig + mu + Ig.*fp(Ig));      % trace of J along the SN branch
stdSN = okSN & (trSN <  0);                % ordinary (one stable, one unstable)
degSN = okSN & (trSN >= 0);                % degenerate (both unstable)

% ====================================================================
%  Codimension-2 points
% ====================================================================
Iwin = Ig(okH);
Ilo  = min(Iwin);  Ihi = max(Iwin);

% --- BT: the Hopf branch terminates on the (low-Ihat) saddle-node fold.
% Equivalently det J -> 0 at the Hopf point, i.e. detH(I,betaH(I)) = 0.
% This holds exactly at the last physical Hopf index; bracket it there
% (robust against the beta -> Inf asymptotes of betaSN on the wide grid).
iH   = find(okH);  jbt = iH(end);
Ibt  = fzero(@(I) detH(I, betaH(I)), [Ig(jbt), Ig(jbt+1)]);
bBT  = betaH(Ibt);  gBT = gamma_of(Ibt,bBT);  R0BT = R0_of(bBT,gBT);

% --- GH: first Lyapunov coefficient l1 = 0 on the Hopf branch ---
l1H  = @(I) l1_kuznetsov(I, betaH(I), mu, f, fp, fpp, fppp);
% locate a sign change of l1 inside the physical Hopf window
Iscan = linspace(1.02*Ilo, 0.98*Ihi, 400);
Lscan = arrayfun(l1H, Iscan);
ksc   = find(diff(sign(Lscan))~=0, 1, 'first');
Igh   = fzero(l1H, [Iscan(ksc), Iscan(ksc+1)]);
bGH   = betaH(Igh);  gGH = gamma_of(Igh,bGH);  R0GH = R0_of(bGH,gGH);

% split Hopf branch into sub/supercritical at GH
subH = okH & (Ig <= Igh);                % l1 > 0  (subcritical)
supH = okH & (Ig >  Igh) & (Ig <= Ibt);  % l1 < 0  (supercritical)

fprintf('\n  Codim-2 points  (R0 , 1/gamma):\n');
fprintf('    BT : R0 = %.5f , gamma = %.5f , 1/gamma = %.4f\n', R0BT, gBT, 1/gBT);
fprintf('    GH : R0 = %.5f , gamma = %.5f , 1/gamma = %.4f\n', R0GH, gGH, 1/gGH);

% ====================================================================
%  GLOBAL curves: manually entered from numerical continuation (AUTO)
% ====================================================================
% These have no closed form; each row is a continued point [gamma, R0].
% Plotted as (R0, 1/gamma) to share the local-bifurcation axes.

% --- Fold of limit cycles (FLC) ---
%   columns: [gamma, R0]
FLC = [ 0.54, 1.309; ...
        0.50 , 1.208   ; ...
        0.45 , 1.116   ; ...
        0.40 , 1.08014 ; ...
        0.30 , 1.053   ; ...
        0.20 , 1.045   ; ...
        0.1428, 1.0375;...
        0.1111 , 1.037   ];
[~,ord] = sort(FLC(:,1));  FLC = FLC(ord,:);    % order by gamma for a clean line

% --- Heteroclinic points ---
%   columns: [gamma, R0]   (gamma=0.45 contributes two branches)
%   gamma=0.20 omitted: heteroclinic/PD status uncertain (oscillations
%   persist to R0~1.2095 where the branch hits an LP).
HET = [ 0.45 , 1.563 ; ...
        0.45 , 1.462 ; ...
        0.40 , 1.36  ; ...
        0.30 , 1.26  ; ...
        0.20, 1.2095;...
        0.1111 , 1.1313 ];

% ====================================================================
%  PLOT  (1/gamma vs R0)
% ====================================================================

%% ── Style ────────────────────────────────────────────────────────────────

set_ggplot_style();
gg = evalin('base', 'gg_colours');

LW = 2;
FN = 'Helvetica';
FS = 12;

%% ── Figure ───────────────────────────────────────────────────────────────

fig = figure('Units', 'inches', 'Position', [1 1 6 5], 'Color', 'w');
ax  = axes(fig);
hold(ax, 'on');
ax.TickDir = 'out';

% transcritical: R0 = 1
yl = [0 9];
plot(ax, [1 1], yl, '-', 'Color', gg(1,:), 'LineWidth', LW, ...
     'DisplayName', 'transcritical');

% saddle-node: the det J = 0 locus has TWO physical folds (n=2 -> cubic
% equilibrium equation). Each fold is coloured by its TYPE (sign of tr J):
%   ordinary (stable node + saddle)  -> olive, solid
%   degenerate (two unstable points) -> olive, dashed
% Disjoint physical runs in Ihat are drawn separately but share one legend
% entry per type. The two types necessarily meet at the BT point.
nStd = plot_sn_segments(ax, stdSN, R0SN, gSN, gg(2,:), '-',  LW, ...
                        'saddle-node (node + saddle)');
nDeg = plot_sn_segments(ax, degSN, R0SN, gSN, gg(2,:), '--', LW, ...
                        'saddle-node (both unstable)');
fprintf('  saddle-node segments: %d ordinary, %d degenerate\n', nStd, nDeg);

% Hopf, sub / supercritical
plot(ax, R0H(subH), 1./gH(subH), '-', 'Color', gg(3,:), 'LineWidth', LW, ...
     'DisplayName', 'Hopf (subcritical)');
plot(ax, R0H(supH), 1./gH(supH), '-', 'Color', gg(4,:), 'LineWidth', LW, ...
     'DisplayName', 'Hopf (supercritical)');

% fold of limit cycles (FLC): global curve, manual continuation points
plot(ax, FLC(:,2), 1./FLC(:,1), '-o', 'Color', gg(5,:), ...
     'MarkerFaceColor', gg(5,:), 'MarkerSize', 6, 'LineWidth', LW, ...
     'DisplayName', 'fold of limit cycles');

% heteroclinic points: global, no connecting line (branches not contiguous)
plot(ax, HET(:,2), 1./HET(:,1), '-ks', 'MarkerFaceColor', 'k', ...
     'LineWidth', LW, 'DisplayName', 'heteroclinic');

% codim-2 markers
plot(ax, R0BT, 1/gBT, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 13, ...
     'HandleVisibility', 'off');
plot(ax, R0GH, 1/gGH, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 13, ...
     'HandleVisibility', 'off');
text(ax, R0BT, 1/gBT, '  BT', 'FontName', FN, 'FontSize', FS, ...
     'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
text(ax, R0GH, 1/gGH, '  GH', 'FontName', FN, 'FontSize', FS, ...
     'FontWeight', 'bold', 'VerticalAlignment', 'top');

xlim(ax, [0.8 1.7]); ylim(ax, yl);
xticks(ax, 0.8 : 0.2 : 1.6);   % 0.8, 1.0, 1.2, 1.4, 1.6
yticks(ax, 1 : 2 : 9);         % 1, 3, 5, 7, 9
xlabel(ax, 'R_0', 'Interpreter', 'tex');
ylabel(ax, 'avg. infectious period (1/\gamma, weeks)', 'Interpreter', 'tex');

legend(ax, 'Interpreter', 'tex', 'Location', 'northeast', ...
       'FontSize', FS, 'Box', 'off');

hold(ax, 'off');
end
% ======================================================================

function nrun = plot_sn_segments(ax, mask, R0SN, gSN, col, style, LW, name)
% Draw one saddle-node TYPE (possibly several disjoint folds in Ihat) in a
% single colour/line-style with one shared legend entry. Each contiguous run
% in the index set is sorted by R0 so the polyline is clean. Returns the
% number of disjoint runs drawn.
    idx = find(mask);
    nrun = 0;
    if isempty(idx), return; end
    brk = [0; find(diff(idx) > 1); numel(idx)];
    nrun = numel(brk) - 1;
    for ks = 1:nrun
        sid = idx(brk(ks)+1 : brk(ks+1));
        [~,ord] = sort(R0SN(sid));      % order for a clean line
        sid = sid(ord);
        if ks == 1
            plot(ax, R0SN(sid), 1./gSN(sid), style, 'Color', col, ...
                 'LineWidth', LW, 'DisplayName', name);
        else
            plot(ax, R0SN(sid), 1./gSN(sid), style, 'Color', col, ...
                 'LineWidth', LW, 'HandleVisibility', 'off');
        end
    end
end
% ======================================================================

function l1 = l1_kuznetsov(I0, beta, mu, f, fp, fpp, fppp)
% Nonlinearity of the (S,I) field about the endemic equilibrium:
%   quadratic:  P = -beta*s*i ;  Q =  beta*s*i - (hpp/2) i^2
%   cubic:      Q = -(hppp/6) i^3      (with s=S-S*, i=I-I*)
% where h(I)=I*f(I), hpp = 2 f' + I f'' , hppp = 3 f'' + I f'''.
    hpp  = 2*fp(I0) + I0*fpp(I0);
    hppp = 3*fpp(I0) + I0*fppp(I0);

    A = [ -(beta*I0+mu), -(mu+f(I0)); ...
           beta*I0,      -fp(I0)*I0 ];
    w = sqrt(A(1,1)*A(2,2) - A(1,2)*A(2,1));   % det>0 at Hopf
    if ~isreal(w) || w<=0, l1 = NaN; return; end

    % right eigenvector q: A q = i w q
    [V,D] = eig(A);  [~,kq] = min(abs(diag(D) - 1i*w));  q = V(:,kq);
    % left  eigenvector p: A' p = -i w p
    [W,Dl] = eig(A.'); [~,kp] = min(abs(diag(Dl) - (-1i*w))); p = W(:,kp);
    p = p / conj(p'*q);                         % <p,q> = p'*q = 1

    B = @(x,y)[ -beta*(x(1)*y(2)+x(2)*y(1)); ...
                 beta*(x(1)*y(2)+x(2)*y(1)) - hpp*x(2)*y(2) ];
    C = @(x,y,z)[ 0; -hppp*x(2)*y(2)*z(2) ];

    qb = conj(q);  E = eye(2);
    g21 = p'*C(q,q,qb) ...
        - 2*(p'*B(q, A\B(q,qb))) ...
        +    p'*B(qb, (2i*w*E - A)\B(q,q));
    l1 = real(g21)/(2*w);
end