function analytic_bifurcations_SIR

% ====================================================================
b      = 0.02;     % baseline per-infection DIM, f(0)=b      [week^-1]
m      = 100;      % fold reduction at high infection, f(Inf)=b/m
a      = -100;     % curvature (a<0 => decreasing DIM)
gamma  = 0.25;        % recovery rate                            [week^-1]
mu     = 1/2600;   % natural mortality                        [week^-1]
Lambda = mu;   % recruitment                              [week^-1]

% plotting window
R0lim = [0.9 5];
Ylim  = [1 25];          % 1/delta

% ---- DIM handles ----
d    = -a*m/(b*(m-1));
f    = @(I) a.*I.^2./(1+d.*I.^2) + b;
fp   = @(I) 2*a.*I./(1+d.*I.^2).^2;
fpp  = @(I) 2*a.*(1-3*d.*I.^2)./(1+d.*I.^2).^3;
fppp = @(I) 24*a*d.*I.*(d.*I.^2-1)./(1+d.*I.^2).^4;

% equilibrium-determined beta and R0 (functions of Ihat, delta)
beta_eq = @(I,dl) mu.*(f(I)+gamma+mu) ./ ...
                  (Lambda + (dl./(mu+dl)).*gamma.*I - I.*(f(I)+gamma+mu));
R0_of   = @(be) be.*Lambda./(mu.*(gamma+mu+b));

% ---- analytic I-derivatives of beta_eq (for the CUSP conditions) ----------
% Along the equilibrium graph beta = N(I)/D(I), a fold (saddle-node) is
% d(beta)/dI = 0 and a cusp (codim-2, where two folds merge) additionally
% satisfies d^2(beta)/dI^2 = 0. Writing N = mu*(f+gamma+mu) and
% D = Lambda + k*I - I*(f+gamma+mu) with k = delta/(mu+delta)*gamma:
cc   = gamma + mu;
Nf   = @(I)    mu.*(f(I)+cc);
Nfp  = @(I)    mu.*fp(I);
Nfpp = @(I)    mu.*fpp(I);
kk   = @(dl)   (dl./(mu+dl)).*gamma;
Df   = @(I,dl) Lambda + kk(dl).*I - I.*(f(I)+cc);
Dfp  = @(I,dl) kk(dl) - f(I) - I.*fp(I) - cc;
Dfpp = @(I,dl) -2*fp(I) - I.*fpp(I);
beta_I  = @(I,dl) (Nfp(I).*Df(I,dl) - Nf(I).*Dfp(I,dl)) ./ Df(I,dl).^2;
beta_II = @(I,dl) (Nfpp(I).*Df(I,dl).^2 - Nf(I).*Dfpp(I,dl).*Df(I,dl) ...
                   - 2*Dfp(I,dl).*(Nfp(I).*Df(I,dl) - Nf(I).*Dfp(I,dl))) ./ Df(I,dl).^3;

% ====================================================================
%  GRID in (Ihat, delta) and Routh-Hurwitz fields
% ====================================================================
nI = 1600; nD = 1600;
Iv = linspace(1e-8, 0.3, nI);
dv = linspace(0.01 , 1 , nD);
[II,DD] = meshgrid(Iv,dv);

B  = beta_eq(II,DD);
[A1,A2,A3,S] = rh_coeffs(II,DD,B,f,fp,mu,gamma);
RR = R0_of(B);
phys = (B>0) & (S>0) & isfinite(RR);

% masked fields for contouring
A3_ns = A3;  A3_ns(~(phys & A2>0)) = NaN;      % SN: node + saddle
A3_uu = A3;  A3_uu(~(phys & A2<0)) = NaN;      % SN: both unstable
A3_all= A3;  A3_all(~phys)         = NaN;      % SN: full (for BT)

% ====================================================================
%  PLOT
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
plot(ax, [1 1], Ylim, '-', 'Color', gg(1,:), 'LineWidth', LW, ...
     'DisplayName', 'transcritical');

% saddle-node, coloured by a2 sign at the fold (matches the SIS convention):
%   a2 > 0 -> node + saddle (ordinary)   -> olive, solid
%   a2 < 0 -> both unstable (degenerate)  -> olive, dashed
% Contours in (Ihat,delta) are mapped to (R0,1/delta).
plot_mapped_contour(Iv,dv,A3_ns, gg(2,:), '-',  LW, beta_eq,R0_of, 'saddle-node (node + saddle)');
plot_mapped_contour(Iv,dv,A3_uu, gg(2,:), '--', LW, beta_eq,R0_of, 'saddle-node (both unstable)');

% ---- BT: on the saddle-node locus where a2 changes sign ----
Psn = parse_contourc(contourc(Iv,dv,A3_all,[0 0]));
R0bt = NaN; BT = [NaN NaN];
for k = 1:numel(Psn)
    Ip = Psn{k}(1,:); Dp = Psn{k}(2,:);
    be = beta_eq(Ip,Dp); [~,a2l] = rh_coeffs(Ip,Dp,be,f,fp,mu,gamma);
    s = find(diff(sign(a2l))~=0);
    for j = s
        x = newton2(@(x) bt_res(x,beta_eq,f,fp,mu,gamma), [Ip(j);Dp(j)]);
        BT = [R0_of(beta_eq(x(1),x(2))), 1/x(2)]; R0bt = BT(1);
    end
end
plot(ax, BT(1),BT(2),'kp','MarkerFaceColor','k','MarkerSize',13,'HandleVisibility','off');
text(ax, BT(1),BT(2),'  BT','FontName',FN,'FontSize',FS,'FontWeight','bold','VerticalAlignment','bottom');

% ---- CP (cusp): codim-2 point where the two saddle-node folds merge --------
% The fold locus is d(beta)/dI = 0; the cusp also has d^2(beta)/dI^2 = 0.
% Solve [beta_I; beta_II] = 0 in (Ihat,delta). Along the a3=0 contour the
% tangent obeys d(delta)/ds ~ beta_II, so delta is extremal exactly at the
% cusp; seed Newton at the delta- (and R0-) turning points of the contour.
cusp_res = @(x) [ beta_I(x(1),x(2)); beta_II(x(1),x(2)) ];
CP = [NaN NaN]; R0cp = NaN; bestres = inf;
for k = 1:numel(Psn)
    Ip = Psn{k}(1,:); Dp = Psn{k}(2,:);
    r0 = R0_of(beta_eq(Ip,Dp));
    cD = find(diff(sign(diff(Dp))) ~= 0) + 1;     % interior delta extrema (cusp)
    cR = find(diff(sign(diff(r0))) ~= 0) + 1;     % interior R0 extrema (backup)
    [~,jd] = max(Dp); [~,jm] = max(r0);
    cand = unique([cD, cR, jd, jm]);              % seed candidates
    for j = cand
        x = newton2(cusp_res, [Ip(j); Dp(j)]);
        r = cusp_res(x);
        if all(isfinite(x)) && x(1) > 0 && x(2) >= dv(1) && x(2) <= dv(end) ...
                && norm(r) < bestres
            bestres = norm(r);
            CP = [R0_of(beta_eq(x(1),x(2))), 1/x(2)]; R0cp = CP(1);
        end
    end
end
plot(ax, CP(1),CP(2),'kp','MarkerFaceColor','k','MarkerSize',13,'HandleVisibility','off');
text(ax, CP(1),CP(2),'  CP','FontName',FN,'FontSize',FS,'FontWeight','bold', ...
     'VerticalAlignment','bottom','HorizontalAlignment','right');

% ---- Hopf: trace the locus as the single-valued graph delta = delta_Hopf(I) -
% The Routh-Hurwitz Hopf condition a1*a2-a3 = 0 (with a1,a2>0) fixes, for each
% endemic level Ihat, a UNIQUE waning rate delta -- the locus is single-valued
% in Ihat (it folds only in delta). We therefore trace it by sweeping Ihat and
% solving the 1-D root delta(Ihat), which returns ONE continuous, correctly
% ordered curve. This avoids contouring a 2-D field (where the near-vertical
% locus fragments into dozens of pieces) and the subsequent stitching.
hopf = trace_hopf(Iv, [dv(1) dv(end)], beta_eq, f, fp, mu, gamma);
Ip = hopf(1,:); Dp = hopf(2,:);
r0 = R0_of(beta_eq(Ip,Dp)); yy = 1./Dp;
L  = arrayfun(@(I,D) l1_3d(I,D,beta_eq,f,fp,fpp,fppp,mu,gamma), Ip, Dp);

% GH: the single l1=0 crossing on the curve, Newton-refined, splits it into a
% subcritical and a supercritical segment that share the GH vertex (no gap).
GH = [NaN NaN]; firstSub = true; firstSup = true;
s = find(diff(sign(L))~=0);
splitj = []; bestgh = inf;
for j = s
    if abs(L(j))<5e4 && abs(L(j+1))<5e4
        x = newton2(@(x) gh_res(x,beta_eq,f,fp,fpp,fppp,mu,gamma), [Ip(j);Dp(j)]);
        res = gh_res(x,beta_eq,f,fp,fpp,fppp,mu,gamma);
        if all(isfinite(x)) && norm(res) < bestgh
            bestgh = norm(res); splitj = j;
            GH = [R0_of(beta_eq(x(1),x(2))), 1/x(2)];
        end
    end
end

if isempty(splitj)
    issub = median(sign(L(isfinite(L)))) > 0;
    [firstSub,firstSup] = draw_crit(ax, r0, yy, issub, gg, LW, firstSub, firstSup);
else
    gx = GH(1); gy = GH(2);
    rA = [r0(1:splitj), gx];     yA = [yy(1:splitj), gy];
    rB = [gx, r0(splitj+1:end)]; yB = [gy, yy(splitj+1:end)];
    Asub = sign(L(splitj)) > 0;     % side containing point j is subcritical?
    [firstSub,firstSup] = draw_crit(ax, rA, yA,  Asub, gg, LW, firstSub, firstSup);
    [firstSub,firstSup] = draw_crit(ax, rB, yB, ~Asub, gg, LW, firstSub, firstSup);
end
plot(ax, GH(1),GH(2),'kp','MarkerFaceColor','k','MarkerSize',13,'HandleVisibility','off');
text(ax, GH(1),GH(2),'  GH','FontName',FN,'FontSize',FS,'FontWeight','bold','VerticalAlignment','top');

fprintf('\n  Codim-2 points (R0 , 1/delta):\n');
fprintf('    BT : R0 = %.5f , delta = %.5f , 1/delta = %.4f\n', BT(1),1/BT(2),BT(2));
fprintf('    GH : R0 = %.5f , delta = %.5f , 1/delta = %.4f\n', GH(1),1/GH(2),GH(2));
fprintf('    CP : R0 = %.5f , delta = %.5f , 1/delta = %.4f\n', CP(1),1/CP(2),CP(2));

xlim(ax, R0lim); ylim(ax, Ylim);
xticks(ax, 1 : 0.2 : 2);
yticks(ax, 10 : 10 : 40);
xlabel(ax, 'R_0', 'Interpreter', 'tex');
ylabel(ax, 'avg. immune period (1/\delta, weeks)', 'Interpreter', 'tex');

legend(ax, 'Interpreter', 'tex', 'Location', 'northwest', ...
       'FontSize', FS, 'Box', 'off');

hold(ax, 'off');
end
% ======================================================================

function [a1,a2,a3,S] = rh_coeffs(I,dl,be,f,fp,mu,gamma)
% Characteristic-polynomial coefficients of the endemic Jacobian.
S  = (gamma+mu+f(I))./be;
a1 = I.*be + I.*fp(I) + dl + 2*mu;
a2 = I.^2.*be.*fp(I) + I.*S.*be.^2 + I.*be.*dl + I.*be.*mu ...
     + I.*dl.*fp(I) + 2*I.*fp(I).*mu + dl.*mu + mu.^2;
a3 = I.*( I.*be.*dl.*fp(I) + I.*be.*fp(I).*mu + S.*be.^2.*dl + S.*be.^2.*mu ...
     - be.*dl.*gamma + dl.*fp(I).*mu + fp(I).*mu.^2 );
end

function r = bt_res(x,beta_eq,f,fp,mu,gamma)   % BT residual [a2; a3]
be = beta_eq(x(1),x(2));
[~,a2,a3] = rh_coeffs(x(1),x(2),be,f,fp,mu,gamma);
r = [a2; a3];
end

function r = gh_res(x,beta_eq,f,fp,fpp,fppp,mu,gamma)
be = beta_eq(x(1),x(2));
[a1,a2,a3] = rh_coeffs(x(1),x(2),be,f,fp,mu,gamma);
r = [a1*a2-a3; l1_3d(x(1),x(2),beta_eq,f,fp,fpp,fppp,mu,gamma)];
end

function l1 = l1_3d(I,dl,beta_eq,f,fp,fpp,fppp,mu,gamma)
% First Lyapunov coefficient at a Hopf point of the 3-D SIRS field via
% Kuznetsov's multilinear formula. Nonlinearity (s,i,r about equilibrium):
%   quadratic P=-beta*s*i ; Q=beta*s*i-(hpp/2)i^2 ; cubic Q=-(hppp/6)i^3
% with h(I)=I*f(I), hpp=2f'+I f'', hppp=3f''+I f'''.  l1<0 supercrit.
be = beta_eq(I,dl); S = (gamma+mu+f(I))/be;
hpp = 2*fp(I)+I*fpp(I);  hppp = 3*fpp(I)+I*fppp(I);
A = [ -(be*I+mu), -be*S,  dl ; ...
       be*I,      -fp(I)*I, 0 ; ...
       0,          gamma,  -(mu+dl) ];
[V,D] = eig(A); ev = diag(D);
[~,k] = min(abs(real(ev)) + (abs(imag(ev))<1e-9)*1e3);  w = abs(imag(ev(k)));
if w < 1e-12, l1 = NaN; return; end
q = V(:,k);
[W,Dl] = eig(A.'); evl = diag(Dl);
[~,kl] = min(abs(evl - conj(ev(k)))); p = W(:,kl);
p = p / conj(p'*q);
Bf = @(x,y) [ -be*(x(1)*y(2)+x(2)*y(1)); ...
               be*(x(1)*y(2)+x(2)*y(1)) - hpp*x(2)*y(2); 0 ];
Cf = @(x,y,z) [ 0; -hppp*x(2)*y(2)*z(2); 0 ];
qb = conj(q); E = eye(3);
g21 = p'*Cf(q,q,qb) - 2*(p'*Bf(q, A\Bf(q,qb))) + p'*Bf(qb, (2i*w*E - A)\Bf(q,q));
l1 = real(g21)/(2*w);
end

function x = newton2(res,x0)
% Damped 2-D Newton with finite-difference Jacobian (toolbox-free).
x = x0(:);
for it = 1:60
    F = res(x); if any(~isfinite(F)), break; end
    h = max(1e-9, 1e-6*abs(x)); J = zeros(2,2);
    for j = 1:2
        xp = x; xp(j) = xp(j)+h(j); J(:,j) = (res(xp)-F)/h(j);
    end
    dx = -J\F; t = 1;
    while t > 1e-4
        Fn = res(x+t*dx);
        if all(isfinite(Fn)) && norm(Fn) < norm(F), break; end
        t = t/2;
    end
    x = x + t*dx;
    if norm(t*dx) < 1e-13, break; end
end
end

function P = parse_contourc(C)
% Parse contourc output matrix into a cell array of [2 x n] polylines.
P = {}; col = 1; k = 1;
while col < size(C,2)
    n = C(2,col);
    P{k} = C(:, col+1:col+n);   %#ok<AGROW>
    k = k + 1; col = col + n + 1;
end
end

function plot_mapped_contour(Iv,dv,Z,color,style,lw,beta_eq,R0_of,name)
% Extract zero contour of Z in (Ihat,delta), map to (R0,1/delta), draw lines
% with the given colour and line style under a single legend entry.
P = parse_contourc(contourc(Iv,dv,Z,[0 0]));
first = true;
for k = 1:numel(P)
    Ip = P{k}(1,:); Dp = P{k}(2,:);
    r0 = R0_of(beta_eq(Ip,Dp)); yy = 1./Dp;
    if first
        plot(r0,yy,style,'Color',color,'LineWidth',lw,'DisplayName',name); first=false;
    else
        plot(r0,yy,style,'Color',color,'LineWidth',lw,'HandleVisibility','off');
    end
end
if first   % nothing drawn: still register legend entry
    plot(NaN,NaN,style,'Color',color,'LineWidth',lw,'DisplayName',name);
end
end

function [firstSub,firstSup] = draw_crit(ax,x,y,issub,gg,lw,firstSub,firstSup)
% Draw one Hopf line coloured by criticality. The first subcritical and the
% first supercritical line drawn carry the legend entry; later ones are
% hidden. Updates and returns the first-seen flags.
if numel(x) < 2, return; end
if issub
    col = gg(3,:); name = 'Hopf (subcritical)';   show = firstSub; firstSub = false;
else
    col = gg(4,:); name = 'Hopf (supercritical)'; show = firstSup; firstSup = false;
end
if show
    plot(ax,x,y,'-','Color',col,'LineWidth',lw,'DisplayName',name);
else
    plot(ax,x,y,'-','Color',col,'LineWidth',lw,'HandleVisibility','off');
end
end

function H = trace_hopf(Iv, drange, beta_eq, f, fp, mu, gamma)
% Trace the Hopf locus as the single-valued graph delta = delta_Hopf(Ihat).
% A coarse pass over Ihat brackets the active range (where a physical Hopf
% root exists); a fine pass then solves the 1-D root delta(Ihat) on that range.
% Returns an ordered [2 x N] = [Ihat; delta] polyline -- one continuous curve.
solve = @(I) hopf_delta(I, drange, beta_eq, f, fp, mu, gamma);
Ic  = linspace(Iv(1), Iv(end), 4000);
act = ~isnan(arrayfun(solve, Ic));
if ~any(act), H = zeros(2,0); return; end
i0 = find(act,1,'first'); i1 = find(act,1,'last');
If = linspace(Ic(max(i0-1,1)), Ic(min(i1+1,numel(Ic))), 4000);
D  = arrayfun(solve, If);
ok = ~isnan(D);
H  = [If(ok); D(ok)];
end

function ds = hopf_delta(I, drange, beta_eq, f, fp, mu, gamma)
% For a fixed Ihat, return the unique delta on the genuine Hopf locus
% (a1*a2-a3 = 0 with a1>0, a2>0 and a physical equilibrium), or NaN if none.
dd = linspace(drange(1), drange(2), 2000).';
be = beta_eq(I,dd);
[a1,a2,a3,S] = rh_coeffs(I,dd,be,f,fp,mu,gamma);
Hv = a1.*a2 - a3;
Hv(~((be>0) & (S>0) & isfinite(be) & a1>0 & a2>0)) = NaN;
sg = sign(Hv); ix = find(sg(1:end-1).*sg(2:end) < 0, 1, 'first');
if isempty(ix), ds = NaN; return; end
ds = fzero(@(x) hopf_scalar(x,beta_eq,f,fp,mu,gamma), [dd(ix) dd(ix+1)]);
    function H = hopf_scalar(dl,beta_eq,f,fp,mu,gamma)
        be0 = beta_eq(I,dl);
        [b1,b2,b3] = rh_coeffs(I,dl,be0,f,fp,mu,gamma);
        H = b1*b2 - b3;
    end
end
