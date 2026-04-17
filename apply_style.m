%% apply_style.m
% Apply ggplot2-style formatting to the current open figure and save as SVG.
%
% Usage: run this script while the target figure is the current figure (gcf).
%   >> apply_style
%
% Requires set_ggplot_style.m to be on the MATLAB path.

set_ggplot_style();          % refresh groot defaults (affects future objects)
fig = gcf;
fig.Color = [1 1 1];

%% ── Axes ─────────────────────────────────────────────────────────────────

all_axes = findall(fig, 'Type', 'axes');

for k = 1:numel(all_axes)
    ax = all_axes(k);

    ax.Box                  = 'off';
    ax.TickDir              = 'out';
    ax.FontName             = 'Helvetica';
    ax.FontSize             = 10;
    ax.LineWidth            = 0.8;
    ax.XColor               = [0.3 0.3 0.3];
    ax.YColor               = [0.3 0.3 0.3];
    ax.ZColor               = [0.3 0.3 0.3];
    ax.XGrid                = 'off';
    ax.YGrid                = 'off';
    ax.TickLabelInterpreter = 'tex';

    % Remove first y-tick mark and label (ggplot2 origin gap)
    yt = ax.YTick;
    if numel(yt) > 1,  ax.YTick = yt(2:end);  end
end

%% ── Lines ────────────────────────────────────────────────────────────────

all_lines = findall(fig, 'Type', 'line');
for k = 1:numel(all_lines)
    all_lines(k).LineWidth = 1.2;
end

%% ── Text ─────────────────────────────────────────────────────────────────

all_text = findall(fig, 'Type', 'text');
for k = 1:numel(all_text)
    all_text(k).FontName = 'Helvetica';
    all_text(k).FontSize = 10;
end

%% ── Legends ──────────────────────────────────────────────────────────────

all_legends = findall(fig, 'Type', 'legend');
for k = 1:numel(all_legends)
    lg             = all_legends(k);
    lg.Box         = 'off';
    lg.FontName    = 'Helvetica';
    lg.FontSize    = 10;
    lg.Interpreter = 'tex';
end

%% ── Save as SVG ──────────────────────────────────────────────────────────

if ~isempty(strtrim(fig.Name))
    fname = strtrim(fig.Name);
else
    fname = sprintf('figure_%d', fig.Number);
end
fname = regexprep(fname, '\.svg$', '', 'ignorecase');  % avoid double extension

print(fig, '-dsvg', fname);
fprintf('Saved -> %s.svg\n', fname);
