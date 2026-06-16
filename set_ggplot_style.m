%% ggplot2-style figure settings for MATLAB
% Matches the aesthetic of the attached plots

function set_ggplot_style()
    % --- Colour palette (matches the lines in the figure) ---
    % Use these in order for multi-line plots
    colours = [
        0.929, 0.490, 0.533;   % pink/salmon   (e.g. beta=1.5 or a=0)
        0.698, 0.718, 0.149;   % olive-green   (e.g. beta=2   or a=10)
        0.196, 0.784, 0.745;   % teal-cyan     (e.g. beta=2.5 or a=100)
        0.557, 0.392, 0.820;   % purple        (e.g. beta=5   or a=1000)
        0.902, 0.498, 0.000;   % warm orange   (e.g. a=10000 or 5th group)
    ];

    % --- Default axes properties ---
    set(groot, 'defaultAxesFontName',        'Helvetica');
    set(groot, 'defaultAxesFontSize',        10);
    set(groot, 'defaultAxesFontWeight',      'normal');
    set(groot, 'defaultAxesColor',           [1 1 1]);          % white background
    set(groot, 'defaultAxesBox',             'off');             % full box
    set(groot, 'defaultAxesLineWidth',       1.5);
    set(groot, 'defaultAxesXColor',          [0.3 0.3 0.3]);
    set(groot, 'defaultAxesYColor',          [0.3 0.3 0.3]);
    set(groot, 'defaultAxesGridColor',       [0.85 0.85 0.85]);
    set(groot, 'defaultAxesGridLineStyle',   '-');
    set(groot, 'defaultAxesGridAlpha',       1);
    set(groot, 'defaultAxesMinorGridAlpha',  0);
    set(groot, 'defaultAxesXGrid',           'off');
    set(groot, 'defaultAxesYGrid',           'off');
    set(groot, 'defaultAxesTickDir',         'out');
    set(groot, 'defaultAxesTickLength',      [0.01 0.01]);

    % --- Default line properties ---
    set(groot, 'defaultLineLineWidth',       1.2);

    % --- Default figure properties ---
    set(groot, 'defaultFigureColor',         [1 1 1]);
    set(groot, 'defaultFigurePaperUnits',    'inches');

    % --- Default text properties ---
    set(groot, 'defaultTextFontName',        'Helvetica');
    set(groot, 'defaultTextFontSize',        10);

    % Store colours in base workspace for convenience
    assignin('base', 'gg_colours', colours);
    disp('ggplot2 style applied. Colour palette stored in ''gg_colours''.');
end


%% Helper: format a single axes to match the style
function format_axes(ax, panel_label, strip_label)
    % ax          - axes handle
    % panel_label - e.g. 'A', 'B', 'C' (pass '' to skip)
    % strip_label - grey header bar text (pass '' to skip)

    if nargin < 2, panel_label = ''; end
    if nargin < 3, strip_label = ''; end

    ax.FontName   = 'Helvetica';
    ax.FontSize   = 10;
    ax.Box        = 'off';
    ax.LineWidth  = 0.6;
    ax.TickDir    = 'out';
    ax.Color      = [1 1 1];
    ax.XColor     = [0.3 0.3 0.3];
    ax.YColor     = [0.3 0.3 0.3];
    ax.XGrid      = 'off';
    ax.YGrid      = 'off';

    % Panel letter (top-left, bold)
    if ~isempty(panel_label)
        xl = ax.XLim; yl = ax.YLim;
        text(ax, xl(1), yl(2), panel_label, ...
            'FontName', 'Helvetica', 'FontSize', 12, ...
            'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
    end

    % Strip label (grey header bar above axes, like ggplot facet strips)
    if ~isempty(strip_label)
        ax_pos = ax.Position;                          % [x y w h] normalised
        strip_h = 0.025;                               % strip height
        annotation(ax.Parent, 'rectangle', ...
            [ax_pos(1), ax_pos(2)+ax_pos(4), ax_pos(3), strip_h], ...
            'FaceColor', [0.85 0.85 0.85], 'EdgeColor', [0.7 0.7 0.7]);
        annotation(ax.Parent, 'textbox', ...
            [ax_pos(1), ax_pos(2)+ax_pos(4), ax_pos(3), strip_h], ...
            'String', strip_label, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontName', 'Helvetica', ...
            'FontSize', 9, 'EdgeColor', 'none', 'BackgroundColor', 'none');
    end
end


%% Helper: add a clean legend matching the style
function h = add_legend(ax, labels, title_str)
    % labels    - cell array of strings
    % title_str - legend title (e.g. '\beta' or 'a')

    h = legend(ax, labels, ...
        'Location',  'eastoutside', ...
        'Box',       'off', ...
        'FontName',  'Helvetica', ...
        'FontSize',  9, ...
        'Title',     title_str);
    h.Title.FontWeight = 'normal';
    h.Title.FontSize   = 9;
end


%% --- Minimal working example ---
% Demonstrates the style on a simple two-panel figure

function demo()
    set_ggplot_style();
    gg = evalin('base', 'gg_colours');

    fig = figure('Units', 'inches', 'Position', [1 1 7 3.5]);

    x = linspace(0, 1000, 500);

    % Panel A
    ax1 = subplot(1, 2, 1);
    hold(ax1, 'on');
    for k = 1:4
        plot(ax1, x, (k*0.05) * (x ./ (x + 100*(5-k)+1)), ...
            'Color', gg(k,:), 'LineWidth', 1.2);
    end
    xlabel(ax1, 'a');
    ylabel(ax1, 'Per-infection mortality rate f(I)');
    format_axes(ax1, 'A', '');

    % Panel B
    ax2 = subplot(1, 2, 2);
    hold(ax2, 'on');
    for k = 1:4
        plot(ax2, x, 1 + (k-1)*0.1 * (x ./ (x + 200)), ...
            'Color', gg(k,:), 'LineWidth', 1.2);
    end
    xlabel(ax2, 'a');
    ylabel(ax2, 'Scaled pop-level mortality');
    format_axes(ax2, 'B', '1/\delta = 0.25 years');
    add_legend(ax2, {'\beta = 1.5', '\beta = 2', '\beta = 2.5', '\beta = 5'}, '\beta');

    % Tight layout
    set(fig, 'PaperPositionMode', 'auto');
end