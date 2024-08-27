function plotInfluenceLines()
    span = [0, 114, 174.8, 288.8];
    span2 = [0, 114, 174.8];
    RAinfluence = [2, -1, 0, 0];
    RBinfluence = [0, 3, 0, 0];
    RDinfluence = [0, 0, 2];
    fs = 25;
    fst = 30;

    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on;
    plot([0, 288.8], [0, 0], '-k', 'LineWidth', 2);
    plot(span, RAinfluence, '--b', 'LineWidth', 2);
    plot(span, RBinfluence, '--r', 'LineWidth', 2);
    plot(span2, RDinfluence, '--g', 'LineWidth', 2);
    hold off;

    grid on; axis padded;
    ylabel('Reaction force [MN]', 'Interpreter', 'latex', 'FontSize', fs);
    xlabel('Location across the bridge span [m]', 'Interpreter', 'latex', 'FontSize', fs);
    title('Influence line plots of the reactions at A and B from hand calculations', 'Interpreter', 'latex', 'FontSize', fst);
    legend('Bridge Deck', 'Reaction at A', 'Reaction at B', 'Reaction at D', ...
           'fontsize', fs, 'interpreter', 'latex');
    set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', fs, ...
             'tickdir', 'out');

    text(0, 0, 'A', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', fs, 'Interpreter', 'latex');
    text(76, 0, 'B', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', fs, 'Interpreter', 'latex');
    text(114, 0, 'C', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', fs, 'Interpreter', 'latex');
    text(174.8, 0, 'D', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', fs, 'Interpreter', 'latex');
    text(212.8, 0, 'E', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', fs, 'Interpreter', 'latex');
    text(288.8, 0, 'F', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', fs, 'Interpreter', 'latex');
end
