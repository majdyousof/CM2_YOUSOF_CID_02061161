function evaluateSafetyFactors(Pcrwithout, Pswithout, Pcrwith, Pswith)
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on;
    plot(1:24, Pcrwithout(:, 1), '--r', 'LineWidth', 2);
    plot(1:24, Pswithout(:, 1), '--b', 'LineWidth', 2);
    plot(1:24, Pcrwith(:, 1), '--g', 'LineWidth', 2);
    plot(1:24, Pswith(:, 1), '--k', 'LineWidth', 2);
    hold off;

    grid on;
    xlabel('Position of Train [m]', 'Interpreter', 'latex', 'FontSize', 25);
    ylabel('Safety Factor', 'Interpreter', 'latex', 'FontSize', 25);
    title('Safety Factors for Buckling and Yielding', 'Interpreter', 'latex', 'FontSize', 30);
    legend('Pcr without O1-O2', 'Ps without O1-O2', ...
           'Pcr with O1-O2', 'Ps with O1-O2', ...
           'fontsize', 25, 'interpreter', 'latex');
    set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 25, ...
             'tickdir', 'out');
end
