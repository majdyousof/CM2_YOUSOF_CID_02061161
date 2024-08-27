function plotResults(areactionwithout, breactionwithout, dreactionwithout, ...
                      areactionwith, breactionwith, dreactionwith, ...
                      Baxialwithout, Baxialwith, Taxialwithout, Taxialwith)
    % Plot influence lines for reactions without O1-O2
    plotInfluenceLinesForReactions(areactionwithout, breactionwithout, dreactionwithout, ...
                                    'without O1-O2', 'N');
    
    % Plot influence lines for reactions with O1-O2
    plotInfluenceLinesForReactions(areactionwith, breactionwith, dreactionwith, ...
                                    'with O1-O2', 'N');

    % Plot axial forces
    plotAxialForces(Baxialwithout, Taxialwithout, 'without O1-O2', 'N');
    plotAxialForces(Baxialwith, Taxialwith, 'with O1-O2', 'N');
end

function plotInfluenceLinesForReactions(A, B, D, titleText, unit)
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on;
    plot(1:24, A, '--b', 'LineWidth', 2);
    plot(1:24, B, '--r', 'LineWidth', 2);
    plot(1:24, D, '--g', 'LineWidth', 2);
    hold off;
    
    grid on;
    xlabel('Position of Train [m]', 'Interpreter', 'latex', 'FontSize', 25);
    ylabel(['Reaction force ', unit], 'Interpreter', 'latex', 'FontSize', 25);
    title(['Influence line plots of the reactions ', titleText], 'Interpreter', 'latex', 'FontSize', 30);
    legend('Reaction at A', 'Reaction at B', 'Reaction at D', ...
           'fontsize', 25, 'interpreter', 'latex');
    set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 25, ...
             'tickdir', 'out');
end

function plotAxialForces(Baxial, Taxial, titleText, unit)
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on;
    plot(1:24, Baxial, 'LineWidth', 2);
    plot(1:24, Taxial, 'LineWidth', 2);
    hold off;
    
    grid on;
    xlabel('Position of Train [m]', 'Interpreter', 'latex', 'FontSize', 25);
    ylabel(['Axial Force ', unit], 'Interpreter', 'latex', 'FontSize', 25);
    title(['Axial Forces ', titleText], 'Interpreter', 'latex', 'FontSize', 30);
    legend('Axial Forces - Without O1-O2', 'Axial Forces - With O1-O2', ...
           'fontsize', 25, 'interpreter', 'latex');
    set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 25, ...
             'tickdir', 'out');
end
