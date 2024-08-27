function [areactionwithout, breactionwithout, dreactionwithout, ...
          Baxialwithout, Taxialwithout, ...
          areactionwith, breactionwith, dreactionwith, ...
          Baxialwith, Taxialwith, Pcrwithout, Pswithout, ...
          Pcrwith, Pswith] = animateTrain(bwithout, bwith)
    areactionwithout = zeros(24, 1);
    breactionwithout = zeros(24, 1);
    dreactionwithout = zeros(24, 1);
    areactionwith = zeros(24, 1);
    breactionwith = zeros(24, 1);
    dreactionwith = zeros(24, 1);
    Baxialwithout = zeros(24, 10);
    Baxialwith = zeros(24, 10);
    Taxialwithout = zeros(24, 9);
    Taxialwith = zeros(24, 9);
    Pcrwith = zeros(24, 2);
    Pcrwithout = zeros(24, 2);
    Pswith = zeros(24, 2);
    Pswithout = zeros(24, 2);

    figure('units', 'normalized', 'outerposition', [0 0 1 1], 'Name', 'Train Animations');
    fprintf('\nAnimating...')
    for i = 1:24
        bwithout = bwithout.Solver(2*i);
        bwith = bwith.Solver(2*i);

        clf;
        subplot(2, 1, 1);
        bwithout = bwithout.Postprocessor();
        title('Without O1-O2', 'Interpreter', 'latex', 'FontSize', 30);
        subplot(2, 1, 2);
        bwith = bwith.Postprocessor();
        title('With O1-O2', 'Interpreter', 'latex', 'FontSize', 30);
        hold on;
        drawnow update;

        bwithout = bwithout.Safety_Factor();
        bwith = bwith.Safety_Factor();

        % Reactions
        areactionwithout(i) = bwithout.reactions(2);
        breactionwithout(i) = bwithout.reactions(22);
        dreactionwithout(i) = bwithout.reactions(48);
        areactionwith(i) = bwith.reactions(2);
        breactionwith(i) = bwith.reactions(22);
        dreactionwith(i) = bwith.reactions(48);

        % Safety Factors
        [Pcrwithout(i,:), Pswithout(i,:)] = getSafetyFactors(bwithout);
        [Pcrwith(i,:), Pswith(i,:)] = getSafetyFactors(bwith);

        % Axial Forces
        Baxialwithout(i, :) = bwithout.Faxial(1:10);
        Baxialwith(i, :) = bwith.Faxial(1:10);
        Taxialwithout(i, :) = bwithout.Faxial(69:77);
        Taxialwith(i, :) = bwith.Faxial(69:77);
    end
    hold off;
    fprintf('\nAnimation OK!\n')
end

function [Pcr, Ps] = getSafetyFactors(b)
    if b.safety.minPcr ~= inf
        Pcr = [b.safety.minPcr, b.safety.minPcr_element];
    else
        Pcr = [1e20, 0];
    end
    if b.safety.minPs ~= inf
        Ps = [b.safety.minPs, b.safety.minPs_element];
    else
        Ps = [1e20, 0];
    end
end
