clear;clc;close all;

% AUTHOR: MAJD YOUSOF
% CIVE50003 COMPUTATIONAL METHODS II

%% Initialising values

% Font sizes
fst = 30;
fs = 25;

E = 2e5; % Youngs Modulus [N/mm^2]
A = 21259.97; % Cross-sectional area of I-section [mm^2]
P = -2000000; % Magnitude of point load (train) [N]
L = 7600; % Length of horizontal bar element [mm]
fy = 355; % Yield strength of steel [N/mm^2]
I = 9.86874e7; % Second moment of area [mm^4]


%% Q1. Influence lines for A, B and D

span = [0,114,174.8,288.8];
span2 = [0,114,174.8];
RAinfluence = [2,-1,0,0];
RBinfluence = [0,3,0,0];
RDinfluence = [0,0,2];

figure('units','normalized','outerposition',[0 0 1 1]);

hold on;
plot([0,288.8],[0,0],'-k',LineWidth=2);
plot(span,RAinfluence,'--b',LineWidth=2);
plot(span,RBinfluence,'--r',LineWidth=2);
plot(span2,RDinfluence,'--g',LineWidth=2);
hold off;

grid on, axis padded;
ylabel('Reaction force [MN]','Interpreter','latex','FontSize',fs);
xlabel('Location across the bridge span [m]','Interpreter','latex','FontSize',fs);
title('Influence line plots of the reactions at A and B from hand calculations','Interpreter','latex','FontSize',fst);
legend('Bridge Deck','Reaction at A','Reaction at B','Reaction at D', ...
    'fontsize',fs,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',fs, ...
                'tickdir','out');

text(0,0,'A','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(76,0,'B','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(114,0,'C','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(174.8,0,'D','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(212.8,0,'E','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(288.8,0,'F','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')

%% Q3, Q5, Q6, Q8. Train animation plot and initialising influence line arrays 

% Calculation prep
b = BRIDGE(E,A,P,L,fy,I);
bwithout = b.Preprocessor();
bwith = b.Preprocessor_Extra();
bwithout = bwithout.Assembler();
bwith = bwith.Assembler();

% Reaction vectors without O1-O2 initialisation
areactionwithout = zeros(24,1);
breactionwithout = zeros(24,1);
dreactionwithout = zeros(24,1);

% Reaction vectors with O1-O2 initialisation
areactionwith = zeros(24,1);
breactionwith = zeros(24,1);
dreactionwith = zeros(24,1);

% Axial force of B1-B10 vectors without O1-O2 initialisation 
Baxialwithout = zeros(24,10);

% Axial force of B1-B10 vectors with O1-O2 initialisation 
Baxialwith = zeros(24,10);

% Axial force of T1-T10 vectors without O1-O2 initialisation 
Taxialwithout = zeros(24,9);

% Axial force of B1-B10 vectors with O1-O2 initialisation 
Taxialwith = zeros(24,9);

% Safety factor vector initialisation
Pcrwith = zeros(24,2);
Pcrwithout = zeros(24,2);
Pswith = zeros(24,2);
Pswithout = zeros(24,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANIMATION START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1],'Name','Train Animations');


% Would be probably more efficient to presolve everything to increase the
% speed of animation, but this will do for now :)
fprintf('\nAnimating...')
for i = 1:24
    
    % Solving FEM
    bwithout = bwithout.Solver(2*i);
    bwith = bwith.Solver(2*i);

    clf;

    % Animation
    subplot(2,1,1);
    bwithout = bwithout.Postprocessor();
    title('Without O1-O2','Interpreter','latex','FontSize',fst)
    subplot(2,1,2);
    bwith = bwith.Postprocessor();
    title('With O1-O2','Interpreter','latex','FontSize',fst)

    hold on;
    drawnow update;

    % Calculating safety factors
    bwithout = bwithout.Safety_Factor();
    bwith = bwith.Safety_Factor();

    % The code below looks like alot, however it is just data retrieval
    % from the class, so doesn't affect speed much

    % Adding reactions for plot
    areactionwithout(i) = bwithout.reactions(2);
    breactionwithout(i) = bwithout.reactions(22);
    dreactionwithout(i) = bwithout.reactions(48);

    areactionwith(i) = bwith.reactions(2);
    breactionwith(i) = bwith.reactions(22);
    dreactionwith(i) = bwith.reactions(48);

    % Adding safety factor values to their respective vectors
    if bwithout.safety.minPcr ~= inf
        Pcrwithout(i,1) = bwithout.safety.minPcr;
        Pcrwithout(i,2) = bwithout.safety.minPcr_element;
    else
        %indicates that all elements have a safety factor of infinity
        Pcrwithout(i,1) = 1e20;
        Pcrwithout(i,2) = 0;
    end

    if bwith.safety.minPcr ~= inf
        Pcrwith(i,1) = bwith.safety.minPcr;
        Pcrwith(i,2) = bwith.safety.minPcr_element;
    else
        %indicates that all elements have a safety factor of infinity
        Pcrwith(i,1) = 1e20;
        Pcrwith(i,2) = 0;
    end

    if bwithout.safety.minPs ~= inf
        Pswithout(i,1) = bwithout.safety.minPs;
        Pswithout(i,2) = bwithout.safety.minPs_element;
    else
        %indicates that all elements have a safety factor of infinity
        Pswithout(i,1) = 1e20;
        Pswithout(i,2) = 0;
    end

    if bwith.safety.minPs ~= inf
        Pswith(i,1) = bwith.safety.minPs;
        Pswith(i,2) = bwith.safety.minPs_element;
    else
        %indicates that all elements have a safety factor of infinity
        Pswith(i,1) = 1e20;
        Pswith(i,2) = 0;
    end

    for b1 = 1:10
        % Adding axial forces B1-B10 for plot
        Baxialwithout(i,b1) = bwithout.Faxial(b1);
        Baxialwith(i,b1) = bwith.Faxial(b1);
    end

    for b2 = 1:9
        % Adding axial forces T1-T10 for plot
        Taxialwithout(i,b2) = bwithout.Faxial(68+b2);
        Taxialwith(i,b2) = bwith.Faxial(68+b2);
    end   
    
end

hold off;
fprintf('\nAnimation OK!\n')

% The next few sections are influence line plots, so ensure you run this
% section BEFORE the plot sections. Cholesky decomposition can be ran
% independent of this however.

%%
%%%%%% Plotting influence lines for reactions WITHOUT O1-O2 INCLUDED %%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);

hold on;
plot([0,288.8],[0,0],'-k',LineWidth=2);
plot(linspace(0,174.800,24),areactionwithout,'--b',LineWidth=2);
plot(linspace(0,174.800,24),breactionwithout,'--r',LineWidth=2);
plot(linspace(0,174.800,24),dreactionwithout,'--g',LineWidth=2);
hold off;

grid on, axis padded;
ylabel('Reaction force [N]','Interpreter','latex','FontSize',fs);
xlabel('Location across the bridge span [m]','Interpreter','latex','FontSize',fs);
title('Influence line plots of the reactions at A and B from FEA without O1-O2','Interpreter','latex','FontSize',fst);
legend('Bridge Deck','Reaction at A','Reaction at B','Reaction at D', ...
    'fontsize',fs,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',fs, ...
                'tickdir','out');

text(0,0,'A','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(76,0,'B','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(114,0,'C','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(174.8,0,'D','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(212.8,0,'E','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(288.8,0,'F','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')

%%
%%%%%%%% Plotting influence lines for reactions WITH O1-O2 INCLUDED %%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);

hold on;
plot([0,288.8],[0,0],'-k',LineWidth=2);
plot(linspace(0,174.800,24),areactionwith,'--b',LineWidth=2);
plot(linspace(0,174.800,24),breactionwith,'--r',LineWidth=2);
plot(linspace(0,174.800,24),dreactionwith,'--g',LineWidth=2);
hold off;

grid on, axis padded;
ylabel('Reaction force [N]','Interpreter','latex','FontSize',fs);
xlabel('Location across the bridge span [m]','Interpreter','latex','FontSize',fs);
title('Influence line plots of the reactions at A and B from FEA with O1-O2 included','Interpreter','latex','FontSize',fst);
legend('Bridge Deck','Reaction at A','Reaction at B','Reaction at D', ...
    'fontsize',fs,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',fs, ...
                'tickdir','out');

text(0,0,'A','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(76,0,'B','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(114,0,'C','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(174.8,0,'D','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(212.8,0,'E','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(288.8,0,'F','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')

%%
%%%%%% Influence lines for B1-B10 axial forces WITHOUT O1-O2 INCLUDED %%%%%
figure('units','normalized','outerposition',[0 0 1 1]);

hold on;

% Bridge deck
plot([0,288.8],[0,0],'-k',LineWidth=2);

% Plots
for i = 1:10
plot(linspace(0,174.800,24),Baxialwithout(:,i),'-',LineWidth=2);
end

hold off;

grid on, axis padded;
ylabel('Axial force [N]','Interpreter','latex','FontSize',fs);
xlabel('Length across the bridge span [m]','Interpreter','latex','FontSize',fs);
title('Influence lines of the axial forces in B1-B10 WITHOUT O1-O2','Interpreter','latex','FontSize',fst);
legend('Bridge Deck','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10', ...
    'fontsize',fs,'interpreter','latex','location','bestoutside');
set(gca,'ticklabelinterpreter','latex','fontsize',fs, ...
                'tickdir','out');

text(0,0,'A','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(76,0,'B','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(114,0,'C','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(174.8,0,'D','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(212.8,0,'E','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(288.8,0,'F','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')

%%
%%%%%%% Influence lines for B1-B10 axial forces WITH O1-O2 INCLUDED %%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);

hold on;

% Bridge deck
plot([0,288.8],[0,0],'-k',LineWidth=2);

% Plots
for i = 1:10
    plot(linspace(0,174.800,24),Baxialwith(:,i),'-',LineWidth=2);
end

hold off;

grid on, axis padded;
ylabel('Axial force [N]','Interpreter','latex','FontSize',fs);
xlabel('Location across the bridge span [m]','Interpreter','latex','FontSize',fs);
title('Influence lines of the axial forces in B1-B10 WITH O1-O2 included','Interpreter','latex','FontSize',fst);
legend('Bridge Deck','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10', ...
    'fontsize',fs,'interpreter','latex','location','bestoutside');
set(gca,'ticklabelinterpreter','latex','fontsize',fs, ...
                'tickdir','out');

text(0,0,'A','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(76,0,'B','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(114,0,'C','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(174.8,0,'D','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(212.8,0,'E','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(288.8,0,'F','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')

%%
%%%%%% Influence lines for T1-T10 axial forces WITHOUT O1-O2 INCLUDED %%%%%
figure('units','normalized','outerposition',[0 0 1 1]);

hold on;

% Bridge deck
plot([0,288.8],[0,0],'-k',LineWidth=2);

% Plots
for i = 1:9
    plot(linspace(0,174.800,24),Taxialwithout(:,i),'-',LineWidth=2);
end

hold off;

grid on, axis padded;
ylabel('Axial force [N]','Interpreter','latex','FontSize',fs);
xlabel('Location across the bridge span [m]','Interpreter','latex','FontSize',fs);
title('Influence lines of the axial forces in T1-T10 WITHOUT O1-O2','Interpreter','latex','FontSize',fst);
legend('Bridge Deck','T2','T3','T4','T5','T6','T7','T8','T9','T10', ...
    'fontsize',fs,'interpreter','latex','location','bestoutside');
set(gca,'ticklabelinterpreter','latex','fontsize',fs, ...
                'tickdir','out');

text(0,0,'A','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(76,0,'B','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(114,0,'C','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(174.8,0,'D','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(212.8,0,'E','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(288.8,0,'F','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')

%%
%%%%%%% Influence lines for T1-T10 axial forces WITH O1-O2 INCLUDED %%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);

hold on;

% Bridge deck
plot([0,288.8],[0,0],'-k',LineWidth=2);

% Plots
for i = 1:9
    plot(linspace(0,174.800,24),Taxialwith(:,i),'-',LineWidth=2);
end

hold off;

grid on, axis padded;
ylabel('Axial force [N]','Interpreter','latex','FontSize',fs);
xlabel('Location across the bridge span [m]','Interpreter','latex','FontSize',fs);
title('Influence lines of the axial forces in T1-T10 WITH O1-O2 included','Interpreter','latex','FontSize',fst);
legend('Bridge Deck','T2','T3','T4','T5','T6','T7','T8','T9','T10', ...
    'fontsize',fs,'interpreter','latex','location','bestoutside');
set(gca,'ticklabelinterpreter','latex','fontsize',fs, ...
                'tickdir','out');

text(0,0,'A','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(76,0,'B','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(114,0,'C','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(174.8,0,'D','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(212.8,0,'E','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')
text(288.8,0,'F','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs,'Interpreter','latex')

%% Q7. Safety factor

% OPTIONAL: adjustable safety factor to see when your bridge fails and why,
% used by for loops commented out.

safetyfactor = nan; % adjust to factor to a number you want as a threshold

fprintf('\nWithout O1-O2:\n');

% Loop can be deployed when a specific safety factor is desired to be
% acheived when traversing the truss bridge to ensure failiure doesn't
% occur

% for i = 1:length(Pcrwithout)
% 
%     if Pcrwithout(i,1) < safetyfactor
%         fprintf('Buckling failiure:\nThe bridge fails through element %g\n with a safety factor of %g\n when the train is in position %g\n', ...
%             Pcrwithout(i,2),Pcrwithout(i,1),i);
%         break
%     elseif Pswithout < safetyfactor
%         fprintf('Squash failiure:\nThe bridge fails through element %g\n with a safety factor of %g\n when the train is in position %g\n', ...
%             Pswithout(i,2),Pswithout(i,1),i);
%         break
%     end
% 
% end

% calculating the absolute minimum safety factors of the bridge
[minPcrwithout, indexPcrwithout] = min(Pcrwithout(:,1));
fprintf('\nBuckling: The minimum safety factor if the bridge is %g\nand is achieved by member %g\nwhen the train is in location %g\n', ...
    minPcrwithout,Pcrwithout(indexPcrwithout,2),indexPcrwithout);

[minPswithout, indexPswithout] = min(Pswithout(:,1));
fprintf('\nSquash: The minimum safety factor if the bridge is %g\nand is achieved by member %g\nwhen the train is in location %g\n ', ...
    minPswithout,Pswithout(indexPswithout,2),indexPswithout);

fprintf('\nWith O1-O2:\n');

% Loop can be deployed when a specific safety factor is desired to be
% acheived when traversing the truss bridge to ensure failiure doesn't
% occur

% for i = 1:length(Pcrwith)
% 
%     if Pcrwith(i,1) < safetyfactor
%         fprintf('Buckling failiure:\nThe bridge fails through element %g\n with a safety factor of %g\n when the train is in position %g\n', ...
%             Pcrwith(i,2),Pcrwith(i,1),i);
%         break
%     elseif Pswith < safetyfactor
%         fprintf('Squash failiure:\nThe bridge fails through element %g\n with a safety factor of %g\n when the train is in position %g\n', ...
%             Pswith(i,2),Pswith(i,1),i);
%         break
%     end
%     
% end

% calculating the absolute minimum safety factors of the bridge
[minPcrwith, indexPcrwith] = min(Pcrwith(:,1));
fprintf('\nBuckling: The minimum safety factor if the bridge is %g\nand is achieved by member %g\nwhen the train is in location %g\n', ...
    minPcrwith,Pcrwith(indexPcrwith,2),indexPcrwith);

[minPswith, indexPswith] = min(Pswith(:,1));
fprintf('\nSquash: The minimum safety factor if the bridge is %g\nand is achieved by member %g\nwhen the train is in location %g\n', ...
    minPswith,Pswith(indexPswith,2),indexPswith);

%% Q3. Cholesky method comparison

b = BRIDGE(E,A,P,L,fy,I);
b = b.Preprocessor();
b = b.Assembler();
b = b.Solver(10); % Feel free to change, numbers from 1-92


figure('units','normalized','outerposition',[0 0 1 1]);
b = b.Postprocessor();
title('Truss Deformation','Interpreter','latex','FontSize',fst);

hold off;

b = b.Safety_Factor();

% Comparison start
b = b.Chol_Comparison(10000); 
% You can change the amount of iterations as you like.

