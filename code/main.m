clear; clc; close all;

% Constants
FONTSIZE_LARGE = 30;
FONTSIZE_MEDIUM = 25;
E = 2e5; % Young's Modulus [N/mm^2]
A = 21259.97; % Cross-sectional area of I-section [mm^2]
P = -2000000; % Magnitude of point load (train) [N]
L = 7600; % Length of horizontal bar element [mm]
fy = 355; % Yield strength of steel [N/mm^2]
I = 9.86874e7; % Second moment of area [mm^4]

% Influence Lines
plotInfluenceLines();

% Initialize Bridge Models
[bwithout, bwith] = initialiseBridgeModels(E, A, P, L, fy, I);

% Train Animation
[areactionwithout, breactionwithout, dreactionwithout, ...
 Baxialwithout, Taxialwithout, ...
 areactionwith, breactionwith, dreactionwith, ...
 Baxialwith, Taxialwith, Pcrwithout, Pswithout, ...
 Pcrwith, Pswith] = animateTrain(bwithout, bwith);

% Plot Results
plotResults(areactionwithout, breactionwithout, dreactionwithout, ...
            areactionwith, breactionwith, dreactionwith, ...
            Baxialwithout, Baxialwith, Taxialwithout, Taxialwith);

% Safety Factor
evaluateSafetyFactors(Pcrwithout, Pswithout, Pcrwith, Pswith);

% Cholesky Method Comparison
compareCholeskyMethod(E, A, P, L, fy, I);
