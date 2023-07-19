% AUTHOR: MAJD YOUSOF
% CIVE50003 COMPUTATIONAL METHODS II

classdef BRIDGE

% This class is used to construct an instance of a section of a truss
% bridge, in which deflections can be represented visually thanks to
% FEA.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONTENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 0. PROPERTIES: 
% Properties of the class BRIDGE, which represents the bridge truss
% structure from A to D
%
% 1. SELF INSTANTIATION: 
% Creates an object of that is the instance of the class BRIDGE, with 
% various inputs taken.
%
% 2. PREPROCESSOR MODULE: 
% Initialises the nodes, elements and nodal degrees of freedom (dofs) 
% within the object.
%
% 3. PREPROCESSOR MODULE WITH O1-O2 ADDED:  
% Initialises the nodes, elements and nodal degrees of freedom (dofs) 
% within the object, including the horizontal bar O1-O2.
%
% 4. ASSEMBLER MODULE: 
% Assembles the global stiffness matrix (k) of the truss structure.
%
% 5. SOLVER MODULE: 
% Solves for the nodal displacements, utilising Cholesky decomposition for
% matrix kFF.
%
% 6. POSTPROCESSOR MODULE
% Visualises the bridge deflections as a MATLAB plot. Checks the for force
% equillibrium to ensure calculations are accurate. ENSURE THAT A FIGURE IS
% INITIALISED FIRST!!! And a title would be nice aswell :)
%
% 7. CHOLESKY COMPARISON MODULE
% Compares the speed that the cholesky decomposition method and the regular
% method for matrix inversion through an allocated number of iterations. It
% also compares the condition and sparsity of the 2 matrices used in the
% calculations. It takes in the amount of iterations as a parameter, and
% uses the load location allocated in the Solver module.
%
% 8. SAFETY FACTOR CALCULATIONS
% Creates arrays corresponding to Pcr, Ps and their respective safety
% factors, whilst filtering out non compression members for Pcr and
% filtering out non tension elements for Ps. Postprocessor must be ran
% beforehand.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0. PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        E % Youngs Modulus [N/mm^2].
        A % Cross-sectional area of I-section [mm^2].
        P % Magnitude of point load (train) [N].
        L % Length of horizontal bar element [mm].
        fy % Yield strength of steel [N/mm^2].
        I % Second moment of area [mm^4].
        loadloc % Location of train load with respect to dofs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OTHER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
        coords % Struct containing arrays relating to old and new
        % node coordinates and number of nodes.
        dofs % Struct containing arrays related to degrees of freedom.
        elements % Elements in the bridge.
        k % Struct containing arrays related to the global stiffness matrix.
        f % Nodal force vector for each dof.
        reactions % 92x1 Force vector containing reaction forces.
        d % 92x1 Displacement vector for nodal displacements.
        Faxial % Array containing the axial force of each element in order.
        maxdifx % Maximum x deflection.
        maxdify % Maximum y deflection.
        safety % Struct containing all variables relating to safety factor.

    end
    
    methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 1. INSTANTIATION MODULE %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        function obj = BRIDGE(E,A,P,L,fy,I)

            obj.E = E;
            obj.A = A;
            obj.P = P;
            obj.L = L;
            obj.fy = fy;
            obj.I = I;

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 2. PREPROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = Preprocessor(obj)
            
            fprintf('\nPreprocessor working...')
            
            % Specifying nodal x-y coordinates, left to right in each
            % section

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.coords.all = [0,    0 ;  % NODE 1 - A, BRIDGE DECK START
                          obj.L,    0 ;  % NODE 2  
                          obj.L*2,  0 ;  % NODE 3  
                          obj.L*3,  0 ;  % NODE 4  
                          obj.L*4,  0 ;  % NODE 5 
                          obj.L*5,  0 ;  % NODE 6 
                          obj.L*6,  0 ;  % NODE 7 
                          obj.L*7,  0 ;  % NODE 8 
                          obj.L*8,  0 ;  % NODE 9 
                          obj.L*9,  0 ;  % NODE 10 
                          obj.L*10, 0 ;  % NODE 11 - B 
                          obj.L*11, 0 ;  % NODE 12 
                          obj.L*12, 0 ;  % NODE 13 
                          obj.L*13, 0 ;  % NODE 14 
                          obj.L*14, 0 ;  % NODE 15 
                          obj.L*15, 0 ;  % NODE 16 - C
                          obj.L*16, 0 ;  % NODE 17 
                          obj.L*17, 0 ;  % NODE 18 
                          obj.L*18, 0 ;  % NODE 19 
                          obj.L*19, 0 ;  % NODE 20 
                          obj.L*20, 0 ;  % NODE 21 
                          obj.L*21, 0 ;  % NODE 22 
                          obj.L*22, 0 ;  % NODE 23 
                          obj.L*23, 0 ;  % NODE 24 - D, BRIDGE DECK END
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          obj.L,     9100 ;  % NODE 25 - TOP OF TRUSS START
                          obj.L*2,   9778 ;  % NODE 26 
                          obj.L*3,  10456 ;  % NODE 27 
                          obj.L*4,  11133 ;  % NODE 28 
                          obj.L*5,  11811 ;  % NODE 29 
                          obj.L*6,  12489 ;  % NODE 30 
                          obj.L*7,  13167 ;  % NODE 31 
                          obj.L*8,  13844 ;  % NODE 32 
                          obj.L*9,  14522 ;  % NODE 33 
                          obj.L*10, 15200 ;  % NODE 34 
                          obj.L*11, 13978 ;  % NODE 35 
                          obj.L*12, 12756 ;  % NODE 36 
                          obj.L*13, 11544 ;  % NODE 37 
                          obj.L*14, 10322 ;  % NODE 38 
                          obj.L*15,  9100 ;  % NODE 39 - O1 
                          obj.L*16,  9100 ;  % NODE 40 - O2
                          obj.L*17,  9867 ;  % NODE 41 
                          obj.L*18, 10633 ;  % NODE 42 
                          obj.L*19, 11400 ;  % NODE 43 
                          obj.L*20, 10633 ;  % NODE 44 
                          obj.L*21,  9867 ;  % NODE 45 
                          obj.L*22,  9100];  % NODE 46 - TOP OF TRUSS END 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Initialising degrees of freedom
            obj.dofs.all = zeros(length(obj.coords.all),2);
            for i = 1:length(obj.coords.all)
                obj.dofs.all(i,:) = [2*i-1, 2*i];
            end
                 
            % Specifying element nodal connectivity (order does not matter)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.elements.all = [1,2 ; % ELEMENT 1 - B1, BRIDGE DECK START
                        2,   3 ; % ELEMENT 2 - B2
                        3,   4 ; % ELEMENT 3 - B3
                        4,   5 ; % ELEMENT 4 - B4
                        5,   6 ; % ELEMENT 5 - B5
                        6,   7 ; % ELEMENT 6 - B6
                        7,   8 ; % ELEMENT 7 - B7
                        8,   9 ; % ELEMENT 8 - B8
                        9,  10 ; % ELEMENT 9 - B9
                        10, 11 ; % ELEMENT 10 - B10
                        11, 12 ; % ELEMENT 11
                        12, 13 ; % ELEMENT 12
                        13, 14 ; % ELEMENT 13
                        14, 15 ; % ELEMENT 14
                        15, 16 ; % ELEMENT 15
                        16, 17 ; % ELEMENT 16
                        17, 18 ; % ELEMENT 17
                        18, 19 ; % ELEMENT 18
                        19, 20 ; % ELEMENT 19
                        20, 21 ; % ELEMENT 20
                        21, 22 ; % ELEMENT 21
                        22, 23 ; % ELEMENT 22
                        23, 24 ; % ELEMENT 23 - BRIDGE DECK END
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        2,  25 ; % ELEMENT 24 - VERTICAL ELEMENTS START
                        3,  26 ; % ELEMENT 25
                        4,  27 ; % ELEMENT 26
                        5,  28 ; % ELEMENT 27
                        6,  29 ; % ELEMENT 28
                        7,  30 ; % ELEMENT 29
                        8,  31 ; % ELEMENT 30
                        9,  32 ; % ELEMENT 31
                        10, 33 ; % ELEMENT 32
                        11, 34 ; % ELEMENT 33
                        12, 35 ; % ELEMENT 34
                        13, 36 ; % ELEMENT 35
                        14, 37 ; % ELEMENT 36
                        15, 38 ; % ELEMENT 37
                        16, 39 ; % ELEMENT 38
                        17, 40 ; % ELEMENT 39
                        18, 41 ; % ELEMENT 40
                        19, 42 ; % ELEMENT 41
                        20, 43 ; % ELEMENT 42
                        21, 44 ; % ELEMENT 43
                        22, 45 ; % ELEMENT 44
                        23, 46 ; % ELEMENT 45 - VERTICAL ELEMENTS END
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        1,  25 ; % ELEMENT 46 - DIAGONAL ELEMENTS START
                        2,  26 ; % ELEMENT 47
                        3,  27 ; % ELEMENT 48
                        4,  28 ; % ELEMENT 49
                        5,  29 ; % ELEMENT 50
                        6,  30 ; % ELEMENT 51
                        7,  31 ; % ELEMENT 52
                        8,  32 ; % ELEMENT 53
                        9,  33 ; % ELEMENT 54
                        10, 34 ; % ELEMENT 55
                        12, 34 ; % ELEMENT 56
                        13, 35 ; % ELEMENT 57
                        14, 36 ; % ELEMENT 58
                        15, 37 ; % ELEMENT 59
                        16, 38 ; % ELEMENT 60
                        16, 40 ; % ELEMENT 61
                        18, 40 ; % ELEMENT 62
                        19, 41 ; % ELEMENT 63
                        20, 42 ; % ELEMENT 64
                        20, 44 ; % ELEMENT 65
                        21, 45 ; % ELEMENT 66
                        22, 46 ; % ELEMENT 67
                        24, 46 ; % ELEMENT 68 - DIAGONAL ELEMENTS END
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        25, 26 ; % ELEMENT 69 - T2, TOP OF TRUSS START
                        26, 27 ; % ELEMENT 70 - T3
                        27, 28 ; % ELEMENT 71 - T4
                        28, 29 ; % ELEMENT 72 - T5
                        29, 30 ; % ELEMENT 73 - T6
                        30, 31 ; % ELEMENT 74 - T7
                        31, 32 ; % ELEMENT 75 - T8
                        32, 33 ; % ELEMENT 76 - T9
                        33, 34 ; % ELEMENT 77 - T10
                        34, 35 ; % ELEMENT 78 
                        35, 36 ; % ELEMENT 79
                        36, 37 ; % ELEMENT 80
                        37, 38 ; % ELEMENT 81
                        38, 39 ; % ELEMENT 82
                        40, 41 ; % ELEMENT 83
                        41, 42 ; % ELEMENT 84
                        42, 43 ; % ELEMENT 85
                        43, 44 ; % ELEMENT 86
                        44, 45 ; % ELEMENT 87
                        45, 46]; % ELEMENT 88 - TOP OF TRUSS END
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Setting the restrained and unrestrained degrees of freedom
            obj.dofs.free = 1:obj.dofs.all(end);
            obj.dofs.free(obj.dofs.free == 1) = [];
            obj.dofs.free(obj.dofs.free == 2) = [];
            obj.dofs.free(obj.dofs.free == 22) = [];
            obj.dofs.free(obj.dofs.free == 48) = [];
            obj.dofs.restrained = [1,2,22,48];
    
            obj.coords.numofnodes = size(obj.coords.all,1);
            obj.elements.numofelements = size(obj.elements.all,1);

            fprintf('\nPreprocessing OK!\n');
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% 3. PREPROCESSOR MODULE WITH O1-O2 ADDED %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = Preprocessor_Extra(obj)

            fprintf('\nPreprocessor working...\n')
           
            % Specifying nodal x-y coordinates, left to right in each
            % section

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.coords.all = [0,    0 ;  % NODE 1 - A, BRIDGE DECK START
                          obj.L,    0 ;  % NODE 2  
                          obj.L*2,  0 ;  % NODE 3  
                          obj.L*3,  0 ;  % NODE 4  
                          obj.L*4,  0 ;  % NODE 5 
                          obj.L*5,  0 ;  % NODE 6 
                          obj.L*6,  0 ;  % NODE 7 
                          obj.L*7,  0 ;  % NODE 8 
                          obj.L*8,  0 ;  % NODE 9 
                          obj.L*9,  0 ;  % NODE 10 
                          obj.L*10, 0 ;  % NODE 11 - B 
                          obj.L*11, 0 ;  % NODE 12 
                          obj.L*12, 0 ;  % NODE 13 
                          obj.L*13, 0 ;  % NODE 14 
                          obj.L*14, 0 ;  % NODE 15 
                          obj.L*15, 0 ;  % NODE 16 - C
                          obj.L*16, 0 ;  % NODE 17 
                          obj.L*17, 0 ;  % NODE 18 
                          obj.L*18, 0 ;  % NODE 19 
                          obj.L*19, 0 ;  % NODE 20 
                          obj.L*20, 0 ;  % NODE 21 
                          obj.L*21, 0 ;  % NODE 22 
                          obj.L*22, 0 ;  % NODE 23 
                          obj.L*23, 0 ;  % NODE 24 - D, BRIDGE DECK END
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          obj.L,     9100 ;  % NODE 25 - TOP OF TRUSS START
                          obj.L*2,   9778 ;  % NODE 26 
                          obj.L*3,  10456 ;  % NODE 27 
                          obj.L*4,  11133 ;  % NODE 28 
                          obj.L*5,  11811 ;  % NODE 29 
                          obj.L*6,  12489 ;  % NODE 30 
                          obj.L*7,  13167 ;  % NODE 31 
                          obj.L*8,  13844 ;  % NODE 32 
                          obj.L*9,  14522 ;  % NODE 33 
                          obj.L*10, 15200 ;  % NODE 34 
                          obj.L*11, 13978 ;  % NODE 35 
                          obj.L*12, 12756 ;  % NODE 36 
                          obj.L*13, 11544 ;  % NODE 37 
                          obj.L*14, 10322 ;  % NODE 38 
                          obj.L*15,  9100 ;  % NODE 39 - O1 
                          obj.L*16,  9100 ;  % NODE 40 - O2
                          obj.L*17,  9867 ;  % NODE 41 
                          obj.L*18, 10633 ;  % NODE 42 
                          obj.L*19, 11400 ;  % NODE 43 
                          obj.L*20, 10633 ;  % NODE 44 
                          obj.L*21,  9867 ;  % NODE 45 
                          obj.L*22,  9100];  % NODE 46 - TOP OF TRUSS END 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Initialising degrees of freedom
            obj.dofs.all = zeros(length(obj.coords.all),2);
            for i = 1:length(obj.coords.all)
                obj.dofs.all(i,:) = [2*i-1, 2*i];
            end
                 
            % Specifying element nodal connectivity (order does not matter)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.elements.all = [1,2 ; % ELEMENT 1 - B1, BRIDGE DECK START 
                        2,   3 ; % ELEMENT 2 - B2
                        3,   4 ; % ELEMENT 3 - B3
                        4,   5 ; % ELEMENT 4 - B4
                        5,   6 ; % ELEMENT 5 - B5
                        6,   7 ; % ELEMENT 6 - B6
                        7,   8 ; % ELEMENT 7 - B7
                        8,   9 ; % ELEMENT 8 - B8
                        9,  10 ; % ELEMENT 9 - B9
                        10, 11 ; % ELEMENT 10 - B10
                        11, 12 ; % ELEMENT 11
                        12, 13 ; % ELEMENT 12
                        13, 14 ; % ELEMENT 13
                        14, 15 ; % ELEMENT 14
                        15, 16 ; % ELEMENT 15
                        16, 17 ; % ELEMENT 16
                        17, 18 ; % ELEMENT 17
                        18, 19 ; % ELEMENT 18
                        19, 20 ; % ELEMENT 19
                        20, 21 ; % ELEMENT 20
                        21, 22 ; % ELEMENT 21
                        22, 23 ; % ELEMENT 22
                        23, 24 ; % ELEMENT 23 - BRIDGE DECK END 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        2,  25 ; % ELEMENT 24 - VERTICAL ELEMENTS START 
                        3,  26 ; % ELEMENT 25
                        4,  27 ; % ELEMENT 26
                        5,  28 ; % ELEMENT 27
                        6,  29 ; % ELEMENT 28
                        7,  30 ; % ELEMENT 29
                        8,  31 ; % ELEMENT 30
                        9,  32 ; % ELEMENT 31
                        10, 33 ; % ELEMENT 32
                        11, 34 ; % ELEMENT 33
                        12, 35 ; % ELEMENT 34
                        13, 36 ; % ELEMENT 35
                        14, 37 ; % ELEMENT 36
                        15, 38 ; % ELEMENT 37
                        16, 39 ; % ELEMENT 38
                        17, 40 ; % ELEMENT 39
                        18, 41 ; % ELEMENT 40
                        19, 42 ; % ELEMENT 41
                        20, 43 ; % ELEMENT 42
                        21, 44 ; % ELEMENT 43
                        22, 45 ; % ELEMENT 44
                        23, 46 ; % ELEMENT 45 - VERTICAL ELEMENTS END 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        1,  25 ; % ELEMENT 46 - DIAGONAL ELEMENTS START
                        2,  26 ; % ELEMENT 47
                        3,  27 ; % ELEMENT 48
                        4,  28 ; % ELEMENT 49
                        5,  29 ; % ELEMENT 50
                        6,  30 ; % ELEMENT 51
                        7,  31 ; % ELEMENT 52
                        8,  32 ; % ELEMENT 53
                        9,  33 ; % ELEMENT 54
                        10, 34 ; % ELEMENT 55
                        12, 34 ; % ELEMENT 56
                        13, 35 ; % ELEMENT 57
                        14, 36 ; % ELEMENT 58
                        15, 37 ; % ELEMENT 59
                        16, 38 ; % ELEMENT 60
                        16, 40 ; % ELEMENT 61
                        18, 40 ; % ELEMENT 62
                        19, 41 ; % ELEMENT 63
                        20, 42 ; % ELEMENT 64
                        20, 44 ; % ELEMENT 65
                        21, 45 ; % ELEMENT 66
                        22, 46 ; % ELEMENT 67
                        24, 46 ; % ELEMENT 68 - DIAGONAL ELEMENTS END
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        25, 26 ; % ELEMENT 69 - T2, TOP OF TRUSS START 
                        26, 27 ; % ELEMENT 70 - T3
                        27, 28 ; % ELEMENT 71 - T4
                        28, 29 ; % ELEMENT 72 - T5
                        29, 30 ; % ELEMENT 73 - T6
                        30, 31 ; % ELEMENT 74 - T7
                        31, 32 ; % ELEMENT 75 - T8
                        32, 33 ; % ELEMENT 76 - T9
                        33, 34 ; % ELEMENT 77 - T10
                        34, 35 ; % ELEMENT 78 
                        35, 36 ; % ELEMENT 79
                        36, 37 ; % ELEMENT 80
                        37, 38 ; % ELEMENT 81
                        38, 39 ; % ELEMENT 82
                        40, 41 ; % ELEMENT 83
                        41, 42 ; % ELEMENT 84
                        42, 43 ; % ELEMENT 85
                        43, 44 ; % ELEMENT 86
                        44, 45 ; % ELEMENT 87
                        45, 46 ; % ELEMENT 88 - TOP OF TRUSS END 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        39, 40]; % ELEMENT 89 - O1-O2 ADDED
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Setting the restrained and unrestrained degrees of freedom
            obj.dofs.free = 1:obj.dofs.all(end);
            obj.dofs.free(obj.dofs.free == 1) = [];
            obj.dofs.free(obj.dofs.free == 2) = [];
            obj.dofs.free(obj.dofs.free == 22) = [];
            obj.dofs.free(obj.dofs.free == 48) = [];
            obj.dofs.restrained = [1,2,22,48];
    
            obj.coords.numofnodes = size(obj.coords.all,1);
            obj.elements.numofelements = size(obj.elements.all,1);

            fprintf('Preprocessing OK!\n');
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% 4. ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = Assembler(obj)
           
            EA = obj.E*obj.A;

            % Initialising stiffness matrix and element length array
            obj.k.kfull = zeros(2*obj.coords.numofnodes);
            obj.elements.length = zeros(obj.elements.numofelements,1);

            % Looping through all elements & build stiffness matrix
            for EL = 1:obj.elements.numofelements

                % Setting element nodal coordinates
                n1 = obj.elements.all(EL,1); n2 = obj.elements.all(EL,2); 
                x1 = obj.coords.all(n1,1); y1 = obj.coords.all(n1,2); % node 1
                x2 = obj.coords.all(n2,1); y2 = obj.coords.all(n2,2); % node 2 
                
                % Creating element nodal dofs array 
                dof11 = obj.dofs.all(n1,1); dof12 = obj.dofs.all(n1,2);
                dof21 = obj.dofs.all(n2,1); dof22 = obj.dofs.all(n2,2);
                eldofs = [dof11 dof12 dof21 dof22]; 
                
                % angle of inclination relative to the positive x axis
                alpha = atan2(y2-y1,x2-x1); 
                c = cos(alpha); s = sin(alpha); % angle parameters
                
                % element length
                obj.elements.length(EL) = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); 
                
                % element axial stiffness
                ke = EA/obj.elements.length(EL); 
                
                % Define 2x4 transformation matrix [T]
                T = [c s 0 0; 0 0 c s];
                
                % Local 2x2 stiffness matrix
                kp = ke*[1 -1; -1 1];
                
                % Transformation to stiffness matrix and insertion to 
                % global [k]
                obj.k.kfull(eldofs, eldofs) = obj.k.kfull(eldofs, ...
                    eldofs) + T'*kp*T;
                
                if mod(EL,20)==0
                    fprintf('\nAssembly working...');
                end

            end
                fprintf('\nAssembly OK!\n');
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% 5. SOLVER MODULE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = Solver(obj,loadloc)

        fprintf('\nSolver working...\n')

            % For readability
            restrained = obj.dofs.restrained;
            free = obj.dofs.free;

            obj.loadloc = loadloc;

            fknown = zeros(2*obj.coords.numofnodes,1);
            fknown(obj.loadloc) = obj.P;
            
            % Splitting full stiffness matrix and prepping for calculation
            obj.k.kRR = obj.k.kfull(restrained,restrained);
            obj.k.kRF = obj.k.kfull(restrained,free);
            obj.k.kFR = obj.k.kfull(free,restrained);
            obj.k.kFF = obj.k.kfull(free,free); 
            fF = fknown(free);

            % Applying Cholesky decomposition
            obj.k.cholkFF = chol(obj.k.kFF);
            

            dR = zeros(length(obj.dofs.restrained),1);
            dFequation = (fF - obj.k.kFR*dR);
            dF = obj.k.cholkFF\(obj.k.cholkFF'\dFequation); % 1st matrix
            % equation

            % Creation of nodal displacement vector
            obj.d = zeros(obj.coords.numofnodes,1);
            obj.d(restrained) = dR';
            obj.d(free) = dF';

            fR = obj.k.kRF*dF + obj.k.kRR*dR; % 2nd matrix equation

            % Creation of nodal force vector
            obj.f = zeros(obj.coords.numofnodes,1);
            obj.f(restrained) = fR';
            obj.f(free) = fF';

            % Creation of vector containing reaction forces at restricted
            % nodes
            obj.reactions = obj.f;
            obj.reactions(obj.loadloc) = obj.reactions(obj.loadloc)-obj.P;
            
            fprintf('Solving OK!\n')
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 6. POSTPROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = Postprocessor(obj)
            
            % visual amplification factors
            amp = 10;
            FS = 25;

            % Initialising new coordinate matrices
            EA = obj.E*obj.A;
            obj.coords.new = zeros(size(obj.coords.all)); 
            obj.coords.amp = zeros(size(obj.coords.all));
            
            %filling in new matrices with nodal displacement vector
            for K = 1:size(obj.coords.all,1)
                for J = 1:size(obj.coords.all,2)

                    fprintf('Element no. %g\n',K)
                    obj.coords.amp(K,J) = obj.coords.all(K, ...
                        J) + obj.d(obj.dofs.all(K,J))*amp;

                    obj.coords.new(K,J) = obj.coords.all(K, ...
                        J) + obj.d(obj.dofs.all(K,J));
                end
            end

            %IMPORTANT: Do not forget to initialise a figure
            hold all; grid on; 
            
            tol = 1e-3;

            % Calculating maximum x deflection
            xmin = min(obj.coords.amp(:,1));
            xmax = max(obj.coords.amp(:,1));
            obj.maxdifx = xmax - xmin;

            % Calculating maximum y deflection
            ymin = min(obj.coords.amp(:,2));
            ymax = max(obj.coords.amp(:,2));
            obj.maxdify = ymax - ymin;

            obj.Faxial = zeros(obj.elements.numofelements,1);
            axis equal padded;
            
            
            % Plotting elements and nodes
            for EL = 1:obj.elements.numofelements

                n1 = obj.elements.all(EL,1); 
                n2 = obj.elements.all(EL,2);
                
                % Plotting original structure
                x1 = obj.coords.all(n1,1); 
                y1 = obj.coords.all(n1,2);
                x2 = obj.coords.all(n2,1); 
                y2 = obj.coords.all(n2,2);

                ke = EA/obj.elements.length(EL); % element axial stiffness
                alpha = atan2(y2-y1,x2-x1);

                plot([x1,x2],[y1,y2],'Color',[0.5 0.5 0.5],'Linewidth',3); 
                
                % Check changes in member lengths and plotting amplified
                % deformed structure
                x1_amp = obj.coords.amp(n1,1); 
                y1_amp = obj.coords.amp(n1,2); % element node 1 - 
                % x,y amplified deformed coordinates
                
                x2_amp = obj.coords.amp(n2,1); 
                y2_amp = obj.coords.amp(n2,2); % element node 2 -
                % x,y amplified deformed coordinates     
                
                x1_new = obj.coords.new(n1,1); 
                y1_new = obj.coords.new(n1,2); % element node 1 -
                % x,y actual deformed coordinates
                
                x2_new = obj.coords.new(n2,1); 
                y2_new = obj.coords.new(n2,2); % element node 2 -
                % x,y actual deformed coordinates 

                % reconstruction of element global dofs
                ux1 = x1_new - x1;
                uy1 = y1_new - y1;
                ux2 = x2_new - x2; 
                uy2 = y2_new - y2; 

                % reconstruction of element local dofs 
                up1 = cos(alpha)*ux1 + sin(alpha)*uy1;
                up2 = cos(alpha)*ux2 + sin(alpha)*uy2;
                dup = up2 - up1;     
               
                % note that this now gives you access to the element force
                % element length has decreased - member in compression
                if dup < -tol 
                    col = 'b'; % blue colour

                % element length as increased - member in tension    
                elseif dup > tol 
                    col = 'r'; % red colour

                % no change in element length    
                else 
                    col = 'k'; % black colour
                end

                plot([x1_amp,x2_amp],[y1_amp,y2_amp],col,'Linewidth',3);
                
                % Calculate the axial force in the element and store it
                obj.Faxial(EL) = ke*dup;
                
                % Plotting nodes last!
                
                % Commented out for speed, feel free to add back in.
                %plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
                %plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');

                plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
                plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                
            end

            % Marking location of train load using green marker
            if mod(obj.loadloc,2) == 0

                % Even dof
                plot(obj.coords.amp(obj.loadloc/2,1), ...
                    obj.coords.amp(obj.loadloc/2,2), ...
                    'ko','Markersize',10,'MarkerFaceColor','g')

                text(obj.coords.amp(obj.loadloc/2,1), ...
                    obj.coords.amp(obj.loadloc/2,2)-1000, ...
                    'Train','HorizontalAlignment','center','VerticalAlignment', ...
                    'top','FontSize',FS,'Interpreter','latex')
            else

                % Odd dof
                plot(obj.coords.amp((obj.loadloc+1)/2,1), ...
                    obj.coords.amp((obj.loadloc+1)/2,2), ...
                    'ko','Markersize',10,'MarkerFaceColor','g')

                text(obj.coords.amp((obj.loadloc+1)/2,1), ...
                    obj.coords.amp((obj.loadloc+1)/2,2)-1000, ...
                    'Train','HorizontalAlignment','center','VerticalAlignment', ...
                    'top','FontSize',FS,'Interpreter','latex')
            end
            
            % Labels
            xlabel('$x$ coordinate [mm]','interpreter','latex','FontSize',FS);
            ylabel('$y$ coordinate [mm]','interpreter','latex','FontSize',FS);
            set(gca,'ticklabelinterpreter','latex','fontsize',FS, ...
                'tickdir','out');

            % Printing out the values of the dofs of each node
            dF = (obj.d(obj.dofs.free))';
            for dof = 1:length(dF)
                fprintf(' The value of the dof %g is %g\n', ...
                    obj.dofs.free(dof),dF(dof));    
            end
            
            % Printing out reactions for each restricted node
            fR = obj.reactions(obj.dofs.restrained);
            for react = 1:length(fR)
                fprintf('The value of the reaction at dof %g is %g\n', ...
                    obj.dofs.restrained(react),fR(react));
            end

            % Equilibrium checks
            fprintf('\nEQUILIBRIUM CHECK\n')

            % Vertical check
            disp(['Total vertical reactions = ', ...
                num2str(obj.reactions(2)+obj.reactions(22)+obj.reactions(48))]);
            disp('');
            disp(['Total applied vertical loads = ',num2str(obj.P)]);
            disp('');
            if abs(obj.reactions(2)+obj.reactions(22)+obj.reactions(48) + obj.P) < 1e-6
                disp('OK!'); 
            end
            disp('');

            % Horizontal check
            disp(['Total horizontal reactions = ', ...
                num2str(obj.reactions(1))]);
            disp('');
            disp('Total applied horizontal loads = 0');
            disp('');
            if abs(obj.reactions(1)) < 1e-6
                disp('OK!'); 
            end
            disp('');
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 7. CHOLESKY COMPARISON MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = Chol_Comparison(obj,numofiterations)
            
            % Initialisation of certain values for comparison
            dR = zeros(length(obj.dofs.restrained),1);
            fknown = zeros(2*obj.coords.numofnodes,1);
            fknown(obj.loadloc) = obj.P;
            fF = fknown(obj.dofs.free);
            dFequation = (fF - obj.k.kFR*dR);

            % initialising calculation timing arrays
            tregular = zeros(numofiterations,1);
            tcholesky = zeros(numofiterations,1);

            % Filling out calculation time arrays
            for i = 1:numofiterations

                if mod(i,(numofiterations/10))==0
                    fprintf('\nCholesky comparison working...');
                end

                % Timing regular calculation
                tic;
                placeholder = obj.k.kFF\dFequation;
                tregular(i) = toc;

                % Timing Cholesky calculation including decomposition 
                tic;
                cholkFF = chol(obj.k.kFF);
                placeholder = cholkFF\(cholkFF'\dFequation);
                tcholesky(i) = toc;

            end
            
            % Reciprocal condition
            cholcond = rcond(cholkFF);
            regularcond = rcond(obj.k.kFF);
            
            tmeanchol = mean(tcholesky);
            tmeanregular = mean(tregular);

            figure('units','normalized','outerposition',[0 0 1 1]);
            
            % Boxplot comparison between computational timings of both
            % methods
            subplot(2,2,[1,3])
            hold on; grid on;
            boxchart([tregular,tcholesky]);
            plot([1,2],[tmeanregular,tmeanchol],'-or');
            hold off;
            
            set(gca,'Yscale','log')
            set(gca,'ticklabelinterpreter','latex','fontsize',20, ...
                'tickdir','out');
            set(gca,'XTickLabel',{'Regular','Cholesky'});

            legend('Boxplots','Mean time interval [s]')
            xlabel('Calculation method','interpreter','latex','FontSize',25);
            ylabel('Time [s] (Logarithmic)','interpreter','latex','FontSize',25);
            title('Time Comparison between methods','Interpreter','latex','FontSize',30);

            % kFF matrix sparsity visualisation
            subplot(2,2,2)
            spy(obj.k.kFF)
            set(gca,'ticklabelinterpreter','latex','fontsize',20, ...
                'tickdir','out');
            title('Regular kFF matrix sparsity','Interpreter','latex','FontSize',30);
            xlabel(sprintf('Number of nonzero elements: %d', ...
                nnz(obj.k.kFF)),Interpreter="latex",FontSize=20);
            % Cholesky matrix sparsity visualisation
            subplot(2,2,4)
            spy(cholkFF)
            set(gca,'ticklabelinterpreter','latex','fontsize',20, ...
                'tickdir','out');
            xlabel(sprintf('Number of nonzero elements: %d', ...
                nnz(cholkFF)),Interpreter="latex",FontSize=20);
            title('Cholesky decomposition of kFF matrix sparsity','Interpreter','latex','FontSize',30);
            

            % printing results of comparison
            fprintf('\n\nNumber of iterations for comparison: %g\n\n',numofiterations)
            fprintf('Mean time taken for Cholesky method: %g\n',tmeanchol)
            fprintf('Reciprocal Condition number for the Cholesky Matrix: %g\n\n',cholcond)
            fprintf('Mean time taken for Regular method: %g\n',tmeanregular)
            fprintf('Reciprocal Condition number for regular kFF Matrix: %g\n',regularcond)
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 8. SAFETY FACTOR CALCULATIONS %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        function obj = Safety_Factor(obj)

            % This is for the GLOBAL truss rather than only B1-B10 and
            % T2-T10
            
            % Setting the 0 tolerance for compression filter
            tol = 1e-4;

            % Calculating values for Pcr and Ps 
            obj.safety.Pcr = obj.E*obj.I*(pi./obj.elements.length).^2;
            obj.safety.Ps = obj.fy*obj.A;
            
            % Calculating safety factor for each element
            obj.safety.safetyfactor_Pcr = obj.safety.Pcr./obj.Faxial;
            obj.safety.safetyfactor_Ps = obj.safety.Ps./obj.Faxial;

            % Filtering out for only values of B1-B10 and T2-T10
            obj.safety.safetyfactor_Pcr = obj.safety.safetyfactor_Pcr;
            obj.safety.safetyfactor_Ps = obj.safety.safetyfactor_Ps;

            % Filtering for Compression Elements
            obj.safety.elementsPcr = find(obj.safety.safetyfactor_Pcr < -tol);
            obj.safety.safetyfactor_Pcr = abs(obj.safety.safetyfactor_Pcr(obj.safety.elementsPcr));

            obj.safety.elementsPs = find(obj.safety.safetyfactor_Ps < -tol);
            obj.safety.safetyfactor_Ps = abs(obj.safety.safetyfactor_Ps(obj.safety.elementsPs));
            
            % Calculating the minimum (critical) safety factor and which
            % element it corresponds to
            [obj.safety.minPcr, indexPcr] = min(obj.safety.safetyfactor_Pcr);
            obj.safety.minPcr_element = obj.safety.elementsPcr(indexPcr);

            [obj.safety.minPs, indexPs] = min(obj.safety.safetyfactor_Ps);
            obj.safety.minPs_element = obj.safety.elementsPs(indexPs);

        end

    end
end

