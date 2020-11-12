function truss = exampleTRUSS_tenbar

% structure of the truss
% |o --- o --- o
%     X  |  X  |
% |o --- o --- *

% numbering of the nodes
% |4 --- 5 --- 6
%     X  |  X  |
% |1 --- 2 --- 3


%% truss data

% dimension of the truss
truss.dimension = 2;

% x and y coordinates of all nodes in the given order
truss.nodeCoordinates = [0 0; ... % node 1
                         1 0; ... % node 2
                         2 0; ... % node 3
                         0 1; ... % node 4
                         1 1; ... % node 5
                         2 1];    % node 6

% number of all fixed nodes in the given order
truss.fixedNodes = [1 4];

% number of the start and end node of all potential bars in the given order
truss.potentialBars = [1 2; ... % bar 1
                       1 5; ... % bar 2
                       2 3; ... % bar 3
                       2 4; ... % bar 4
                       2 5; ... % bar 5
                       2 6; ... % bar 6
                       3 5; ... % bar 7
                       3 6; ... % bar 8
                       4 5; ... % bar 9
                       5 6];    % bar 10

% possible load scenarios
% for each load scenario: indicate the force applied to all nodes in the given order
% stack the resulting (n_nodes x dimension)-sized matrices to obtain a tensor
truss.loadCases(:,:,1) = [0  0; ... % force applied in node 1
                          0  0; ... % force applied in node 2
                          0 -1; ... % force applied in node 3
                          0  0; ... % force applied in node 4
                          0  0; ... % force applied in node 5
                          0  0];    % force applied in node 6
truss.loadCases(:,:,2) = [0  0; ... % force applied in node 1
                          0  0; ... % force applied in node 2
                          -0.2 -1; ... % force applied in node 3
                          0  0; ... % force applied in node 4
                          0  0; ... % force applied in node 5
                          0  0];    % force applied in node 6
truss.loadCases(:,:,3) = [0  0; ... % force applied in node 1
                          0  0; ... % force applied in node 2
                          0.2 -1; ... % force applied in node 3
                          0  0; ... % force applied in node 4
                          0  0; ... % force applied in node 5
                          0  0];    % force applied in node 6

% parameters
truss.compliance_max = 100;
truss.displacement_max = inf;
truss.stresses_max = 1;
truss.barDiameters_max = 10;
truss.barDiameters_min = 0; % is enforced only for realized bars
truss.barDiameters_start = 1;


%% options for the MPVC algorithm

options.objectiveGradient = true;
options.constraintsJacobian = true;
options.NLPsolver = 'fmincon';
options.slacks = false;
% options.algorithm = 'relaxation';
% options.relaxation = 'scholtes';


%% optimize truss

figure

options.algorithm = 'direct';
[truss, information] = optimizeTRUSS(truss, options);
maxVio = max([information.maxVio_diameter...
              information.maxVio_displacement...
              information.maxVio_equilibrium...
              information.maxVio_compliance...
              information.maxVio_stress]);
% information.message
disp([options.algorithm, '   ', num2str(truss.volume_opt), '   ', num2str(information.iterations), '   ', num2str(maxVio)])
% truss.barDiameters_opt'

subplot(2,5,1)
title('direct')
hold on
plotTRUSS_unloadedStructure(truss)


options.algorithm = 'relaxation';
relaxations = {'scholtes', 'steffensen', 'schwartz', 'kadrani'};
for relax = 1:4
    options.relaxation = relaxations{relax};
    [truss, information] = optimizeTRUSS(truss, options);
    % information.message
    maxVio = max([information.maxVio_diameter...
                  information.maxVio_displacement...
                  information.maxVio_equilibrium...
                  information.maxVio_compliance...
                  information.maxVio_stress]);
    disp([options.algorithm, '   ', options.relaxation, '   ', num2str(truss.volume_opt), '   ', num2str(information.iterations), '   ', num2str(maxVio)])
    % truss.barDiameters_opt'
    
    subplot(2,5,relax+1)
    title([options.algorithm ' ' options.relaxation])
    hold on
    plotTRUSS_unloadedStructure(truss)
end

options.algorithm = 'relaxation_posLB';
relaxations = {'scholtes', 'steffensen', 'schwartz', 'kadrani'};
for relax = 1:4
    options.relaxation = relaxations{relax};
    [truss, information] = optimizeTRUSS(truss, options);
    % information.message
    maxVio = max([information.maxVio_diameter...
                  information.maxVio_displacement...
                  information.maxVio_equilibrium...
                  information.maxVio_compliance...
                  information.maxVio_stress]);
    disp([options.algorithm, '   ', options.relaxation, '   ', num2str(truss.volume_opt), '   ', num2str(information.iterations), '   ', num2str(maxVio)])
    % truss.barDiameters_opt'
    
    subplot(2,5,relax+6)
    title(['relaxation\_posLB' ' ' options.relaxation])
    hold on
    plotTRUSS_unloadedStructure(truss)
end


%% plot initial structure and optimized structure

% plot initial structure
% figure
% plotTRUSS_initialStructure(truss)

% plot optimized but unloaded truss structure
% plotTRUSS_unloadedStructure(truss)

% plot optimized and loaded truss structure
% plotTRUSS_loadedStructure(truss)
