function truss = exampleTRUSS_fivebar

% structure of the truss
% |o --- o 
%     X  |
% |o --- *

% numbering of the nodes
% |3 --- 4 
%     X  |  
% |1 --- 2 


%% truss data

% dimension of the truss
truss.dimension = 2;

% x and y coordinates of all nodes in the given order
truss.nodeCoordinates = [0 0; ... % node 1
                         1 0; ... % node 2
                         0 1; ... % node 3
                         1 1];    % node 4

% number of all fixed nodes in the given order
truss.fixedNodes = [1 3];

% number of the start and end node of all potential bars in the given order
truss.potentialBars = [1 2; ... % bar 1
                       1 4; ... % bar 2
                       2 3; ... % bar 3
                       2 4; ... % bar 4
                       3 4];    % bar 5

% possible load scenarios
% for each load scenario: indicate the force applied to all nodes in the given order
% stack the resulting (n_nodes x dimension)-sized matrices to obtain a tensor
truss.loadCases(:,:,1) = [0     0; ... % force applied in node 1
                          0    -1; ... % force applied in node 2
                          0     0; ... % force applied in node 3
                          0     0];    % force applied in node 4
truss.loadCases(:,:,2) = [0     0; ... % force applied in node 1
                          0.2  -1; ... % force applied in node 2
                          0     0; ... % force applied in node 3
                          0     0];    % force applied in node 4
truss.loadCases(:,:,3) = [0     0; ... % force applied in node 1
                          -0.2 -1; ... % force applied in node 2
                          0     0; ... % force applied in node 3
                          0     0];    % force applied in node 4

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


%% optimize truss usung different algorithms

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
    %information.message
    maxVio = max([information.maxVio_diameter...
                  information.maxVio_displacement...
                  information.maxVio_equilibrium...
                  information.maxVio_compliance...
                  information.maxVio_stress]);
    disp([options.algorithm, '   ', options.relaxation, '   ', num2str(truss.volume_opt), '   ', num2str(information.iterations), '   ', num2str(maxVio)])
    %truss.barDiameters_opt'
    
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