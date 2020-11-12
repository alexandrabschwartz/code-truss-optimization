function truss = exampleTRUSS_cantilever(width, height, n_horizontalNodes, n_verticalNodes)

% The input parameters for this function to optimize a rectangular truss
% fixed to a wall on the left side are
    % width              = width of the truss
    % height             = height of the truss
    % n_horizontalNodes  = number of nodes along the width of the truss
    % n_verticalNodes    = number of nodes along the height of the truss

    
% structure of the truss
% |o --- o --- o --- o
%     X  |  X  |  X  |
% |o --- o --- o --- o
%     X  |  X  |  X  |
% |o --- o --- o --- o
%     X  |  X  |  X  |
% |o --- o --- o --- *


% example: numbering of the nodes in case n_horizontal = 4, n_vertical = 3
%  9 --- 10 --- 11 --- 12
%     X   |  X   |  X  |
%  5 ---  6 ---  7 --- 8
%     X   |  X   |  X  |
%  1 ---  2 ---  3 --- 4

%% truss data

% dimension of the truss
truss.dimension = 2;

% x and y coordinates of all nodes in the given order
n_nodes = n_horizontalNodes * n_verticalNodes;

distance_x = width/(n_horizontalNodes-1); % horizontal distance of the nodes
distance_y = height/(n_verticalNodes-1); % vertical distance of the nodes
nodeCoordinates(:,1) = repmat((0:n_horizontalNodes-1)'*distance_x, n_verticalNodes, 1);
nodeCoordinates(:,2) = reshape(repmat((0:n_verticalNodes-1)*distance_y, n_horizontalNodes,1), n_nodes, 1);

truss.nodeCoordinates = nodeCoordinates;

% all nodes with node_x = 0 are fixed
truss.fixedNodes = find(nodeCoordinates(:,1) == 0);

% indices of the start and end nodes of all potential bars in the truss
truss = setupTRUSS_computePotentialBars(truss);

% load case: force pulling down at the bottom right corner of the truss
n_loadCases = 3;
loadCases = zeros(n_nodes, truss.dimension, n_loadCases);
loadPoint = find((nodeCoordinates(:,1) == width) & (nodeCoordinates(:,2) == 0));
loadCases(loadPoint, :, 1) = [0    -1];
loadCases(loadPoint, :, 2) = [0.2  -1];
loadCases(loadPoint, :, 2) = [-0.2 -1];
truss.loadCases = loadCases;

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