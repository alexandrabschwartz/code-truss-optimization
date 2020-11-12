function truss = exampleTRUSS_hook(size)

% The input parameter for this function is the width = height of the
% truss. Within these dimensions, nodes are spaced with distance 1 and
% the upper right quadrant of the truss is cut out.

% structure of the truss
%  _     _     _    
%  o --- o --- o  
%     X  |  X  | 
%  o --- o --- o
%     X  |  X  |   
%  o --- o --- o     size = 5
%     X  |  X  |   
%  o --- o --- o --- o --- o --- o 
%     X  |  X  |  X  |  X  |  X  | 
%  o --- o --- o --- o --- o --- o 
%     X  |  X  |  X  |  X  |  X  | 
%  o --- o --- o --- o --- o --- *

%  _     _     _     _
%  o --- o --- o --- o 
%     X  |  X  |  X  |
%  o --- o --- o --- o 
%     X  |  X  |  X  | 
%  o --- o --- o --- o    size = 6
%     X  |  X  |  X  |  
%  o --- o --- o --- o --- o --- o --- o
%     X  |  X  |  X  |  X  |  X  |  X  |
%  o --- o --- o --- o --- o --- o --- o
%     X  |  X  |  X  |  X  |  X  |  X  |
%  o --- o --- o --- o --- o --- o --- o
%     X  |  X  |  X  |  X  |  X  |  X  |
%  o --- o --- o --- o --- o --- o --- *


%% truss data

% dimension of the truss
truss.dimension = 2;

% x and y coordinates of all nodes
width = size;
height = size;
n_horizontalNodes = size + 1;
n_verticalNodes = size + 1;
n_nodes = n_horizontalNodes * n_verticalNodes;
distance_x = 1; % horizontal distance of the nodes
distance_y = 1; % vertical distance of the nodes
nodeCoordinates(:,1) = repmat((0:n_horizontalNodes-1)'*distance_x, n_verticalNodes, 1);
nodeCoordinates(:,2) = reshape(repmat((0:n_verticalNodes-1)*distance_y, n_horizontalNodes,1), n_nodes, 1);
nodeCoordinates((nodeCoordinates(:,1) > size/2) & (nodeCoordinates(:,2) > size/2),:) = [];

truss.nodeCoordinates = nodeCoordinates;
n_nodes = length(nodeCoordinates);

% indices of those nodes, which are fixed to the ceiling
truss.fixedNodes = find(nodeCoordinates(:,2) == height);

% indices of the start and end nodes of all potential bars in the truss
truss = setupTRUSS_computePotentialBars(truss);
n_bars = length(truss.potentialBars);
nodes_lowerRight = find((nodeCoordinates(:,1) > size/2) & (nodeCoordinates(:,2) <= size/2));
nodes_upperLeft = find((nodeCoordinates(:,1) <= size/2) & (nodeCoordinates(:,2) > size/2));
potentialBars = zeros(0,2);
for bar = 1:n_bars
    startNode = truss.potentialBars(bar,1);
    endNode = truss.potentialBars(bar,2);
    if ismember(startNode, nodes_upperLeft) && ismember(endNode, nodes_lowerRight)
    elseif ismember(endNode, nodes_upperLeft) && ismember(startNode, nodes_lowerRight)
    else
       potentialBars = [potentialBars; startNode endNode];
    end
end
truss.potentialBars = potentialBars;


% possible load scenarios
% for each load scenario: indicate the force applied to all nodes in "node order"
% stack the resulting (n_nodes x dimension)-sized matrices on top of each other
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
truss.stresses_max = 100;
truss.barDiameters_max = 1;
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


% options.algorithm = 'relaxation';
% relaxations = {'scholtes', 'steffensen', 'schwartz', 'kadrani'};
% for relax = 1:4
%     options.relaxation = relaxations{relax};
%     [truss, information] = optimizeTRUSS(truss, options);
%     % information.message
%     maxVio = max([information.maxVio_diameter...
%                   information.maxVio_displacement...
%                   information.maxVio_equilibrium...
%                   information.maxVio_compliance...
%                   information.maxVio_stress]);
%     disp([options.algorithm, '   ', options.relaxation, '   ', num2str(truss.volume_opt), '   ', num2str(information.iterations), '   ', num2str(maxVio)])
%     % truss.barDiameters_opt'
%     
%     subplot(2,5,relax+1)
%     title([options.algorithm ' ' options.relaxation])
%     hold on
%     plotTRUSS_unloadedStructure(truss)
% end
% 
% options.algorithm = 'relaxation_posLB';
% relaxations = {'scholtes', 'steffensen', 'schwartz', 'kadrani'};
% for relax = 1:4
%     options.relaxation = relaxations{relax};
%     [truss, information] = optimizeTRUSS(truss, options);
%     % information.message
%     maxVio = max([information.maxVio_diameter...
%                   information.maxVio_displacement...
%                   information.maxVio_equilibrium...
%                   information.maxVio_compliance...
%                   information.maxVio_stress]);
%     disp([options.algorithm, '   ', options.relaxation, '   ', num2str(truss.volume_opt), '   ', num2str(information.iterations), '   ', num2str(maxVio)])
%     % truss.barDiameters_opt'
%     
%     subplot(2,5,relax+6)
%     title(['relaxation\_posLB' ' ' options.relaxation])
%     hold on
%     plotTRUSS_unloadedStructure(truss)
% end


%% plot initial structure and optimized structure

% plot initial structure
% figure
% plotTRUSS_initialStructure(truss)

% plot optimized but unloaded truss structure
% plotTRUSS_unloadedStructure(truss)

% plot optimized and loaded truss structure
% plotTRUSS_loadedStructure(truss)