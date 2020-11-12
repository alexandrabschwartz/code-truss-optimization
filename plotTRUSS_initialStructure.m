function plotTRUSS_initialStructure(truss)

% This function plots the initial unloaded truss.
% The function needs a truss struct with the fields
    % truss.dimension
    % truss.nodeCoordinates
    % truss.fixedNodes
    % truss.potentialBars
    % truss.barDiameters_start (optional)
% If no value for barDiameters_start is provided, a diameter of 1 is used for all bars instead.
    

%% parameters

tol_barDiameter = 10^-4;
scale_maxDiameter = 10;

%% read truss data

if ~isfield(truss, 'dimension') || isempty(truss.dimension)
    error('Plot Groundstructure: Truss dimension (= 2 or 3) is missing')
else 
    dimension = truss.dimension;
end

if ~isfield(truss, 'nodeCoordinates') || isempty(truss.nodeCoordinates)
    error('Plot Groundstructure: Node coordinates missing')
else 
    nodeCoordinates = truss.nodeCoordinates;
    n_nodes = size(nodeCoordinates, 1);
end

if ~isfield(truss, 'fixedNodes') || isempty(truss.fixedNodes)
    error('Plot Groundstructure: Fixed nodes are missing')
else
    fixedNodes = truss.fixedNodes;
end

if ~isfield(truss, 'potentialBars') || isempty(truss.potentialBars)
    disp('Plot Groundstructure: Potential bars are missing, using all possible bars instead')
    truss.potentialBars = TRUSS_potentialBars(truss);
end
potentialBars = truss.potentialBars;
n_bars = size(potentialBars,1);

        
if ~isfield(truss,'barDiameters_start') || isempty(truss.barDiameters_start)
    % no initial values for the diameters are provided, use default value 1 instead
    disp('Plot Groundstructure: Initial diameters for the bars are missing, using a diameter of 1 for all possible bars instead')
    barDiameters = ones(n_bars,1);
else
    % use initial diameters
    barDiameters = truss.barDiameters_start;
    if length(barDiameters) == 1
        barDiameters = barDiameters*ones(n_bars,1);
    end
end


%% determine extremal values

minXCoordinate = min(nodeCoordinates(:,1));
maxXCoordinate = max(nodeCoordinates(:,1));
minYCoordinate = min(nodeCoordinates(:,2));
maxYCoordinate = max(nodeCoordinates(:,2));

if dimension == 3
    minZCoordinate = min(nodeCoordinates(:,3));
    maxZCoordinate = max(nodeCoordinates(:,3));
end
    
diameter_max = max(barDiameters);
diameter_scale = scale_maxDiameter/diameter_max;
diameter_tol = tol_barDiameter * diameter_max;

%% plot initial truss

if dimension == 2
    for bar = 1:n_bars
        if barDiameters(bar) >= diameter_tol
            barCoordinates_x = [nodeCoordinates(potentialBars(bar,1),1) nodeCoordinates(potentialBars(bar,2),1)];
            barCoordinates_y = [nodeCoordinates(potentialBars(bar,1),2) nodeCoordinates(potentialBars(bar,2),2)];

            plot(barCoordinates_x, barCoordinates_y, 'k-', 'lineWidth', barDiameters(bar,1)*diameter_scale);
            hold on
        end
    end

    for node = 1:n_nodes
        if ismember(node, fixedNodes)
            plot(nodeCoordinates(node,1),nodeCoordinates(node,2),'db', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
            hold on
        else 
            plot(nodeCoordinates(node,1),nodeCoordinates(node,2),'og', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
            hold on
        end
    end
    axis('equal')
    axis('off')
    axis([minXCoordinate-0.5 maxXCoordinate+0.5 minYCoordinate-0.5 maxYCoordinate+0.5]);
        
else
    for bar = 1:n_bars
        if barDiameters(bar) >= diameter_tol
            barCoordinates_x = [nodeCoordinates(potentialBars(bar,1),1) nodeCoordinates(potentialBars(bar,2),1)];
            barCoordinates_y = [nodeCoordinates(potentialBars(bar,1),2) nodeCoordinates(potentialBars(bar,2),2)];
            barCoordinates_z = [nodeCoordinates(potentialBars(bar,1),3) nodeCoordinates(potentialBars(bar,2),3)];

            plot3(barCoordinates_x, barCoordinates_y, barCoordinates_z, 'k-', 'lineWidth', barDiameters(bar,1)*diameter_scale);
            hold on
        end
    end

    for node = 1:n_nodes
        if ismember(node, fixedNodes)
            plot(nodeCoordinates(node,1), nodeCoordinates(node,2), nodeCoordinates(node,3), 'db', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
            hold on
        else 
            plot(nodeCoordinates(node,1), nodeCoordinates(node,2), nodeCoordinates(node,3), 'og', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
            hold on
        end
    end
    axis('equal')
    axis('off')
    axis([minXCoordinate-0.5 maxXCoordinate+0.5 minYCoordinate-0.5 maxYCoordinate+0.5  minZCoordinate-0.5 maxZCoordinate+0.5]);
end

hold off