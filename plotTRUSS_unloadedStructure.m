function plotTRUSS_unloadedStructure(truss)

% This function plots the optimized, but unloaded truss.
% The function needs a truss struct with the fields
    % truss.dimension
    % truss.nodeCoordinates
    % truss.fixedNodes
    % truss.potentialBars
    % truss.barDiameters_opt 
    
    
%% parameters

tol_barDiameter = 10^-4;
scale_maxDiameter = 10;


%% read truss data

if ~isfield(truss, 'dimension') || isempty(truss.dimension)
    error('Plot unloaded truss: Truss dimension (= 2 or 3) is missing')
else 
    dimension = truss.dimension;
end

if ~isfield(truss, 'nodeCoordinates') || isempty(truss.nodeCoordinates)
    error('Plot unloaded truss: Node coordinates missing')
else 
    nodeCoordinates = truss.nodeCoordinates;
    n_nodes = size(nodeCoordinates, 1);
end

if ~isfield(truss, 'fixedNodes') || isempty(truss.fixedNodes)
    error('Plot unloaded truss: Fixed nodes are missing')
else
    fixedNodes = truss.fixedNodes;
end

if ~isfield(truss, 'potentialBars') || isempty(truss.potentialBars)
    error('Plot unloaded truss: Potential bars are missing')
else
    potentialBars = truss.potentialBars;
    n_bars = size(potentialBars,1);
end

        
if ~isfield(truss,'barDiameters_opt') || isempty(truss.barDiameters_opt)
    error('Optimized diameters for the bars are missing')
else
    barDiameters = truss.barDiameters_opt;
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


%% plot unloaded truss

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