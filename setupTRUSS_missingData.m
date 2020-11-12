function truss = setupTRUSS_missingData(truss)

% Given a truss structure, this function checks the input data for
% completeness and consistency. If possible, missing data is inserted using
% default values.
% Additionally, properties such as the length of all potential bars
% and angles between bars and coordinate axes at their start and end nodes
% are computed and stored in
    % truss.barLengths
    % truss.barAngles

% The truss should be provided as a struct with the following fields:
    % truss.dimension           = 2 or 3
    % truss.nodeCoordinates     given in a (n_nodes x dimension) matrix
    % truss.fixedNodes          given in a list corresponding to the order of nodeCoordinates
    % truss.potentialBars       given in a (n_bars x 2) matrix, where the number of the start/endpoints
    %                           corresponds to the order of nodeCoordinates
    % truss.loadCases           given in a (n_nodes x dimension x n_loadCases) matrix, where each 
    %                           	(n_nodes x dimension) matrix corresponds to one load case and
    %                           	contains the forces applied to each node in the same order as nodeCoordinates
    % truss.compliance_max      a positive number, to be used for all load cases
    % truss.stresses_max        to be used for all load cases, either a positive number to be used for all
    %                           	potential bars or a list, where the order corresponds to potentialBars
    % truss.displacement_max    optional, to be used for all load cases, a positive number
    % truss.barDiameters_max    optional, either a positive number to be used for all potential bars or
    %                               a list, where the order corresponds to potentialBars
    % truss.barDiameters_min    optional, only for realized bars, either a positive number to be used for all potential bars or
    %                               a list, where the order corresponds to potentialBars
    % truss.barDiameters_start  optional, either a positive number to be used for all potential bars or
    %                               a list, where the order corresponds to potentialBars
    

%% parameters

barDiameters_start_default = 1; % default initial value for bar diameters


%% check truss data for completeness and consistency, filling in missing information

% node coordinates need to be provided
if ~isfield(truss, 'nodeCoordinates') || isempty(truss.nodeCoordinates) 
    error('nodeCoordinates are missing')
elseif ~ismember(size(truss.nodeCoordinates, 2), [2, 3])
    error('nodeCoordinates have illegal dimension')
end
n_nodes = size(truss.nodeCoordinates,1);


% dimension of the truss needs to be 2 or 3
if ~isfield(truss, 'dimension') || isempty(truss.dimension) 
    truss.dimension = size(truss.nodeCoordinates, 2);
elseif truss.dimension ~= size(truss.nodeCoordinates, 2)
    disp(['truss dimension is inconsistent with nodeCoordinates, using ' num2str(size(truss.nodeCoordinates, 2))])
    truss.dimension = size(truss.nodeCoordinates, 2);
end
dimension = truss.dimension;


% truss needs to have some fixed nodes
if ~isfield(truss, 'fixedNodes') || isempty(truss.fixedNodes) 
    error('fixedNodes are missing')
elseif any(~ismember(truss.fixedNodes(:), 1:n_nodes))
    error('fixedNodes are not a subset of the truss nodes')
end


% truss needs to have some potential bars
if ~isfield(truss, 'potentialBars') || isempty(truss.potentialBars) 
    disp('potentialBars are missing, computing all potential bars')
    truss = setupTRUSS_computePotentialBars(truss);
elseif any(~ismember(truss.potentialBars(:), 1:n_nodes))
    error('endpoints of some potentialBars are not truss nodes')
end
n_bars = size(truss.potentialBars,1);


% truss needs to have at least one load case
if ~isfield(truss, 'loadCases') || isempty(truss.loadCases) 
    error('load cases are missing')
elseif (size(truss.loadCases,1) ~= n_nodes) || (size(truss.loadCases,2) ~= dimension) 
    error('loadCases have the wrong dimension')
end

% compliance in each load case should be bounded
if ~isfield(truss, 'compliance_max') || isempty(truss.compliance_max) 
    error('compliance_max is missing')
elseif truss.compliance_max <= 0
    error('compliance_max has an illegal value')
end

% displacement of free nodes can be bounded
if ~isfield(truss, 'displacement_max') || isempty(truss.displacement_max) 
    % nodal displacement is not bounded
    truss.displacement_max = inf;
elseif truss.compliance_max <= 0
    error('displacement_max has an illegal value')
end


% stress on realized bars in each load case should be bounded
if ~isfield(truss, 'stresses_max') || isempty(truss.stresses_max)
    error('stresses_max is missing')
elseif ~ismember(length(truss.stresses_max),[1, n_bars])
    error('stresses_max has an illegal dimension')
elseif any(truss.stresses_max <= 0) 
    error('stresses_max has an illegal value')
elseif length(truss.stresses_max) == 1
    disp('only one value for stresses_max given, using this value for all potential bars');
    truss.stresses_max = truss.stresses_max * ones(n_bars,1);
end


% bar diameters can be bounded above
if ~isfield(truss, 'barDiameters_max') || isempty(truss.barDiameters_max)
    truss.barDiameters_max = inf(n_bars,1);
elseif ~ismember(length(truss.barDiameters_max),[1,n_bars])
    error('barDiameters_max has an illegal dimension')
elseif any(truss.barDiameters_max <= 0)
    error('barDiameters_max has an illegal value')
elseif length(truss.barDiameters_max) == 1
    disp('only one value for barDiameters_max given, using this value for all potential bars');
    truss.barDiameters_max = truss.barDiameters_max * ones(n_bars,1);
end

% positive lower bound on diameter of *realized* bars can be given
if ~isfield(truss, 'barDiameters_min') || isempty(truss.barDiameters_min)
    truss.barDiameters_min = zeros(n_bars,1);
elseif ~ismember(length(truss.barDiameters_min),[1,n_bars])
    error('barDiameters_min has an illegal dimension')
elseif length(truss.barDiameters_min) == 1
    disp('only one value for barDiameters_min given, using this value for all potential bars');
    truss.barDiameters_min = truss.barDiameters_min * ones(n_bars,1);
end
if any(truss.barDiameters_min < 0) || any(truss.barDiameters_min > truss.barDiameters_max)
    error('some barDiameters_min are not in [0, barDiameters_max]')
end


% intitial value for diameter of potential bars can be given, otherwise default is used
if ~isfield(truss, 'barDiameters_start') || isempty(truss.barDiameters_start) 
    disp(['no initial value diameter of bars given, using default value ' num2str(barDiameters_start_default)])
    truss.barDiameters_start = barDiameters_start_default * ones(n_bars,1);
elseif ~ismember(length(truss.barDiameters_start),[1,n_bars]) 
    error('barDiameters_start has an illegal dimension')
elseif any(truss.barDiameters_start < 0)
    error('barDiameters_start has an illegal value')
elseif length(truss.barDiameters_start) == 1
    disp('only one value for barDiameters_start given, using this value for all potential bars');
    truss.barDiameters_start = truss.barDiameters_start * ones(n_bars,1);
end


% compute length of all potential bars, i.e. the Euclidean distance between start and end nodes
truss.barLengths = computeBarLengths(truss);

% for each bar and each corresponding free startNode or endNode, 
% determine -cos of the intermediate angles between the bar and the x/y/z-axis
% and save them in a (n_nodes x dimension x n_bars) sized tensor
truss.barAngles = computeBarAngles(truss);

end



%% auxiliary functions computing properties of the truss

function barLengths = computeBarLengths(truss)
    % for each bar compute the corresponding length, i.e. the euclidean distance between startNode and endNode
    nodeCoordinates = truss.nodeCoordinates;
    potentialBars = truss.potentialBars;
    startNodes = potentialBars(:,1);
    endNodes = potentialBars(:,2);
    
    barLengths = vecnorm(nodeCoordinates(startNodes,:) - nodeCoordinates(endNodes,:),2,2);
end 


function barAngles = computeBarAngles(truss) 
    % For each bar and each corresponding free startNode or endNode, 
    % determine -cos of the intermediate angles between the bar and the x-axis, y-axis, (and z-axis).
    dimension = truss.dimension;
    nodeCoordinates = truss.nodeCoordinates;
    n_nodes = size(nodeCoordinates,1);
    fixedNodes = truss.fixedNodes;
    potentialBars = truss.potentialBars;
    barLengths = truss.barLengths;
    n_bars = size(potentialBars,1);
    startNodes = potentialBars(:,1);
    endNodes = potentialBars(:,2);
    
    barAngles = zeros(n_nodes, dimension, n_bars);
    for bar = 1:n_bars
        startNode = startNodes(bar);
        endNode = endNodes(bar);
        
        if ~ismember(startNode, fixedNodes)
            % startNode of bar is not fixed
            barAngles(startNode, :, bar) = - (nodeCoordinates(endNode,:) - nodeCoordinates(startNode,:))/barLengths(bar);
        end
        
        if ~ismember(endNode, fixedNodes)
            % endNode of bar is not fixed
            barAngles(endNode, :, bar) = - (nodeCoordinates(startNode,:) - nodeCoordinates(endNode,:))/barLengths(bar);
        end
    end
end
