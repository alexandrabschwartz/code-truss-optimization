function truss = setupTRUSS_computePotentialBars(truss)

% Given a truss structure with free and fixed nodes, this function computes
% all potential bars, where
    % - at least one end node is free
    % - no other node lies on this bar.
% The function needs a truss struct with (at least) the fields
    % truss.dimension
    % truss.nodeCoordinates
    % truss.fixedNodes
% and returns a truss with the new or modified field
    % truss.potentialBars
    
%% check for necessary input data    

if ~isfield(truss, 'dimension') || isempty(truss.nodeCoordinates)
    error('truss dimension is missing')
end
dimension = truss.dimension;

if ~isfield(truss, 'nodeCoordinates') || isempty(truss.nodeCoordinates)
    error('nodeCoordinates are missing')
end
nodeCoordinates = truss.nodeCoordinates;
n_nodes = size(nodeCoordinates, 1);

if ~isfield(truss, 'fixedNodes') || isempty(truss.fixedNodes)
    error('fixedNodes are missing')
end
fixedNodes = truss.fixedNodes;
% n_fixedNodes = length(fixedNodes);



%% determine potential bars
n_bars = 0;
potentialBars = zeros(n_nodes^2,2);

for startNode = 1:n_nodes-1
    for endNode = startNode+1:n_nodes

        % check if both nodes are fixed
        if ismember(startNode,fixedNodes) && ismember(endNode,fixedNodes)
            % both nodes are fixed, hence no bar is necessary

        else % check, if there would be other nodes on a bar between these two nodes
            if all(nodeCoordinates(startNode,:) == nodeCoordinates(endNode,:))
                error ('There are two nodes in the same spot')
            else
                % the nodes have at least one different coordinate
                dim_diff = find(nodeCoordinates(startNode,:) ~= nodeCoordinates(endNode,:), 1);
                % find all nodes where the same coordinate takes values between startNode and endNode
                criticalNodes = find((nodeCoordinates(:,dim_diff) > min(nodeCoordinates(startNode,dim_diff),nodeCoordinates(endNode,dim_diff)))...
                        & (nodeCoordinates(:,dim_diff) < max(nodeCoordinates(startNode,dim_diff),nodeCoordinates(endNode,dim_diff))));
                    n_criticalNodes = length(criticalNodes);
                % other dimensions
                dim_others = setdiff(1:dimension, dim_diff);
                % compare the other coordiantes of the critical nodes with the potential bar
                nodeComparison = (nodeCoordinates(criticalNodes,dim_others) == repmat(nodeCoordinates(startNode,dim_others),n_criticalNodes,1) + ...
                        repmat(nodeCoordinates(criticalNodes,dim_diff) - nodeCoordinates(startNode,dim_diff),1,length(dim_others)) ./...
                        repmat(nodeCoordinates(endNode,dim_diff) - nodeCoordinates(startNode,dim_diff),n_criticalNodes,length(dim_others)) .* ...
                        repmat(nodeCoordinates(endNode,dim_others) - nodeCoordinates(startNode,dim_others),n_criticalNodes,1));
                if any(nodeComparison(:,1) .* nodeComparison(:,end) == 1)
                      % there would be another node on this bar, so no bar is possible
                else % add the bar to the list of potential bars
                    n_bars = n_bars + 1;
                    potentialBars(n_bars,:) = [startNode endNode];
                end 
            end
        end
    end
end

truss.potentialBars = potentialBars(1:n_bars,:);


            
