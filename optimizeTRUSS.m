function [truss, information] = optimizeTRUSS(truss, options)

% Given a truss structure with nodes, potential bars and fixed nodes, this 
% function tries to minimize the total volume of the truss, while satisfying
% certain physical constraints such as bounds on deformation, stress or the
% diameter of bars.

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
    
% If you want to specify the MPVC algorithm or NLP solver, decide on the use
% of gradients or slack variables for the vanishing constraints, additionally
% provide an options struct with the fields
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.solver     = 'fmincon' or 'snopt'
    % options.algorithm  = 'direct' or 'relaxation' or 'relaxation_posLB'
    % options.relaxation = 'scholtes' or 'steffensen' or 'schwartz' or 'kadrani'
    % options.slacks     = true or false
    
% The function returns
    % truss                             the truss struct including optimal diameters and displacements
    % information.message               final exit message of the NLP solver
    % information.iterations            number of NLPs solved    
    % information.maxVio_diameter       maximum violation of constraints on the diameter of bars 
    %                                      (>= 0, >= barDiameters_min, <= barDiameters_max)
    % information.maxVio_displacement   maximum violation of displacement constraints
    % information.maxVio_equilibrium    maximum violation of the force equilibrium
    % information.maxVio_compliance     maximum violation of compliance constraints
    % information.maxVio_stress         maximum violation of stress constraints
    
    
%% parameters

E = 1;                   % Young's modulus
bar_density = 1;         % density of each bar
bar_tolerance = 10^-4;   % tolerance to decide if a bar's diameter is positive in the optimum

%% check truss data for completeness and consistency, insert missing information

truss = setupTRUSS_missingData(truss);


%%  extract information from the truss

% dimension of the truss (= 2 or 3)
dimension = truss.dimension;

% number of nodes and their coordinates in x/y/z-direction
nodeCoordinates = truss.nodeCoordinates;
n_nodes = size(nodeCoordinates,1);

% identify fixed and free nodes
fixedNodes = truss.fixedNodes;
n_fixedNodes = length(fixedNodes);
freeNodes = setdiff(1:n_nodes, fixedNodes);
n_freeNodes = n_nodes - n_fixedNodes;

% number of potential bars and their start and end nodes
potentialBars = truss.potentialBars;
n_bars = size(potentialBars,1);

% number of load cases and applied forces
loadCases = truss.loadCases;
n_loadCases = size(loadCases,3);

% upper bound on compliance, to be used in all load cases
compliance_max = truss.compliance_max;

% upper bound on the displacement of the free nodes, to be used in all load cases
displacement_max = truss.displacement_max;

% upper bound on stress for each bar, to be used in all load cases
stresses_max = truss.stresses_max;

% upper bound on diameter for each bar
barDiameters_max = truss.barDiameters_max;

% nonslender bars, i.e. bars with positive lower bound on the diameter
barDiameters_min = truss.barDiameters_min;
nonslenderBars = find(barDiameters_min > 0);
n_nonslenderBars = length(nonslenderBars);

% initial values for the diameters of potential bars
barDiameters_start = truss.barDiameters_start;

% length of all potential bars
barLengths = truss.barLengths;

% -cos of the angle between each bar and the x/y/z axes at the corresponding free end- and startnodes
barAngles = truss.barAngles;


%% define the MPVC

% variable X = (barDiameters, vectorized displacements)
% barDiameters are the diameters of all potential bars
% displacements are the displacements of all free nodes in x/y/z-direction for all loadcases

% displacements are a tensor of dimension (n_freeNodes x dimension x n_loadcases)
% displacements can be vectorized by reshape(displacements, n_freeNodes * dimension * n_loadCases, 1)
% vectorized displacements can be tensorized by reshape(displacements, n_freeNodes, dimension, n_loadCases)
% vectorized displacements are of the form [displacements in loadCase 1;
%                                           displacements in loadCase 2;
%                                           ...]
% displacements in loadCase i are of the form [displacements in x-direction;
%                                              displacements in y-direction;
%                                              displacements in z-direction]
% within one direction (x/y/z), the free nodes are ordered as usual


% number of variables in the MPVC
n_displacements = n_freeNodes * dimension * n_loadCases; 
n_x = n_bars + n_displacements; % total number of variables
problem.dimension = n_x;

% objective
problem.objective = @objective;

% box constraints
% diameters have to be nonegative and cannot exceed barDiameter_max
% positive lower bounds barDiameter_min on *realized* bars are handled in 
% the nonlinear constraints 
diameters_l = zeros(n_bars, 1);
diameters_u = barDiameters_max;

% displacement of free nodes can be bounded
displacements_l = -displacement_max * ones(n_displacements, 1);
displacements_u =  displacement_max * ones(n_displacements, 1);

problem.xl = [diameters_l; displacements_l];
problem.xu = [diameters_u; displacements_u];

% linear constraints
% compliance cannot exceed compliance_max in all load cases
forces = zeros(n_loadCases, n_freeNodes*dimension, n_loadCases);
for loadcase = 1:n_loadCases
    forces(loadcase, :, loadcase) = reshape(loadCases(freeNodes, :, loadcase), n_freeNodes*dimension, 1)';
end
A = [zeros(n_loadCases,n_bars) reshape(forces, n_loadCases, n_displacements)];

problem.A = A;
problem.bl = -inf(n_loadCases,1);
problem.bu = compliance_max*ones(n_loadCases,1);

% nonlinear constraints
% force equilibrium has to be satisfied
problem.nlcons = @nlcons;
problem.cl = reshape(loadCases(freeNodes,:,:),n_displacements,1); 
problem.cu = problem.cl;

% vanishing constraints
% stress on realized bars cannot exceed stress_max
% optional: positive lower bounds barDiameter_min have to be satisfied in realized bars
problem.vancons = @vancons;

% initial value
stiffnessMatrix_start = computeStiffnessMatrix(truss, barDiameters_start);
nodeDisplacements_start = zeros(n_freeNodes*dimension,n_loadCases);
for loadcase = 1:n_loadCases
    force = reshape(loadCases(freeNodes,:,loadcase),n_freeNodes*dimension,1);
    nodeDisplacements_start(:,loadcase) = stiffnessMatrix_start \ force;
end
problem.x_start = [barDiameters_start; reshape(nodeDisplacements_start,n_displacements,1)];



%% solve the MPVC 

[X_opt, f_opt, MPVCinformation] = solveMPVC(problem, options);
truss.volume_opt = f_opt;
truss.barDiameters_opt = X_opt(1:n_bars,1);
truss.nodeDisplacements_opt = reshape(X_opt(n_bars+1:end), n_freeNodes, dimension, n_loadCases);
% note that only free nodes have displacements (per load case)

%% compute return values

information.iterations = MPVCinformation.iterations;
information.message = MPVCinformation.message;

information.maxVio_diameter = max([max(truss.barDiameters_opt - truss.barDiameters_max, 0);...
                                   max(0 - truss.barDiameters_opt, 0);...
                                   max((truss.barDiameters_min - truss.barDiameters_opt).*(truss.barDiameters_opt >= bar_tolerance), 0)]);
                               
information.maxVio_displacement = max([abs(truss.nodeDisplacements_opt(:))-truss.displacement_max; 0]);                              
                               
information.maxVio_compliance = 0;
information.maxVio_equilibrium = 0;
stiffnessMatrix_opt = computeStiffnessMatrix(truss, truss.barDiameters_opt);
for loadcase = 1:n_loadCases
    force = reshape(loadCases(freeNodes,:,loadcase),n_freeNodes*dimension,1);
    displacement = reshape(truss.nodeDisplacements_opt(:,:,loadcase), n_freeNodes*dimension,1);
    information.maxVio_equilibrium = max(norm(stiffnessMatrix_opt * displacement - force,inf),information.maxVio_equilibrium);
    information.maxVio_compliance = max(force' * displacement  - truss.compliance_max, information.maxVio_compliance);
end

information.maxVio_stress = 0;
for loadcase = 1:n_loadCases
    displacement = reshape(truss.nodeDisplacements_opt(:,:,loadcase), n_freeNodes*dimension,1);
    for Bar = 1:n_bars
        angles_bar = reshape(barAngles(freeNodes,:,Bar),n_freeNodes*dimension,1);
        stress_bar = E/barLengths(Bar) * angles_bar' * displacement;  
        information.maxVio_stress = max((abs(stress_bar) - truss.stresses_max(Bar))*(truss.barDiameters_opt(Bar) >= bar_tolerance) , information.maxVio_stress);
    end
end


%% objective function and nonlinear/vanishing constraints for the MPVC

function [f,Df] = objective(X)
    % minimize the total volume of the truss
    barDiameters = X(1:n_bars);
    
    f = bar_density * barLengths' * barDiameters;
    
    if nargout > 1
        Df = [bar_density*barLengths' zeros(1,n_displacements)];
    end
end


function [c, Dc] = nlcons(X)
    % ensure force equilibrium
    barDiameters = X(1:n_bars);
    nodeDisplacements = reshape(X(n_bars+1:end), n_freeNodes, dimension, n_loadCases);
    
    stiffnessMatrix = computeStiffnessMatrix(truss, barDiameters);
        
    c = zeros(n_freeNodes*dimension, n_loadCases);
    for scenario = 1:n_loadCases
        displacement = reshape(nodeDisplacements(:, :, scenario), n_freeNodes * dimension, 1);
        c(:,scenario) = stiffnessMatrix * displacement;
    end
    c = reshape(c,n_displacements,1);

    if nargout > 1
        Dc = zeros(n_displacements, n_x);
        for scenario = 1:n_loadCases
            displacement = reshape(nodeDisplacements(:, :, scenario), n_freeNodes * dimension, 1);
            Dc((scenario-1)*dimension*n_freeNodes+1:scenario*dimension*n_freeNodes, n_bars + ((scenario-1)*dimension*n_freeNodes+1:scenario*dimension*n_freeNodes)) = stiffnessMatrix;
            
            for bar = 1:n_bars
                angles_bar = reshape(barAngles(freeNodes,:,bar),n_freeNodes*dimension,1);
                Dc((scenario-1)*dimension*n_freeNodes+1:scenario*dimension*n_freeNodes, bar) = E/barLengths(bar) * angles_bar * (angles_bar') * displacement;
            end
        end
    end
end


function [G, H, DG, DH] = vancons(X)
    % bound the stress on realized bars
    % optional: ensure positive lower bounds on the diameter of realized bars
     
    % G = [stress(in loadCase 1)^2 - stressMax^2                    H = [barDiameters;              % for loadCase 1
    %      stress(in loadCase 2)^2 - stressMax^2                         barDiameters;              % for loadCase 2
    %      ...                                                           ...                        % repeat for each loadCase 
    %      barDiameters_min(nonslender) - barDiameters(nonslender)]      barDiameters(nonslender)]  % only once, only nonslenderBars
    
    barDiameters = X(1:n_bars);
    nodeDisplacements = reshape(X(n_bars+1:end), n_freeNodes, dimension, n_loadCases);
    
    if nargout == 2
        H = [repmat(barDiameters, n_loadCases, 1); barDiameters(nonslenderBars)];

        G = zeros(n_bars*n_loadCases+n_nonslenderBars, 1);
        for scenario = 1:n_loadCases
            displacement = reshape(nodeDisplacements(:, :, scenario), n_freeNodes * dimension, 1);
            for bar = 1:n_bars
                angles_bar = reshape(barAngles(freeNodes,:,bar),n_freeNodes*dimension,1);
                stress_bar = E/barLengths(bar) * angles_bar' * displacement;
                G((scenario-1)*n_bars + bar,1) = stress_bar^2 - stresses_max(bar)^2;
            end
        end
        G(n_bars*n_loadCases+1:end,1) = barDiameters_min(nonslenderBars) - barDiameters(nonslenderBars);
    
    else 
        I_nonslender = eye(n_bars);
        I_nonslender = I_nonslender(nonslenderBars, :);
        
        H = [repmat(barDiameters, n_loadCases, 1); barDiameters(nonslenderBars)];
        DH = [repmat([eye(n_bars, n_bars) zeros(n_bars, n_displacements)], n_loadCases, 1);...
              I_nonslender                zeros(n_nonslenderBars, n_displacements)];
        
        G = zeros(n_bars*n_loadCases+n_nonslenderBars, 1);
        DG = zeros(n_bars * n_loadCases, n_x);
        for scenario = 1:n_loadCases
            displacement = reshape(nodeDisplacements(:, :, scenario), n_freeNodes * dimension, 1);
            for bar = 1:n_bars
                angles_bar = reshape(barAngles(freeNodes,:,bar),n_freeNodes*dimension,1);
                stress_bar = E/barLengths(bar) * angles_bar' * displacement;
                Dstress_bar = E/barLengths(bar) * angles_bar';
                G((scenario-1)*n_bars + bar,1) = stress_bar^2 - stresses_max(bar)^2;
                DG((scenario-1)*n_bars + bar, n_bars + ((scenario-1)*n_freeNodes*dimension+1:scenario*n_freeNodes*dimension)) = 2 * stress_bar * Dstress_bar ;          
            end
        end
        G(n_bars*n_loadCases+1:end,1) = barDiameters_min(nonslenderBars) - barDiameters(nonslenderBars);
        DG = [DG; -I_nonslender zeros(n_nonslenderBars, n_displacements)];
    end
end



%% auxiliary functions computing properties of the truss

function stiffnessMatrix = computeStiffnessMatrix(truss, barDiameters)
    % computes the stiffness matrix for a given truss configuration
    % this matrix is used to compute the nodal displacements in the force equilibrium
    dimension = truss.dimension;
    n_bars = size(truss.potentialBars,1);
    barAngles = truss.barAngles;
    barLengths = truss.barLengths;
    freeNodes = setdiff(1:n_nodes, truss.fixedNodes);
    
    stiffnessMatrix = zeros(n_freeNodes * dimension, n_freeNodes * dimension);
    for bar = 1:n_bars
        angles_bar = reshape(barAngles(freeNodes,:,bar), n_freeNodes * dimension, 1);
        stiffnessMatrix = stiffnessMatrix + E/barLengths(bar) * barDiameters(bar) * angles_bar * (angles_bar');
    end
end


end