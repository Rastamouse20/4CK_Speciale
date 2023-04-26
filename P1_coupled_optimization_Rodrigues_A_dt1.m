%% Clear everything from memory and start from scratch
clear all; close all; clc;
addpath('topopt\')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                  PREPROCESSING                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the problem parameters
Emax        = 210e9;        % [Pa = N/m^2]
Emin        = 1e-9*Emax;    % [Pa = N/m^2]
nu          = 0.3;          % [-]       Poisson
CTE         = 15e-6;        % [K^-1]    Thermal exspantion 
kappamax    = 80;           % [w/(m*k)] Thermal conductivity
kappamin    = 1e-9*kappamax;% [w/(m*k)] Thermal conductivity
heat_rate   = 1;            % [W]       Heat rate in the geometry
Ptotal      = 1e3;          % [N]       Mechanical force
%% Domain size and discretisation
Ly  = 0.48;         % [m]
Lx  = 0.72;         % [m]
t   = 0.01;         % [m]
elemy = 30; nodey = elemy+1; %%% 30/240 
elemx = 60; nodex = elemx+1; %%% 60/120
voidON = 0;             % Passiveh = 1 active passive circle domain otherwise = 0
solidON = 1;            % Passives = 1 active passive solid domain otherwise = 0
%% Volume constraint
vtype = 2;         % Maximum volume type: 1 = specific volume; 2 = volume fraction
Vmax = 6.0e-04;    % [m^3]
volfrac = 0.4;     % [1/100 %]
%% Optimisation parameters
maxiter = 100; tol = 1e-2; mvlim = 0.2;
filter = 0; rmin = 1.2;     % Filter on/off[1/0], filter radius [*elementSize]
p = 3;                      % [-] Penalty factor for SIMP
%% Compute the nodal coordinates and the element connectivity
[a,b,nodalCoordinates,elementConnectivity] = quadgen(Lx,Ly,elemx,elemy);
numberOfNodes = size(nodalCoordinates,1);       % extracts the number of nodes based
                                                % on the number of rows in connectivity table
numberOfElements = size(elementConnectivity,1); % extracts the number of elements based
                                                % on the number of rows in connectivity table
numberOfDOFs    = 2*numberOfNodes;  % MECHANICAL each node has 2 DOFs: x and y displacement
numberOfDOFs_t  = 1*numberOfNodes;  % THERMAL each node has 1 DOFs: [K ]
%% Passive domain
centerElement = zeros(numberOfElements,3);
v=1; s=1; solid = 0; void = 0;
for ii = 1:numberOfElements
    centerElement(ii,:) = [ii mean(nodalCoordinates(elementConnectivity(ii,2:end),2:end))]; % [Ele_number, Ele_x_pos, Ele_y_pos]  
    if (voidON), end
    if (solidON), if centerElement(ii,2) < 2*Lx/elemx || centerElement(ii,2) > Lx - 2*Lx/elemx
            solid(s) = centerElement(ii);   s = s+1; end, end
end
%% Specify which DOFs are fixed
% Mechanical
fixedDOFs = [1:2*nodey numberOfDOFs-2*nodey:numberOfDOFs];  % Fixed DOFs
allDOFs = 1:numberOfDOFs;                                   % all DOFs
freeDOFs = setdiff(allDOFs,fixedDOFs);                      % extracts only free DOFs by removing the fixed DOFs
% Thermal
fixedDOFs_t = [1]; %[(nx-1)*ny+1:1:nx*ny];                   % Fixed Thermal DOFs 
allDOFs_t = 1:numberOfDOFs_t;                               % all DOFs
freeDOFs_t = setdiff(allDOFs_t,fixedDOFs_t);                % extracts only free DOFs by removing the fixed DOFs
%% Define the DOFs with loads
% Mechanical
forceDOFs = 2*(elemx/2*nodey)+2;        % define the DOFs that have loads
forceValues = -Ptotal;                  % define the force values in the same order as the DOFs above
    % Define the global force vector
f = sparse(numberOfDOFs,1);             % the force vector is set to 0 initially
f(forceDOFs) = forceValues;             % inserts all forces at once into the right global DOFs
% Thermal
forceDOFs_t = [1:1:round(nodey/2)];                                         % define the DOFs that have loads
forceValues_t = [1 2*ones(1,round(nodey/2)-2) 1]*heat_rate/(elemy/2*2);     % define the force values in the same order as the DOFs above
    % Define the global forve vector
f_t = sparse(numberOfDOFs_t,1);                                             % the force vector is set to 0 initially
f_t(forceDOFs_t) = forceValues_t;                                           % inserts all forces at once into the right global DOFs
%% Precompute reference stiffness matrix and DOFs for assembly
% Mechanical
k0 = stiffnessMatrix(a,b,t,1,nu);
nodes = elementConnectivity(:,2:5)';
dofs = reshape([2*nodes(:)-1 2*nodes(:)]',8,numberOfElements);
iK = reshape(kron(dofs,ones(8,1)),64*numberOfElements,1);
jK = reshape(kron(dofs,ones(1,8)),64*numberOfElements,1);
% Thermal
k0_t = thermalstiffnessMatrix(a,b,t,1);
nodes_t = elementConnectivity(:,2:5)';
dofs_t = nodes_t;
iK_t = reshape(kron(dofs_t,ones(4,1)),16*numberOfElements,1);
jK_t = reshape(kron(dofs_t,ones(1,4)),16*numberOfElements,1);
%% Prepare filter
if (filter); filt = filterPrep(elemx,elemy,nodes,rmin); end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                 INITIAL DESIGN                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define design variables and allowable volume
if (vtype == 1)
    if (Vmax < 0); error('Vmax is lower than zero! \n');
    elseif (Vmax > Lx*Ly); error('Vmax is higher than allowed by the domain! \n'); end
    fprintf('Vmax is equivalent to volfrac = %3.2f percent. \n \n',100*(Vmax/(Lx*Ly*t)));
elseif (vtype == 2); Vmax = volfrac*Lx*Ly*t;
    fprintf('Volume fraction gives Vmax = %3.2e. \n \n',Vmax); end
rho = volfrac*ones(numberOfElements,1);

if (voidON); rho(void)=0; end; if (solidON); rho(solid)=1; end              % Filter
if (filter); rhof = filterApply(rho,filt); else; rhof = rho; end 
if (voidON); rhof(void)=0; end; if (solidON); rhof(solid)=1; end            % Filter

elementStiffness = Emin+(Emax-Emin)*rhof.^p;                                % define initial stiffness
elementConductivity = kappamin+(kappamax-kappamin)*rhof.^p;                 % define initial conductivity
%% Initialise MMA parameters - Rasmus
rhoold1 = rho; rhoold2 = rho; low = 0; upp = 0;
numberOfConstraints = 1;
a0 = 1; ai = 0*ones(numberOfConstraints,1);
c = 1000*ones(numberOfConstraints,1);
d = 1*ones(numberOfConstraints,1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                 OPTIMISATION                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimisation loop
change = Inf;
for iter = 1:maxiter

    %% Assemble the global stiffness matrix
    % Mechanical
    sK = reshape(k0(:)*elementStiffness(:)',64*numberOfElements,1);
    globalK = sparse(iK,jK,sK); globalK = (globalK+globalK')/2;
    % Termal
    sK_t = reshape(k0_t(:)*elementConductivity(:)',16*numberOfElements,1);
    globalK_t = sparse(iK_t,jK_t,sK_t); globalK_t = (globalK_t+globalK_t')/2;
    %% Solve the Heat Transfer system for the free nodes
%     freeK_t = globalK_t(freeDOFs_t,freeDOFs_t);     % extracts only rows and columns for the free nodes,
%                                                     % i.e. removes the rows and columns for fixed  nodes
%     u_t = zeros(numberOfDOFs_t,1);                  % the vector of unknowns is set to zero initially
%     u_t(freeDOFs_t) = freeK_t\f_t(freeDOFs_t);      % solves the system
    u_t = ones(numberOfDOFs_t,1)*1;                  % the vector of unknowns is set to zero initially %% DELETE!!!
    %% Thermal load vector
    alpha = [CTE; CTE; 0];                  % Thermal expansion
    A_e0 = CouplingMatrix(CTE,a,b,t,nu);    % retrieve coupling matrix
    dfdt = sparse(numberOfDOFs,numberOfNodes);  % EXPERIMENTAL
    % Force calculation
    for e = 1:numberOfElements
        nodes_e = elementConnectivity(e,2:5);             % Extract nodes_e for element e
        dofs_e = reshape([2*nodes_e-1; 2*nodes_e],8,1);   % Compute DOF numbers
        A_e = elementStiffness(e) * A_e0;                 % Ae, coupling matrix
        f(dofs_e) = f(dofs_e) + A_e * u_t(nodes_e);       % Force
        dfdt(dofs_e,nodes_e) = dfdt(dofs_e,nodes_e) + A_e;  % EXPERIMENTAL Jacobian
    end
    %% Solve the system for the free nodes
    freeK = globalK(freeDOFs,freeDOFs);             % extracts only rows and columns for the free nodes,
                                                    % i.e. removes the rows and columns for fixed  nodes
    u = zeros(numberOfDOFs,1);                      % the vector of unknowns is set to zero initially
    u(freeDOFs) = freeK\f(freeDOFs);                % solves the system
    %% Solve for adjoint variable(1)
    lmda1 = zeros(numberOfNodes,1);
    rhs = 2*dfdt'*u;
    lmda1(freeDOFs_t) = globalK_t(freeDOFs_t,freeDOFs_t)\rhs(freeDOFs_t);
    %% Compute the compliance and sensitivity %% CHANGE THIS ONE
    temp = reshape(sum((u(dofs')*k0).*u(dofs'),2),numberOfElements,1);
    C = sum(temp.*elementStiffness);
    dC = zeros(numberOfElements,1);
    for e = 1:numberOfElements
        nodes_e = elementConnectivity(e,2:5);             % Extract nodes_e for element e
        dofs_e = reshape([2*nodes_e-1; 2*nodes_e],8,1);   % Compute DOF numbers
        
        dC(e,1) = lmda1(nodes_e)'*(0 - (kappamax-kappamin)*p*rhof(e)^(p-1)*k0_t*u_t(nodes_e)) ...
                    - u(dofs_e)'*(Emax-Emin)*p*rhof(e)^(p-1)*k0*u(dofs_e) ...
                    + 2*((Emax-Emin)*p*rhof(e)^(p-1)*A_e0*u_t(nodes_e))'*u(dofs_e);
    end
    % dC = -temp.*(p*(Emax-Emin)*rhof(:).^(p-1)); % CHANGE THIS ONE ONLY
    
    if (voidON); dC(void)=0; end; if (solidON); dC(solid)=0; end            % Filter!!!
    if (filter); dC = filterApply(dC,filt); end
    if (voidON); dC(void)=0; end; if (solidON); dC(solid)=0; end            % Filter!!!
    %% Compute the volume and sensitivity
    Ve = 4*a*b*t; V = sum(Ve*rhof);
    dV = Ve*ones(numberOfElements,1);
    if (voidON); dV(void)=0; end; if (solidON); dV(solid)=0; end            % Filter!!!
    if (filter); dV = filterApply(dV,filt); end
    if (voidON); dV(void)=0; end; if (solidON); dV(solid)=0; end            % Filter!!!
    %% Plot the current design
    plot2DElements(1,nodalCoordinates,elementConnectivity,u,[],rhof,1); 
    % Convergence check
    fprintf('Iter. %3i - C = %4.3e, V = %4.3e, change = %3.2e \n',iter,C,V,change); 
    if (change < tol || iter == maxiter); break; end
    %% Scaling
    if (iter == 1); C0 = 20; end
    f0 = 1+C/C0; df0 = dC/C0;
    f1 = V/Vmax - 1; df1 = dV/Vmax;
    %% Call MMA 
    rholow = max(0,rho-mvlim); rhoupp = min(1,rho+mvlim);
    [rhonew,~,~,~,~,~,~,~,~,low,upp] =  mmasub(numberOfConstraints,numberOfElements,...
                                        iter,rho,rholow,rhoupp,rhoold1,rhoold2,...
                                        f0,df0,f1,df1',low,upp,a0,ai,c,d);
    change = max(abs(rho-rhonew));
    rhoold2 = rhoold1; rhoold1 = rho; rho = rhonew;

    if (voidON); rho(void)=0; end; if (solidON); rho(solid)=1; end          % Filter!!!
    if (filter); rhof = filterApply(rhonew,filt); else rhof = rhonew; end
    if (voidON); rhof(void)=0; end; if (solidON); rhof(solid)=1; end        % Filter!!!
    
    elementStiffness = Emin + (Emax-Emin)*rhof.^p;
    elementConductivity = kappamin + (kappamax-kappamin)*rhof.^p;

    p = min([p+0.25 3]);                                    % Increasin Penalty Factor
end
%% Measure of non-discreteness
Mnd = full(4*sum(rhof(:).*(1-rhof(:)))/numberOfElements*100)
figure(2); histogram(rho,10); hold on; histogram(rhof,10); legend('rho','rhof');