%% Clear everything from memory and start from scratch
clear all; close all; clc;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                  PREPROCESSING                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the problem parameters
Emax = 70e9;     % [Pa = N/m^2]
Emin = 1e-6*Emax;% [Pa = N/m^2]
p = 3.0;         % [-]
nu = 0.3;        % [-]
Ptotal = 1000;   % [N]
%% Domain size and discretisation
Ly = 0.2;        % [m]
Lx = 2*Ly;       % [m]
t = 25e-3;       % [m]
elemy = 40; nodey = elemy+1;
elemx = elemy*Lx/Ly; nodex = elemx+1; 
%% Volume constraint
vtype = 2;         % Maximum volume type: 1 = specific volume; 2 = volume fraction
Vmax = 6.0e-04;   % [m^3]
volfrac = 0.3;     % [1/100 %]
%% Optimisation parameters
maxiter = 100; tol = 1e-2; mvlim = 0.2;
filter = 1; rmin = 1.4;
%% Compute the nodal coordinates and the element connectivity
[a,b,nodalCoordinates,elementConnectivity] = quadgen(Lx,Ly,elemx,elemy);
numberOfNodes = size(nodalCoordinates,1); % extracts the number of nodes based
                                          % on the number of rows in connectivity table
numberOfDOFs = 2*numberOfNodes; % each node has 2 DOFs: x and y displacement
numberOfElements = size(elementConnectivity,1); % extracts the number of elements based
                                                % on the number of rows in connectivity table                           
%% Specify which DOFs are fixed
fixedDOFs = 1:2*nodey; % Fixed DOFs
allDOFs = 1:numberOfDOFs; % all DOFs
freeDOFs = setdiff(allDOFs,fixedDOFs); % extracts only free DOFs by removing the fixed DOFs
%% Define the DOFs with loads
forceDOFs = 2*((nodex-1)*nodey+1); % define the DOFs that have loads
forceValues = -Ptotal; % define the force values in the same order as the DOFs above
% Define the global force vector
f = sparse(numberOfDOFs,1); % the force vector is set to 0 initially
f(forceDOFs) = forceValues; % inserts all forces at once into the right global DOFs
%% Precompute reference stiffness matrix and DOFs for assembly
k0 = stiffnessMatrix(a,b,t,1,nu);
nodes = elementConnectivity(:,2:5)';
dofs = reshape([2*nodes(:)-1 2*nodes(:)]',8,numberOfElements);
iK = reshape(kron(dofs,ones(8,1)),64*numberOfElements,1);
jK = reshape(kron(dofs,ones(1,8)),64*numberOfElements,1);
%% Prepare filter
if (filter); filt = filterPrep(elemx,elemy,nodes,rmin); end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                 INITIAL DESIGN                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define design variables and allowable volume
if (vtype == 1)
    if (Vmax < 0)
        error('Vmax is lower than zero! \n');
    elseif (Vmax > Lx*Ly)
        error('Vmax is higher than allowed by the domain! \n');
    end
    fprintf('Vmax is equivalent to volfrac = %3.2f percent. \n \n',100*(Vmax/(Lx*Ly*t)));
elseif (vtype == 2)
    Vmax = volfrac*Lx*Ly*t;
    fprintf('Volume fraction gives Vmax = %3.2e. \n \n',Vmax); 
end
rho = volfrac*ones(numberOfElements,1);
if (filter); rhof = filterApply(rho,filt); else; rhof = rho; end
elementStiffness = Emin+(Emax-Emin)*rhof.^p;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                 OPTIMISATION                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimisation loop
change = Inf;
for iter = 1:maxiter
    %% Assemble the global stiffness matrix
    sK = reshape(k0(:)*elementStiffness(:)',64*numberOfElements,1);
    globalK = sparse(iK,jK,sK); globalK = (globalK+globalK')/2;   
    %% Solve the system for the free nodes
    freeK = globalK(freeDOFs,freeDOFs); % extracts only rows and columns for the free nodes,
                                          % i.e. removes the rows and columns for fixed  nodes
    u = zeros(numberOfDOFs,1); % the vector of unknowns is set to zero initially
    u(freeDOFs) = freeK\f(freeDOFs); % solves the system
    %% Compute the compliance and sensitivity
    temp = reshape(sum((u(dofs')*k0).*u(dofs'),2),numberOfElements,1);
    C = sum(temp.*elementStiffness);
    dC = -temp.*(p*(Emax-Emin)*rhof(:).^(p-1));
    if (filter); dC = filterApply(dC,filt); end
    %% Compute the volume and sensitivity
    Ve = 4*a*b*t; V = sum(Ve*rhof);
    dV = Ve*ones(numberOfElements,1);
    if (filter); dV = filterApply(dV,filt); end
    %% Plot the current design
    plot2DElements(1,nodalCoordinates,elementConnectivity,u,[],rhof,1); 
    %% Convergence check
    fprintf('Iter. %3i - C = %4.3e, V = %4.3e, change = %3.2e \n',iter,C,V,change); 
    if (change < tol || iter == maxiter); break; end
    %% OC update
    rhonew = rho; rholow = max(0,rho-mvlim); rhoupp = min(1,rho+mvlim);
    s1 = 0; s2 = 1e12; ocfac = rho.*((-dC./dV).^0.5);
    while (s2-s1)/(s1+s2) > 1e-6
        smid = 0.5*(s1+s2);
        gre = ocfac/smid;
        rhonew = max(rholow,min(rhoupp,gre));
        if (filter); rhof = filterApply(rhonew,filt); else rhof = rhonew; end
        Vnew = sum(Ve*rhof);
        if Vnew > Vmax; s1 = smid; else; s2 = smid; end
    end
    change = max(rho-rhonew); rho = rhonew;
    elementStiffness = Emin + (Emax-Emin)*rhof.^p;
end
%% Measure of non-discreteness
Mnd = full(4*sum(rhof(:).*(1-rhof(:)))/numberOfElements*100)
figure(2); histogram(rho,10); hold on; histogram(rhof,10); legend('rho','rhof');