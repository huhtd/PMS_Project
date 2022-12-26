% 2D FEM analysis of the waveguide discontinuities 

        % 
        % The goal of this analysis is to solve for the electric and magnetic fields and the Poynting vector in a 2D waveguide with a discontinuity. The waveguide is assumed to have an even mode of a specific order. The permittivity and distance between the waveguide plates are specified as input parameters, as well as the frequency of the waveguide. The input and output branches of the waveguide have specified permittivities, and the permittivity of the waveguide discontinuity is also specified.
        %
        % Inputs:
        % - mesh: ASCII data for the mesh of the waveguide, generated using COMSOL
        % - boundary nodes: ASCII data for the boundary nodes of the waveguide, generated using COMSOL
        % - domains: ASCII data for the domains of the waveguide, generated using COMSOL
        % - d: distance of the waveguide plates
        % - f: frequency of the waveguide
        % - Order: order of the even mode
        % - eps1: permittivity of the input branch of the waveguide
        % - eps2: permittivity of the waveguide discontinuity
        % - eps3: permittivity of the output branch of the waveguide
        % - solver: flag to specify the solution method to be used:
        % 1 - direct sparse matrix solver
        % 2 - GMRES iterative solver
        %
        % Outputs:
        % - Ez: electric field in the waveguide
        % - Hx: x-component of the magnetic field in the waveguide
        % - Hy: y-component of the magnetic field in the waveguide
        % - Smod: magnitude of the Poynting vector in the waveguide
        % - S11_input: input reflection coefficient
        % - S21_output: output transmission coefficient
        %



% Load ASCII data: mesh, boundary nodes, domains (COMSOL is used as a mesh
% generator)
load_data;

% Parameters definition
% **************************************************************
d=0.02;                     % Distance of the waveguide plates 
f=10e9;                     % Frequency
Order=0;                    % Order of the even mode
eps1=1.;                    % Permittivity of the input branch of the waveguide
eps2=4.;                    % Permittivity of the waveguide discontinuity    
eps3=1.;                    % Permittivity of the output branch of the waveguide
imj=sqrt(-1);
mu0=4*pi*1.e-7;             % Mag. permeability of vacuum
eps0=8.85e-12;              % Diel. permittivity of vacuum
c=1/sqrt(mu0*eps0);
omega=2*pi*f;               % Angular frequency
k0=omega*sqrt(mu0*eps0);    % Free space wave-vector
ky=(2*Order+1)*pi/d;        % Even modes
kx=sqrt(k0^2-ky^2);         % Guiding mode propagation constant
% E0=1.;                                % Amplitude of the input signal
E0=sqrt(2)*sqrt(2*omega*mu0/(kx*d));    % Normalized input signal:Pin=1
solver=1;                   % 1 - sparse matrix solver; 2 - GMRES iterative solver
% *************************************************************************

[x_min_bc,y_min_bc,y_max_bc,x_max_bc]=find_boundaries(Nn,x_no,y_no);
        % This function finds the minimum and maximum x and y coordinates of the boundary nodes
        %
        % Inputs:
        % - Nn: number of nodes in the mesh
        % - x_no: x-coordinates of the nodes
        % - y_no: y-coordinates of the nodes
        %
        % Outputs:
        % - x_min_bc: minimum x-coordinate of the boundary nodes
        % - y_min_bc: minimum y-coordinate of the boundary nodes
        % - y_max_bc: maximum y-coordinate of the boundary nodes
        % - x_max_bc: maximum x-coordinate of the boundary nodes
        %

boundary_nodes;             % This function excludes the nodes between two different materials from the set of boundary nodes
        % This function excludes the nodes between two different materials from the set of boundary nodes
        %
        % Inputs:
        % - Nn: number of nodes in the mesh
        % - x_no: x-coordinates of the nodes
        % - y_no: y-coordinates of the nodes
        % - x_min_bc: minimum x-coordinate of the boundary nodes
        % - y_min_bc: minimum y-coordinate of the boundary nodes
        % - y_max_bc: maximum y-coordinate of the boundary nodes
        % - x_max_bc: maximum x-coordinate of the boundary nodes
        % - material: material domain for each node
        %
        % Outputs:
        % - boundary_nodes: set of boundary nodes
        %

figure(1);
set_figure_1;
mesh_plot;
bcs_plot;
axis on;

asm_matrix_wp;  % Matrix assembly % This function assembles the global stiffness matrix for the waveguide

def_ports;      % Ports definition

def_bcs;        % Definition of the PEC boundary conditions

switch solver
            % ***************************
    case 1  % direct sparse matrix solver
            % ***************************
        spparms('autoamd',0); %we do not want autoamd, we want colamd
        permut=colamd(A); % Matrix reordering (bandwidth reduction)
        Vp = A (permut,permut) \ b(permut); % Direct solution
        Vc(permut)=Vp; % Mapping to the original numbering  
            % **********************
    case 2  % GMRES iterative solver
            % **********************
        M1=sparse(Nn,Nn);
        M=diag(A);
        for p=1:Nn
           M1(p,p)=M(p); 
        end
        [Vc,flag,relres,iter]=gmres(A,b,[],1e-7,Nn,M1);
end

Ez=Vc; % Electric field

figure(2);
V=real(Vc);
set_figure_1;
find_min_max;
field_plot;
geometry_plot;
title('Real part of Ez (V/m)');

% Magnetic field computation: Hx, Hy
magnetic_field;
figure(3);
set_figure_1;
field_plot;
geometry_plot;
quiver(x_no,y_no,real(Hx),real(Hy),1,'k'); % Plot real part of the magnetic field

% Poynting vector
        % This function computes the Poynting vector in the waveguide
        %
        % Inputs:
        % - Nn: number of nodes in the mesh
        % - x_no: x-coordinates of the nodes
        % - y_no: y-coordinates of the nodes
        % - tri_no: list of nodes for each element
        % - Vc: nodal voltages
        % - Hx: x-component of the magnetic field
        % - Hy: y-component of the magnetic field
        % - mu0: mag. permeability of vacuum
        % Outputs:
        % - Sx: x-component of the Poynting vector
        % - Sy: y-component of the Poynting vector
        % - Smod: magnitude of the Poynting vector
        %
Poynting_field;
figure(4);
set_figure_1;
V=abs(Smod);    % Magnitude of the Poynting vector
Vc=Smod;        % Poynting vector
find_min_max;   % This function finds the minimum and maximum values of the Poynting vector
min_field=0;    %minimum value of the Poynting vector
field_plot;
geometry_plot;

% Sparameters
S11_input;
        % This function computes the input reflection coefficient
        %
        % Inputs:
        % - port1: set of nodes for port 1
        % - Vc: nodal voltages
        % - E0: normalized input signal
        % - mu0: mag. permeability of vacuum
        % - kx: guiding mode propagation constant
        %
        % Outputs:
        % - S11: input reflection coefficient
        %
S21_output;
        % This function computes the output transmission coefficient
        %
        % Inputs:
        % - port2: set of nodes for port 2
        % - Vc: nodal voltages
        % - E0: normalized input signal
        % - mu0: mag. permeability of vacuum
        % - kx: guiding mode propagation constant
        %
        % Outputs:
        % - S21: output transmission coefficient
return;
