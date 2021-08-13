% SEM 2D Elastic wave code.
%
% by Haipeng Li, Hongbo Han, and Wancong Zhao. 
% at University of Science and Technology of China
% Contact: haipengl@mail.ustc.edu.cn


clear; clc; close all;
addpath(genpath(pwd), '-begin'); warning('off');

%%%%%%% Model Size Pars %%%%%%%%
NDIM = 2;       % 2-D simulation
xmin = 0;       % xmin
zmin = 0;       % zmin
xmax = 4000;    % xmax
zmax = 2000;    % zmax
nelem_x = 80;   % number of spectral elements along X
nelem_z = 40;   % number of spectral elements along Z
P     = 4;      % polynomial degree (inside each element, along each direction)
NGLLX = P + 1;  % number of GLL integration points in x direction of an element
NGLLZ = NGLLX;  % number of GLL integration points in z direction of an element
ngnod = 4;      % use 4-node elements to describe the geometry of the mesh, other option like 9-node

%%%%%%% Time Evolution Pars %%%%%%%%
dt = 5e-4;
nt = 2001;
t  = linspace(0, (nt-1)*dt, nt);
time        = dt * nt;
threshold   = 1.e+20; % displacement threshold above which we consider that the code became unstable
output_step = 100;

%%%%%%%% PML Pars %%%%%%%%
use_pml_or_not = true;
pml_x_thick  = 300;
pml_z_thick  = 300;
free_surface = 1;

%%%%%%%% Source & Receiver Pars %%%%%%%%
%% Source time function
f0 = 15;
t0 = 1.2 / f0;
amplitude   = 1.e10;
src_wavelet = amplitude * Source_Time_Function(t,'ricker',f0, t0);

%% Source type: choose between 'vector_source' and 'moment_tensor_source'
src_type = 'vector_source'; 
% for vector_source
force_angle = 90; 
% for moment_tensor_source
strike      = 45;
dip         = 30; 
rake        = 104; 
Moment_Tensor = Moment_Tensor_Solution(strike, dip, rake); 
%% Source location
src_n       = 1;
src_xz      = zeros(2, src_n);
src_xz(1,:) = 2000;
src_xz(2,:) = 2000;

%%%%%%%%%%%%%%%%% Mesh Regular Model %%%%%%%%%%%%%%%%%
% Mesh information
% NSPEC = nelem_x * nelem_z;                                 % number of spectral elements of the mesh
% NGLOB = (nelem_x*(NGLLX-1) + 1) * (nelem_z*(NGLLZ-1) + 1); % number of unique grid points of the mesh
% NPGEO = (nelem_x+1)*(nelem_z+1);                           % total number of geometrical points that describe the geometry 
% % Geometrical points that describe the mesh:
% [coorg, knods] = Mesh_Geometry(xmin, xmax, zmin, zmax, nelem_x, nelem_z, 'q4'); % only q4 at present.  coorg and knods: coordinates and numbering of control nodes
% % Calculation points that describe the mesh:
% [ibool] = Mesh_Global(nelem_x, nelem_z, NGLOB, NGLLX, NSPEC); % ibool: numbering of all the calculation points
% % for plot
% geo_x  = [xmin, xmax, xmax, xmin, xmin];
% geo_z  = [zmin, zmin, zmax, zmax, zmin];
% geo_xz = [geo_x;geo_z];
% %% Receiver location
% rec_n       = 21;
% rec_xz      = zeros(2, rec_n);
% rec_xz(1,:) = linspace(xmin, xmax, rec_n);
% rec_xz(2,:) = 2000;

%%%%%%%%%%%%%%%%% Mesh Irregular Model %%%%%%%%%%%%%%%%%
[coorg, knods, ibool] = Read_Gmsh('./mesh/mymesh.msh', NGLLX, NGLLZ);
NSPEC = size(ibool,3);                       % number of spectral elements of the mesh
NGLOB = size(coorg,2);                       % number of unique grid points of the mesh
NPGEO = size(knods,1) * size(knods,2);       % total number of geometrical points that describe the geometry
% for plot
topo_z = load('./mesh/topo_data.dat');topo_z = topo_z(1:3:end);
topo_z = smooth(topo_z, 4)';
topo_x = linspace(xmin, xmax, length(topo_z));
geo_x  = [topo_x, xmax, xmin, xmin];
geo_z  = [topo_z, zmin, zmin, topo_z(1)];
geo_xz = [geo_x;geo_z];

%% Receiver location
rec_n       = length(topo_z);
rec_xz      = zeros(2, rec_n);
rec_xz(1,:) = topo_x;
rec_xz(2,:) = topo_z;

%%%%%%%%%%%%%%%%% Assemble Global Jacobian %%%%%%%%%%%%%%%%%
[coord, xix, xiz, gammax, gammaz, jacobian] = Assemble_Global_Jacobian(ibool,coorg,knods,NGLLX,NGLLZ);

%%%%%%%%%%%%%%%%% Model Medium Pars %%%%%%%%%%%%%%%%%
% Case 1: isotropic media
rho = ones(NGLOB, 1) * 2700;
vp  = ones(NGLOB, 1) * 5000;
vs  = vp / 1.732;
mu     = rho .* vs .* vs;
lambda = rho .* vp .* vp - 2.0 * mu;
lambdaplus2mu = lambda + 2.0 * mu;
c11 = lambdaplus2mu;
c13 = lambda;
c33 = lambdaplus2mu;
c44 = mu;

% Case 2: anisotropic media zinc
% scale_aniso = 1e9;
% rho = ones(NGLOB, 1) * 7100;
% c11 = ones(NGLOB, 1) * 165.0 * scale_aniso;
% c13 = ones(NGLOB, 1) * 50.0  * scale_aniso;
% c33 = ones(NGLOB, 1) * 62.0  * scale_aniso;
% c44 = ones(NGLOB, 1) * 39.6  * scale_aniso;
% vp  = sqrt(c33./rho);
% vs  = sqrt(c44./rho);


%%%%%%%%%%%%%%%%% PML damping profiles %%%%%%%%%%%%%%%%%
[pml_dx, pml_dz, pml_dxx, pml_dzz] = Set_PML( vp, xmin, xmax, zmin, zmax, NSPEC, NGLOB, NGLLX, NGLLZ, pml_x_thick, pml_z_thick, free_surface, coord, ibool);
if use_pml_or_not==false  % if pml is not used, set the dampping profiles to zeros.
    pml_dx(:) = 0.0; pml_dz(:) = 0.0; pml_dxx(:) = 0.0; pml_dzz(:) = 0.0;
end 

%%%%%%%%%%%%%%%%% Build the mass matrix %%%%%%%%%%%%%%%%%
[M, Cx, Cz, Cxx, Cxz, Czz] = Assemble_Global_M_and_C(NSPEC, NGLOB, NGLLX, NGLLZ, jacobian, ibool, rho, pml_dx, pml_dz);

%%%%%%%%%%%%%%%%% Build the stiffness matrix, with elastic tensor and PML %%%%%%%%%%%%%%%%%
[Kx1, Kx2, Kx3, Kx4, Kx5, Kz1, Kz2, Kz3, Kz4, Kz5] = Assemble_Global_K(NGLOB, NGLLX, NGLLZ, jacobian, gammax, gammaz, xix, xiz, ibool, c11, c13, c33, c44, pml_dxx, pml_dzz);
                                                                            
%%%%%%%% Constant value %%%%%%%%
a0 = 1 / dt^2; a1 = 1 / (2 * dt); a2 = 2 * a0; a3 = 1 / a2;

A1  =  a0 * M + a1 * 2 * Cx;
A2  =  a0 * M + a1 * ( Cx + Cz );
A3  =  a0 * M + a1 * 2 * Cz;
AL  =  2 * a1 * M;
% invert the these matrices here once and for all
A1_inv = zeros(NGLOB,1);
A2_inv = zeros(NGLOB,1);
A3_inv = zeros(NGLOB,1);
AL_inv = zeros(NGLOB,1);
M_inv  = zeros(NGLOB,1);
for i = 1:NGLOB
    A1_inv(i) = 1.0 / A1(i);
    A2_inv(i) = 1.0 / A2(i);
    A3_inv(i) = 1.0 / A3(i);
    AL_inv(i) = 1.0 / AL(i);
    M_inv(i)  = 1.0 / M(i);
end

%%%%%%%% Wavefield plot Pars %%%%%%%%
xx = linspace(xmin, xmax, nelem_x*(NGLLX-1) + 1);
zz = linspace(zmin, zmax, nelem_z*(NGLLZ-1) + 1);
[xGrid, zGrid] = meshgrid(xx,zz);
tri_used = Plot_tri_used(coord);
plot_ratio= max(coord(1,:)) / max(coord(2,:));

%%%%%%%% Relocate S & R %%%%%%%%
[src_xz, src_node, ~] = Find_Nearest_Node(src_xz,coord);
[rec_xz, rec_node, ~] = Find_Nearest_Node(rec_xz,coord);

%%%%%%%% Check the simulation is OK %%%%%%%%
Check_Stability_Before(dt, f0, xmin, xmax, zmin, zmax, nelem_x, nelem_z, NGLLX, NGLLZ, vp, vs);

%%%%%%%% Set source arrays, this is how to add source in SEM %%%%%%%%
[sourcearray, ispec_source] = Set_Source_Array(src_node, src_type, force_angle, Moment_Tensor, ibool, xix, xiz, gammax, gammaz, NDIM, NGLLX, NGLLZ);

% plot the velocity model
figure(1)
trisurf(tri_used, coord(1,:)/1000, (coord(2,:) - geo_xz(2,end))/1000, vp(:)); view(2); shading interp; hold on;
scatter3(rec_xz(1,:)/1000, (rec_xz(2,:)-geo_xz(2,end))/1000,vp(rec_node),100,'rv', 'filled');hold on;
scatter3(src_xz(1,:)/1000, (src_xz(2,:)-geo_xz(2,end))/1000,vp(src_node),400,'bp', 'filled'); 
xlabel('Distance (km)'); ylabel('Depth (km)'); pbaspect([1 1/plot_ratio 1]); set(gca, 'FontSize', 16); 
hcb = colorbar; title(hcb,{'V_{s} (m/s)'}, 'FontName','Times New Roman','FontSize',10);
set(gcf,'position',[50,50,1000,600])
print(gcf, './results/vp.png', '-dpng', '-r300');

%%%%%%%% Initial arrays %%%%%%%%
seism = zeros(NDIM,rec_n,nt); % seismogram recorded at the receiver
displ_now       = zeros(NGLOB, NDIM); % global displacement vector
displ_old_term1 = zeros(NGLOB, NDIM); % global displacement vector 
displ_old_term2 = zeros(NGLOB, NDIM); % global displacement vector 
displ_old_term3 = zeros(NGLOB, NDIM); % global displacement vector 
displ_now_term1 = zeros(NGLOB, NDIM); % global displacement vector 
displ_now_term2 = zeros(NGLOB, NDIM); % global displacement vector 
displ_now_term3 = zeros(NGLOB, NDIM); % global displacement vector 
displ_new_term1 = zeros(NGLOB, NDIM); % global displacement vector 
displ_new_term2 = zeros(NGLOB, NDIM); % global displacement vector 
displ_new_term3 = zeros(NGLOB, NDIM); % global displacement vector 
p_x             = zeros(NGLOB, NDIM);  
p_z             = zeros(NGLOB, NDIM);    
energy          = zeros(nt,1);  

% start of the time loop
for it = 1 : nt
    % compute current time
    time = (it-1)*dt;
    % compute maximum of norm of displacement from time to time and display it in order to monitor the simulation
    if (mod(it,output_step) == 0 || it == 1 || it == nt)
        Check_Stability_After(displ_now', NGLOB, it, dt, nt, threshold);
        Plot_Wavefield(displ_now', tri_used, geo_xz, rec_xz, src_xz, coord, time, plot_ratio, src_node, rec_node, it);
    end
    energy(it) = sum(displ_now(:,1).*displ_now(:,1) + displ_now(:,2).*displ_now(:,2)); 
    % Time evlolution scheme. 
    displ_new_term1(:,1)  =  A1_inv.* (M .* p_x(:,1) - Kx1 * displ_now(:,1) - ( Cxx - a2 * M ) .* displ_now_term1(:,1) - ( a0 * M - a1 * 2 * Cx ) .* displ_old_term1(:,1));
    displ_new_term2(:,1)  =  A2_inv.* (              - Kx2 * displ_now(:,2) - ( Cxz - a2 * M ) .* displ_now_term2(:,1) - ( a0 * M - a1 * (Cx+Cz)) .* displ_old_term2(:,1));
    displ_new_term3(:,1)  =  A3_inv.* (M .* p_x(:,2) - Kx3 * displ_now(:,1) - ( Czz - a2 * M ) .* displ_now_term3(:,1) - ( a0 * M - a1 * 2 * Cz ) .* displ_old_term3(:,1));
    p_x(:,1)              =  AL_inv.* (              - Kx4 * displ_now(:,1) - Cx .* p_x(:,1) + 2 * a1 * M .* p_x(:,1));
    p_x(:,2)              =  AL_inv.* (              - Kx5 * displ_now(:,1) - Cz .* p_x(:,2) + 2 * a1 * M .* p_x(:,2));
    
    displ_new_term1(:,2)  =  A1_inv.* (M .* p_z(:,1) - Kz1 * displ_now(:,2) - ( Cxx - a2 * M ) .* displ_now_term1(:,2) - ( a0 * M - a1 * 2 * Cx ) .* displ_old_term1(:,2));
    displ_new_term2(:,2)  =  A2_inv.* (              - Kz2 * displ_now(:,1) - ( Cxz - a2 * M ) .* displ_now_term2(:,2) - ( a0 * M - a1 * (Cx+Cz)) .* displ_old_term2(:,2));
    displ_new_term3(:,2)  =  A3_inv.* (M .* p_z(:,2) - Kz3 * displ_now(:,2) - ( Czz - a2 * M ) .* displ_now_term3(:,2) - ( a0 * M - a1 * 2 * Cz ) .* displ_old_term3(:,2));
    p_z(:,1)              =  AL_inv.* (              - Kz4 * displ_now(:,2) - Cx .* p_z(:,1) + 2 * a1 * M .* p_z(:,1) );
    p_z(:,2)              =  AL_inv.* (              - Kz5 * displ_now(:,2) - Cz .* p_z(:,2) + 2 * a1 * M .* p_z(:,2));
   
    displ_old_term1 = displ_now_term1;
    displ_old_term2 = displ_now_term2;
    displ_old_term3 = displ_now_term3;

    displ_now_term1 = displ_new_term1;
    displ_now_term2 = displ_new_term2;
    displ_now_term3 = displ_new_term3;
    displ_now = displ_now_term1 + displ_now_term2 + displ_now_term3;
    displ_now(src_node,1) = displ_now(src_node,1) + sin(force_angle*pi/180) * src_wavelet(it);
    displ_now(src_node,2) = displ_now(src_node,2) + cos(force_angle*pi/180) * src_wavelet(it);

    
    % add the source here
%     for j = 1:NGLLZ
%         for i = 1:NGLLX
%             iglob = ibool(i,j,ispec_source);
%             displ_now(iglob,:) = displ_now(iglob,:) - sourcearray(:,i,j)' * src_wavelet(it);
%         end
%     end
    
%     displ_new_x = dt * dt * M_inv .* (force_x - (Kx1+Kx3) * displ_now_x - Kx2*displ_now_z) + 2 * displ_now_x - displ_old_x;
%     displ_new_z = dt * dt * M_inv .* (force_z - (Kz1+Kz3) * displ_now_z - Kz2*displ_now_x) + 2 * displ_now_z - displ_old_z;
% 
%     displ_old_x = displ_now_x;
%     displ_old_z = displ_now_z;
% 
%     displ_now_x = displ_new_x;
%     displ_now_z = displ_new_z;
    
    % record a seismogram to check that the simulation went well
    seism(:,:,it) = displ_now(rec_node,:)';
end

figure
Plot_Trace(rec_xz, t, seism, 1.2, './results/seismgram.png');



