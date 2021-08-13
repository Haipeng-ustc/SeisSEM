clear;clc;close all;

dx   = 100;  % approximated element size
xmin = 0.0;
xmax = 4000;
zmin = 0.0;

topo_z = load('./mesh/topo_data.dat');
% topo_z = [2000*ones(1,15),topo_z, 2000*ones(1,15)];
topo_z = topo_z(1:3:end);
topo_z = smooth(topo_z, 4)';
% topo_z = [4000, 4000];
topo_x = linspace(xmin, xmax, length(topo_z));
geo_x = [topo_x, xmax, xmin];
geo_z = [topo_z, zmin, zmin];
geo_n = length(geo_x);
figure 
scatter(geo_x/1000, geo_z/1000-2,'ro'); hold on;
plot(geo_x/1000, geo_z/1000-2, 'b-', 'LineWidth',2.0); hold on;
plot([xmin/1000, zmin/1000], [geo_x(1)/1000-2  , geo_z(1)/1000-2], 'b-', 'LineWidth',2.0);
xlabel('Distance (m)'); ylabel('Depth (m)'); set(gca,'FontSize', 16);
pbaspect([1 1/1.6443 1])

height = max(geo_z) - min(geo_z);
width  = max(geo_x) - min(geo_x);

fp = fopen('./mesh/mymesh.geo','w');
fprintf(fp, '//Define geometry\n', dx);
fprintf(fp, 'dx=%f;\n', dx);
fprintf(fp, 'height=%f;\n', height);
fprintf(fp, 'width=%f;\n', width);

for i =1 : geo_n
    fprintf(fp, 'Point(%d) = { %f, %f, 0.0, dx};\n', i, geo_x(i), geo_z(i));
end 

% Topography interface, Line 1
for i = 1 : geo_n - 2
    if i ==1
        fprintf(fp, '\nLine(1) = {%d, ', i);
    elseif i < length(geo_x) - 2
        fprintf(fp, '%d, ', i);
    else 
        fprintf(fp, '%d};\n',i);
    end 
end
fprintf(fp, 'Line(2) = {%d, %d};\n', geo_n - 2, geo_n - 1);
fprintf(fp, 'Line(3) = {%d, %d};\n', geo_n - 1, geo_n    );
fprintf(fp, 'Line(4) = {%d, %d};\n', geo_n    , 1        );

fprintf(fp, '\nTransfinite Line {4, 2} = Ceil(height/dx) Using Progression 1;\n');
fprintf(fp, 'Transfinite Line {1, 3} = Ceil(width/dx) Using Progression 1;\n');

fprintf(fp, '\nLine Loop(1) = {1, 2, 3, 4};\n');
fprintf(fp, 'Plane Surface(2) = {1};\n');
fprintf(fp, 'Transfinite Surface {2};\n');
fprintf(fp, 'Recombine Surface {2};\n');

fprintf(fp, '\n// quads mesh;\n');
fprintf(fp, 'Mesh.SubdivisionAlgorithm = 1;\n');
% 4-node element, it should be the same with polynomial degree (inside each element, along each direction) 
fprintf(fp, 'Mesh.ElementOrder = %d;\n', 4); 
fclose(fp);
% More advanced meshing is under construction






