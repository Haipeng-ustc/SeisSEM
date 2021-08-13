function Check_Stability_Before(dt, f0, xmin, xmax, zmin, zmax, nelem_x, nelem_z, NGLLX, NGLLZ, vp, vs)



% define some values
CFL        = 0.6;   % stability number, sqrt(3.0/8.0)
sample_num = 10;    % samping points in the average wavelength

% average grid size
x_average = (xmax - xmin) /(nelem_x*(NGLLX-1) + 1);
z_average = (zmax - zmin) /(nelem_z*(NGLLZ-1) + 1);
dx        = max(x_average, z_average);

% velocity.
v_min = min(min(vs(vs>0)), min(vp(:)));
v_max = max(vp(:));

% maximum frequency.
fmax = v_min / dx / sample_num;
if  f0 > fmax
    fprintf('\n');
    fprintf('Error: The frequency of source wavelet could cause numerical dispersion.\n');
    fprintf('Set frequency no more than %.2f Hz.\n', fmax);
    error('Source requency parameters error');
end

% maximum time step.
dtmax = CFL* dx/v_max;
if dt > dtmax
    fprintf('\n');
    fprintf('Error: The time step is too large and the simulation could become unstable.\n');
    fprintf('Set time step no more than %f s\n', dtmax);
    error('Time step error');
end



% we need NGLLX == NGLLZ. 
if NGLLX ~= NGLLZ
    error('NGLLX is not same as NGLLZ.\n');       
end 
    


end