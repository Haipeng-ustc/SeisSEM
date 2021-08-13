
function [x,w,h,hw]=GetGLL(ngll)

% x: Gauss-Lobatto-Legendre points,
% w: weights
% h: derivatives of Lagrange polynomials, H_ij = h'_i(xgll(j))
% hw: h * w
    
name = sprintf('gll_library/%s_%0.2u.tab','gll',ngll);
% get full path
pathstr = fileparts(mfilename('fullpath'));
name = fullfile(pathstr,name);

if ~exist(name,'file')
    error(sprintf('Data file %s does not exist',name))
end

fid=fopen(name);
data=fscanf(fid,'%f',[ngll,ngll+2]);
fclose(fid);

x=data(:,1);
w=data(:,2);
h=data(:,3:end)';  % tranpose is Very important

hw = zeros(ngll, ngll);
for k1 = 1:ngll
    for k2 = 1:ngll
        hw(k2,k1) = w(k2) * h(k2,k1);
    end
end

end