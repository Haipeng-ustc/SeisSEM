
function tri_used = Plot_tri_used(coord)

xx_grid_fine = coord(1,:);
zz_grid_fine = coord(2,:);
zz_elev_fine = coord(2,:);

tri      = delaunay(xx_grid_fine, zz_grid_fine);
tri_used = zeros(size(tri));

use_count = 0;
for ii=1:size(tri, 1)
    flag = 0;
    dis1 = sqrt((coord(1,tri(ii,1)) - coord(1,tri(ii,2)))^2 +  (coord(2,tri(ii,1)) - coord(2,tri(ii,2)))^2);
    dis2 = sqrt((coord(1,tri(ii,2)) - coord(1,tri(ii,3)))^2 +  (coord(2,tri(ii,2)) - coord(2,tri(ii,3)))^2);
    dis3 = sqrt((coord(1,tri(ii,1)) - coord(1,tri(ii,3)))^2 +  (coord(2,tri(ii,1)) - coord(2,tri(ii,3)))^2);
    if max([dis1, dis2, dis3]) < 50
        use_count = use_count + 1;
        tri_used(use_count, :) = tri(ii, :);
    end 
end

tri_used(use_count+1:end, :) = [];

end