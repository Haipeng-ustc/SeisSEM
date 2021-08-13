%FindNearestNode finds the nearest mesh node to a requested location (2D)

function [xz_out,iglob,dist] = FindNearestNode(xz_in, coord)


n     = size(xz_in,2);
dist  = zeros(n,1);
iglob = zeros(n,1);

for k=1:n
  [dist(k),iglob(k)] = min( (coord(1,:)-xz_in(1,k)).^2+(coord(2,:)-xz_in(2,k)).^2 );
end

dist  = sqrt(dist);
xz_out = coord(:,iglob); 

end 
