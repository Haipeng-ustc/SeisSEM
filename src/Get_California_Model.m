function [vp_model, vs_model, rho_model] = Get_California_Model(depth_in_m)
    if depth_in_m >= 0.0
        depth_in_m = 0.0;
    else 
        depth_in_m = - depth_in_m;
    end 

    thick = [ 0.2, 0.03, 10.3, 11.66, 10.36, 300] * 1000;
    vp    = [3.81, 2.50, 6.10,  6.51,  6.90, 7.5] * 1000;
    vs    = [1.94, 1.07, 3.53,  3.71,  3.93, 4.5] * 1000;
    rho   = [0.92, 2.11, 2.74,  2.83,  2.92, 3.0] * 1000;
    depth_layer = zeros(size(thick));
    for i = 1 : length(thick)
        depth_layer(i) = sum(thick(1:i));
    end
    layer_index = min(find(depth_in_m < depth_layer(:)));
    vp_model  = vp(layer_index);
    vs_model  = vs(layer_index);
    rho_model = rho(layer_index);

end