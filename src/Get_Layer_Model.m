function [vp_model, vs_model, rho_model] = Get_Layer_Model(depth_in_m)
    depth_in_m = depth_in_m - 2000;
    if depth_in_m >= 0.0
        depth_in_m = 0.0;
    else 
        depth_in_m = - depth_in_m;
    end 

    thick = [0.2, 0.5,   1.0,   1.5,   2.0] * 1000;
    vp    = [3.5, 4.5, 5.5,  6.5,  7.0] * 1000;
    vs    = vp/1.732;
    rho   = [0.92, 2.11, 2.74,  2.83,  2.92] * 1000;
    depth_layer = zeros(size(thick));
    for i = 1 : length(thick)
        depth_layer(i) = sum(thick(1:i));
    end
    layer_index = min(find(depth_in_m < depth_layer(:)));
    vp_model  = vp(layer_index);
    vs_model  = vs(layer_index);
    rho_model = rho(layer_index);

end