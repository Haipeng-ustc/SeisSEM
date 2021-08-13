function  Plot_Wavefield(displ, tri, geo_xz, rec_xz, src_xz, coord, time, plot_ratio, src_node, rec_node, it)

    rec_xz = rec_xz/1000;
    src_xz = src_xz/1000;
    geo_xz = geo_xz/1000;
    coord  = coord/1000;
    coord(2,:)  = coord(2,:)  - geo_xz(2,end);
    rec_xz(2,:) = rec_xz(2,:) - geo_xz(2,end);
    src_xz(2,:) = src_xz(2,:) - geo_xz(2,end);
    geo_xz(2,:) = geo_xz(2,:) - geo_xz(2,end);

    Fig = figure(2);        
    max_value = max(abs(displ(:))) * 5;

    subplot(1,2,1)
%     trisurf(tri, coord(1,:), coord(2,:), displ(1,:)); view(2); shading interp; hold on; box on
    scatter(coord(1,:), coord(2,:), 10, displ(1,:), 'filled'); grid on; hold on;

    plot(geo_xz(1,:), geo_xz(2,:), 'k-', 'LineWidth',2.0); 
    xlabel('Distance (km)'); ylabel('Depth (km)');  colormap('redblue'); 
    pbaspect([1 1/plot_ratio 1]); xlim([0, max(coord(1,:))]);
    caxis([-0.1 0.1]*max_value); set(gca, 'FontSize', 16);
    title(sprintf('Displacement-x at %.2f s', time))
    set(gcf,'position',[50,50,1800,600]);
    colorbar
    hold off
    

    subplot(1,2,2)
%   trisurf(tri, coord(1,:), coord(2,:), displ(2,:)); view(2); shading interp; hold on;box on
    scatter(coord(1,:), coord(2,:), 10, displ(2,:), 'filled'); grid on; hold on;
    plot(geo_xz(1,:), geo_xz(2,:), 'k-', 'LineWidth',2.0);
    xlabel('Distance (km)'); ylabel('Depth (km)');  colormap('redblue'); 
    pbaspect([1 1/plot_ratio 1]); xlim([0, max(coord(1,:))]);
    caxis([-0.1 0.1]*max_value); set(gca, 'FontSize', 16);
    title(sprintf('Displacement-z at %.2f s', time))
    set(gcf,'position',[50,50,1800,600]);
    colorbar
    hold off

%     figname = sprintf('./California_Observation_data_new/snap_%d.png',it);
%     print(gcf, figname, '-dpng', '-r300');
end



