function  Plot_Wavefield_For_DB(displ, tri, geo_xz, rec_xz, src_xz, coord, time, plot_ratio, src_node, rec_node, it, seism,t,vp)

rec_xz = rec_xz/1000;
src_xz = src_xz/1000;
geo_xz = geo_xz/1000;
coord  = coord/1000;
coord(2,:)  = coord(2,:)  - geo_xz(2,end);
rec_xz(2,:) = rec_xz(2,:) - geo_xz(2,end);
src_xz(2,:) = src_xz(2,:) - geo_xz(2,end);
geo_xz(2,:) = geo_xz(2,:) - geo_xz(2,end);

h = figure(2);
max_value = max(abs(displ(:))) * 5;

set(h, 'Position', [10, 10, 800, 1000]);       % <— ‘Height’ Increased

subplot(2,1,1)
trisurf(tri, coord(1,:), coord(2,:), displ(2,:)); view(2); shading interp; hold on; box on
% plot(geo_xz(1,:), geo_xz(2,:), 'k-', 'LineWidth',0.5);
% scatter3(rec_xz(1,1:5:end), (rec_xz(2,1:5:end)),vp(rec_node(1:5:end)),100, 'kv', 'filled');hold on;
% scatter3(src_xz(1,      :), (src_xz(2,      :)),vp(src_node),         400,'rp', 'filled'); 
set(gca,'FontName','Times New Roman','FontSize',18);
xticks([0 1 2 3 4]);
% xticklabels({'0.5','0'})
yticks([1.5 2.0]);
yticklabels({'0.5','0'})

% xlabel('Distance (km)'); 
ylabel('Depth (km)');  
colormap('redblue');
pbaspect([1 1/8 1]);
caxis([-1 1]*1e7)
ylim([1.5, 2])
set(gca, 'FontSize', 16);
title(sprintf('Displacement Wavefield', time));
set(gca, 'FontName', 'Times New Roman');
% set(gca, 'Ydir', 'reverse');
hold off
hsp1 = get(gca, 'Position');                     % Get 'Position' for (2,1,1)

subplot(2,1,2)
imagesc(rec_xz(1,:), t, squeeze(seism(2,:,:))');
xlabel('Distance (km)'); ylabel('Time (s)');  colormap('redblue');
set(gca,'FontName','Times New Roman','FontSize',18,'Layer','top','XTick',[0 1 2 3 4], 'YTick',[0 0.5 1 1.5]);
caxis([-1 1]*1e7)
pbaspect([1 1/plot_ratio 1]);
set(gca, 'FontSize', 16);
title(sprintf('Recorded Seismogram'))
box on
hold off
hsp2 = get(gca, 'Position');                     % Get 'Position' for (2,1,2)
set(gca, 'Position', [hsp2(1:3)  2*hsp1(4)])    % Use 2*(2,1,1) Height for (2,1,2)
set(gca, 'FontName', 'Times New Roman')
figname = sprintf('./db/snap_%d.png',it);
print(gcf, figname, '-dpng', '-r300');
end



