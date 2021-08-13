function Plot_Trace(rec_xz, t, seism, factor, figpath)
    
    rec_x   = rec_xz(1,:);
    seism_x = squeeze(seism(1,:,:));
    seism_z = squeeze(seism(2,:,:));

    [rec_n, nt] = size(seism_x);
   
    subplot(121)
    ampli = factor * max(abs( seism_x(:) ));
    if ampli>0
        offset = max(abs(diff(rec_x)));
        ampli = offset/ampli;
    end
    plot(t, seism_x * ampli + reshape(repmat(rec_x,1,nt), [rec_n, nt]), '-k', 'LineWidth', 1.5);
    title('Displacement-x'); xlabel('Time (s)'); ylabel('Distance (m)');
    set(gca, 'FontSize', 16);

    subplot(122)
    ampli = factor*max(abs( seism_z(:) ));
    if ampli>0
        offset = max(abs(diff(rec_x)));
        ampli = offset/ampli;
    end
    plot(t, seism_z*ampli + reshape(repmat(rec_x,1,nt), [rec_n, nt]), '-k', 'LineWidth', 1.5);
    title('Displacement-z'); xlabel('Time (s)'); ylabel('Distance (m)');
    set(gca, 'FontSize', 16);
    set(gcf,'position',[50,50,1600,800]);
    print(gcf, figpath, '-dpng', '-r300');

end