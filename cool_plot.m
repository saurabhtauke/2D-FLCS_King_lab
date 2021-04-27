function  cool_plot(T1, T2, F_ILT)  

        theCM = parula;
        theCM(1,:) = 1;
            
    pcolor(T1, T2, F_ILT(1:length(T2), 1:length(T1)));
        colormap(theCM);
        shading interp;
        colorbar;
 
	axis square
	h = gca;
	set(h, 'XScale', 'log', 'YScale', 'log')
	xlabel('\tau_1 (ns)', 'FontSize', 14,'FontWeight','bold')
	ylabel('\tau_2 (ns)', 'FontSize', 14,'FontWeight','bold')
	set(gca, 'FontSize', 14,'FontWeight','bold');
    set(gca, 'TickDir','in')
    set(gca, 'Box','on')
    set(gca, 'LineWidth',1.5)
    xlim ([min(T1)-min(T1)/50 max(T1)+max(T1)/50])
    ylim ([min(T2)-min(T2)/50 max(T2)+max(T2)/50])
   
end