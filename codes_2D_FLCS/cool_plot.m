function  cool_plot(T1, T2, F_ILT)  
% plots a pseudo color optimized for 2D-ILT

        theCM = parula;
        theCM(1,:) = 1;
            
    pcolor(T1, T2, F_ILT(1:length(T2), 1:length(T1)));
        colormap(theCM);
%         caxis([0.00 1]);
        warning('check: custom coloraxis limits set')
        shading interp;
%         colorbar;
        
	axis square
    title('2D-ILT')
	h = gca;
	set(h, 'XScale', 'log', 'YScale', 'log')
	xlabel('\tau_1 (ns)', 'FontSize', 14,'FontWeight','bold')
	ylabel('\tau_2 (ns)', 'FontSize', 14,'FontWeight','bold')
	set(gca, 'FontSize', 14,'FontWeight','bold');
    set(gca, 'TickDir','in')
    set(gca, 'Box','on')
    set(gca, 'LineWidth',1.5)
    xticks([1E-1, 1E0, 1E1])
    yticks([1E-1, 1E0, 1E1])
    xlim ([min(T1)-min(T1)/50 max(T1)+max(T1)/50])
    ylim ([min(T2)-min(T2)/50 max(T2)+max(T2)/50])
   
end