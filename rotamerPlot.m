function rotamerPlot(fileName, bin, cutoff)
    l = csvread(fileName, 0, 4);
    l = l(:, 3:5); % use phi, psi, and chi_1
    binNum = 360/bin+1;
    bins = linspace(0, 360, binNum);
    bins = bins(2:end);
    index = 1;
    figure;
    for i = bins
        for j = bins
            t = l( (l(:, 1) >= i-bin) & (l(:, 1) < i) & (l(:, 2) >= j-bin) & (l(:, 2) < j), :); 
            chi_1 = t(:, 3); % get the column representing chi_1
            ax = subplot(size(bins, 2), size(bins, 2), index);
            index = index + 1;
            %histogram(chi_1);
            %figure
            [bincounts] = histc(chi_1, 0:10:360);
            if max(bincounts) < cutoff
                uplimit = cutoff;
            else
                uplimit = max(bincounts);
            end
            if uplimit > 2*cutoff
                bar(0:10:360, bincounts, 'r');
            else
                bar(0:10:360, bincounts, 'histc');
            end
            axis(ax, [0 360 0 uplimit]);
            set(ax, 'XTick', [], 'XTickLabel', []);
            xlabel(ax, strcat('(', int2str(i), ',', int2str(j), ')'));
        end
    end
    figure;
    plot3k(l, 'Labels', {'Rotamer Density', 'Phi', 'Psi', 'Chi'});
end