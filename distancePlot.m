function distancePlot(fileName, bin, cutoff)
    l = csvread(fileName, 0, 4);
    a = l(:, 1:2);% dist_1, distance between CA and side-chain mass center; dist_2, distance between CA and block mass center
    b = l(:, 3:4); % phi, psi
    l = [b a]; % use phi, psi, dist_1, dist_2
    l = l(:, 1:4); % use phi, psi, and chi_1
    binNum = 360/bin+1;
    bins = linspace(0, 360, binNum);
    bins = bins(2:end);
    index = 1;
    figure;
    for i = bins
        for j = bins
            t = l( (l(:, 1) >= i-bin) & (l(:, 1) < i) & (l(:, 2) >= j-bin) & (l(:, 2) < j), :); 
            dist_1 = t(:, 3); % get the column representing dist_1, which is the distance between CA and the side-chain mass center
            ax = subplot(size(bins, 2), size(bins, 2), index);
            index = index + 1;
            [bincounts] = histc(dist_1, 0:0.5:8);
            norm_1 = bincounts/max(bincounts);
            bar(0:0.5:8, norm_1, 'r');
            %histogram(ax, norm_1, 'Normalization', 'probability');
            axis(ax, [0 8 0 inf]);
            set(ax, 'XTick', 0:1:8);
            %set(ax, 'XTick', [], 'XTickLabel', []);
            %xlabel(ax, strcat('(', int2str(i), ',', int2str(j), ')'));
        end
    end
    figure;
    plot3k(l(:, 1:3), 'Labels', {'Side-Chain Mass Distribution', 'Phi', 'Psi', 'Distance'});
    index = 1;
    figure;
        for i = bins
        for j = bins
            t = l( (l(:, 1) >= i-bin) & (l(:, 1) < i) & (l(:, 2) >= j-bin) & (l(:, 2) < j), :); 
            dist_1 = t(:, 3); % get the column representing dist_1, which is the distance between CA and the side-chain mass center
            ax = subplot(size(bins, 2), size(bins, 2), index);
            index = index + 1;
            histogram(ax, dist_1);
            axis(ax, [0 8 0 cutoff]);
            set(ax, 'XTick', [], 'XTickLabel', []);
            xlabel(ax, strcat('(', int2str(i), ',', int2str(j), ')'));
        end
        end
        index = 1;
    figure;
    for i = bins
        for j = bins
            t = l( (l(:, 1) >= i-bin) & (l(:, 1) < i) & (l(:, 2) >= j-bin) & (l(:, 2) < j), :); 
            dist_2 = t(:, 4); % get the column representing dist_1, which is the distance between CA and the side-chain mass center
            ax = subplot(size(bins, 2), size(bins, 2), index);
            index = index + 1;
            [bincounts] = histc(dist_2, 0:0.5:8);
            norm_1 = bincounts/max(bincounts);
            bar(0:0.5:8, norm_1, 'r');
            %histogram(ax, norm_1, 'Normalization', 'probability');
            axis(ax, [0 8 0 inf]);
            set(ax, 'XTick', 0:1:8);
            %set(ax, 'XTick', [], 'XTickLabel', []);
            %xlabel(ax, strcat('(', int2str(i), ',', int2str(j), ')'));
        end
    end    
    figure;
    plot3k([l(:, 1:2) l(:, 4)], 'Labels', {'Block Mass Distribution', 'Phi', 'Psi', 'Distance'});
    index = 1;
    figure;
        for i = bins
        for j = bins
            t = l( (l(:, 1) >= i-bin) & (l(:, 1) < i) & (l(:, 2) >= j-bin) & (l(:, 2) < j), :); 
            dist_2 = t(:, 4); % get the column representing dist_1, which is the distance between CA and the side-chain mass center
            ax = subplot(size(bins, 2), size(bins, 2), index);
            index = index + 1;
            histogram(ax, dist_2);
            axis(ax, [0 8 0 cutoff]);
            set(ax, 'XTick', [], 'XTickLabel', []);
            xlabel(ax, strcat('(', int2str(i), ',', int2str(j), ')'));
        end
    end
end