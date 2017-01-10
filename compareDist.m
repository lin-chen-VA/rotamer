function rotamerPlot(res_1, res_2, phiAngleRange, psiAngleRange, xRange)
    % compare the side-chain distance distribution of proteins from X-Ray
    % and EM
    % For example:
    %       compareDist('X-RayR_1.5/ARG.csv', 'EM_3.5_4.0/ARG.csv', [240 300], [120 180], 2:0.01:7)
    phiLowAngle = phiAngleRange(1);
    phiHighAngle = phiAngleRange(2);
    psiLowAngle = psiAngleRange(1);
    psiHighAngle = psiAngleRange(2);
    l = getMatrix(res_1);
    l_2 = getMatrix(res_2);
    t = l( (l(:, 1) >= phiLowAngle) & (l(:, 1) < phiHighAngle) &(l(:, 2) >= psiLowAngle) & (l(:, 2) < psiHighAngle), :);     
    t_2 = l_2( (l_2(:, 1) >= phiLowAngle) & (l_2(:, 1) < phiHighAngle) &(l_2(:, 2) >= psiLowAngle) & (l_2(:, 2) < psiHighAngle), :);  
    %t_2 = l_2( (l_2(:, 1) >= lowAngle) & (l_2(:, 1) < highAngle), :); 
    dist_1 = t(:, 3); % get the column representing dist_1, which is the distance between CA and the side-chain mass center
    dist_2 = t_2(:, 3);
    histogram(dist_1, xRange, 'Normalization', 'pdf', 'FaceColor', 'r'); % red color for the first residue
    hold on;
    histogram(dist_2, xRange, 'Normalization', 'pdf', 'FaceColor', 'b'); % blue color for the second residue
    figure;
    dist_1 = t(:, 4); % get the column representing dist_1, which is the distance between CA and the side-chain mass center
    dist_2 = t_2(:, 4);
    histogram(dist_1, xRange, 'Normalization', 'pdf', 'FaceColor', 'r'); % red color for the first residue
    hold on;
    histogram(dist_2, xRange, 'Normalization', 'pdf', 'FaceColor', 'b'); % blue color for the second residue
end

function l = getMatrix(fileName)
    l = csvread(fileName, 0, 4);
    a = l(:, 1:2);
    b = l(:, 3:4);
    l = [b a]; % use phi, psi, dist, dist'
    l = l(:, 1:4); % use phi, psi, and chi_1
end