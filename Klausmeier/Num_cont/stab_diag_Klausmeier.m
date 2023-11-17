%% Stability diagram for the Klausmeier model
% This script loads numerical continuation data and visualises the Busse
% balloon.
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023

clear; close all;
% L_col = [20, 30, 40:20:100, 150, 200:100:400];
L_col = [20,40,60,100,200,400];
% L_col = 40;
f = figure;

col = lines;

cd Pattern_Stability\
eckhaus_data = importeckhaus_klausmeier('Eckhaus_stab_boundary/b.eckhaus_Ac');
cd ..\Pattern_generation\
homoc_data = importhomoclinic_klausmeier('b.homoclinic');
fold_data = importhomoclinic_klausmeier('b.fold');
hopf_data = importhomoclinic_klausmeier('b.hopf');
cd ../



%pattern_onset_dat1 = load('output/Multi_species_pattern_onset1')
%load('output/Multi_species_pattern_onset')

hold on 
grid on
exclude_index = find(eckhaus_data(:,1)==1);
plot_index = 1:length(eckhaus_data(:,1));
plot_index = setdiff(plot_index,exclude_index);
p(1) = plot(eckhaus_data(plot_index,1),eckhaus_data(plot_index,2), '.', 'displayname', 'Eckhaus stability boundary', 'color', col(2,:));

exclude_index = find(homoc_data(:,1)~=-2);
plot_index = 1:length(homoc_data(:,1));
plot_index = setdiff(plot_index,exclude_index);
p(2) = plot(homoc_data(plot_index,2),homoc_data(plot_index,3), '.', 'displayname', 'Homoclinic', 'color', col(1,:));

exclude_index = find(fold_data(:,1)~=-3);
plot_index = 1:length(fold_data(:,1));
plot_index_fold = setdiff(plot_index,exclude_index);
p(3) = plot(fold_data(plot_index_fold,3),fold_data(plot_index_fold,2), '.', 'displayname', 'Fold', 'color', col(1,:));

exclude_index = find(hopf_data(:,1)~=2);
plot_index = 1:length(hopf_data(:,1));
plot_index_hopf = setdiff(plot_index,exclude_index);
plot_index_hopf_remove = zeros(1,length(plot_index_hopf));
for ll = 1:length(plot_index_hopf)
    [mindist, mindistind] = min(abs(hopf_data(plot_index_hopf(ll),2) - fold_data(plot_index_fold,2)));
    if mindist <1e-2 && hopf_data(plot_index_hopf(ll),3) - fold_data(plot_index_fold(mindistind),3) < 0
        plot_index_hopf_remove(ll) = 1;
    end
end
plot_index_hopf(plot_index_hopf_remove==1) = [];
p(4) = plot(hopf_data(plot_index_hopf,3),hopf_data(plot_index_hopf,2), '.', 'displayname', 'Hopf', 'color', col(1,:));

% xlabel('Rainfall, $A$', 'interpreter','latex')
xlabel('Bifurcation parameter, $A$', 'interpreter','latex')
% set(gca, "XTick", [0,0.1798,3.45567])
% set(gca, "XTickLabel", ["0","$A_b^L$","$A_b^U$"], "TickLabelInterpreter","latex")
ylabel('Migration speed, $c$', 'interpreter','latex')
xlim([0,3.5])
ylim([0,4.1])
pbaspect([1 1 1])

cd Pattern_generation\
for ll=1:length(L_col)
    contour_data = importwavelength_klausmeier(['Wavelength_cont/b.wavelength_contour_',num2str(L_col(ll))]);
    exclude_index = find(contour_data(:,1)~=-2);
    fold_ind0 = find(contour_data(:,2)==5);
    [~,amaxind] = max(contour_data(fold_ind0,3));
    fold_ind = fold_ind0(amaxind);
    if ~isempty(fold_ind)
        exclude_index = [exclude_index', fold_ind:length(contour_data(:,1))];
    end
    plot_index = 1:length(contour_data(:,1));
    plot_index = setdiff(plot_index,exclude_index);
    plot(contour_data(plot_index,3),contour_data(plot_index,4), '.', 'color', 'k', 'Markersize', 0.5);
    textind1 = find(abs(contour_data(plot_index,3)-1.85) <1e-2);
    [~,textind] = min(contour_data(plot_index(textind1),4));
    text(contour_data(plot_index(textind1(textind)),3),contour_data(plot_index(textind1(textind)),4)-0.06,num2str(L_col(ll)))
    if L_col(ll) == 20 % marker on specific contour
        [~,markerind] = min(abs(contour_data(plot_index,3)-1.4));
         plot(contour_data(plot_index(markerind),3),contour_data(plot_index(markerind),4), 'o', 'color', 'g', 'Markersize', 2);
         % [~,markerind] = min(abs(contour_data(plot_index,3)-1.6));
         % plot(contour_data(plot_index(markerind),3),contour_data(plot_index(markerind),4), 'o', 'color', 'g', 'Markersize', 2);
    end
end

% legend(p, 'location', 'northwest')

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[0 0 8 8])

