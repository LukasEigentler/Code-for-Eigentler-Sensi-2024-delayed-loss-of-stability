clear; 
close all;
L_col = [11, 13, 15:5:30, 40];
% L_col = 15;
f = figure;

col = lines;

cd Pattern_Stability\
eckhaus_data = importeckhaus_klausmeier('Eckhaus_stab_boundary/b.eckhaus_deltac');
cd ..\Pattern_generation\
homoc_data = importhomoclinic_klausmeier('b.homoclinic');
fold_data = importhomoclinic_klausmeier('b.fold');
fold2_data = importhomoclinic_klausmeier('b.fold2');
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
plot_index_homoc = setdiff(plot_index,exclude_index);

% plot_index_homoc_remove11 = find(homoc_data(plot_index_homoc,3)>1.34);
% plot_index_homoc_remove12 = find(homoc_data(plot_index_homoc,3)<1.3868);
% plot_index_homoc_remove13 = find(homoc_data(plot_index_homoc,2)<330.5);
% plot_index_homoc_remove1 = intersect(intersect(plot_index_homoc_remove11,plot_index_homoc_remove12),plot_index_homoc_remove13);
plot_index_homoc_remove1=[];
plot_index_homoc_remove21 = find(homoc_data(plot_index_homoc,3)<0.182);
plot_index_homoc_remove22 = find(homoc_data(plot_index_homoc,2)>274.8);
plot_index_homoc_remove23 = find(homoc_data(plot_index_homoc,2)<311.4);
plot_index_homoc_remove2 = intersect(intersect(plot_index_homoc_remove21,plot_index_homoc_remove22),plot_index_homoc_remove23);
plot_index_homoc_remove31 = find(homoc_data(plot_index_homoc,3)>1.36);
plot_index_homoc_remove32 = find(homoc_data(plot_index_homoc,3)<1.3868);
plot_index_homoc_remove33 = find(homoc_data(plot_index_homoc,2)>330.5);
plot_index_homoc_remove3 = intersect(intersect(plot_index_homoc_remove31,plot_index_homoc_remove32),plot_index_homoc_remove33);
plot_index_homoc_remove = union(union(plot_index_homoc_remove1,plot_index_homoc_remove2),plot_index_homoc_remove3);
plot_index_homoc(plot_index_homoc_remove) = [];

p(2) = plot(homoc_data(plot_index_homoc,2),homoc_data(plot_index_homoc,3), '.', 'displayname', 'Homoclinic', 'color', col(1,:));

exclude_index = find(fold_data(:,1)~=-3);
plot_index = 1:length(fold_data(:,1));
plot_index_fold = setdiff(plot_index,exclude_index);
p(3) = plot(fold_data(plot_index_fold,3),fold_data(plot_index_fold,2), '.', 'displayname', 'Fold', 'color', col(1,:));

exclude_index = find(fold2_data(:,1)~=-4);
plot_index = 1:length(fold2_data(:,1));
plot_index_fold2 = setdiff(plot_index,exclude_index);
p(3) = plot(fold2_data(plot_index_fold2,3),fold2_data(plot_index_fold2,2), '.', 'displayname', 'Fold', 'color', col(1,:));

exclude_index = find(hopf_data(:,1)~=2);
plot_index = 1:length(hopf_data(:,1));
plot_index_hopf = setdiff(plot_index,exclude_index);
plot_index_hopf_remove11 = find(hopf_data(plot_index_hopf,2)<1.34);
plot_index_hopf_remove12 = find(hopf_data(plot_index_hopf,2)>0.901);
plot_index_hopf_remove1 = intersect(plot_index_hopf_remove11,plot_index_hopf_remove12);
plot_index_hopf_remove21 = find(hopf_data(plot_index_hopf,2)<0.12);
plot_index_hopf_remove22 = find(hopf_data(plot_index_hopf,3)>311.4);
plot_index_hopf_remove2 = intersect(plot_index_hopf_remove21,plot_index_hopf_remove22);
plot_index_hopf_remove = union(plot_index_hopf_remove1,plot_index_hopf_remove2);
% plot_index_hopf_remove = zeros(1,length(plot_index_hopf));
% for ll = 1:length(plot_index_hopf)
%     [mindist, mindistind] = min(abs(hopf_data(plot_index_hopf(ll),2) - fold_data(plot_index_fold,2)));
%     if mindist <1e-2 && hopf_data(plot_index_hopf(ll),3) - fold_data(plot_index_fold(mindistind),3) < 0
%         plot_index_hopf_remove(ll) = 1;
%     end
% end
plot_index_hopf(plot_index_hopf_remove) = [];
p(4) = plot(hopf_data(plot_index_hopf,3),hopf_data(plot_index_hopf,2), '.', 'displayname', 'Hopf', 'color', col(1,:));

% xlabel('Max. mussel growth rate, $\delta$', 'interpreter','latex')
xlabel('Bifurcation parameter, $\delta$', 'interpreter','latex')
% set(gca, "XTick", [0,243.98,334.74])
% set(gca, "XTickLabel", ["0","$A_b^L$","$A_b^U$"], "TickLabelInterpreter","latex")
ylabel('Migration speed, $c$', 'interpreter','latex')
xlim([230,350])
ylim([0,1.4])
pbaspect([1 1 1])

cd Pattern_generation\
for ll=1:length(L_col)
    contour_data = importwavelength_klausmeier(['Wavelength_cont/b.wavelength_contour_',num2str(L_col(ll))]);
    exclude_index = find(contour_data(:,1)~=-2);
    fold_ind0 = find(contour_data(:,2)==5);
    [~,amaxind] = max(contour_data(fold_ind0,3));
    fold_ind = fold_ind0(amaxind);
    if length(fold_ind0)>1
        exclude_index = [exclude_index', fold_ind:length(contour_data(:,1))];
    end
    plot_index = 1:length(contour_data(:,1));
    plot_index = setdiff(plot_index,exclude_index);
    plot(contour_data(plot_index,3),contour_data(plot_index,4), '.', 'color', 'k', 'Markersize', 0.5);
    % textind1 = find(abs(contour_data(plot_index,3)-275) <1e-1);
    textind1 = find(contour_data(plot_index,4)>0.35);
    % textind3 = intersect(textind1,textind2);
    [~,textind] = min(abs(contour_data(plot_index(textind1),3)-300));
    text(contour_data(plot_index(textind1(textind)),3),contour_data(plot_index(textind1(textind)),4)-0.02,num2str(L_col(ll)))
    if L_col(ll) == 15
        [~,markerind] = min(abs(contour_data(plot_index,3)-280));
         plot(contour_data(plot_index(markerind),3),contour_data(plot_index(markerind),4), 'o', 'color', 'g', 'Markersize', 2);
         % [~,markerind] = min(abs(contour_data(plot_index,3)-285));
         % plot(contour_data(plot_index(markerind),3),contour_data(plot_index(markerind),4), 'o', 'color', 'g', 'Markersize', 2);
    end
end

% legend(p, 'location', 'northwest')

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[0 0 8 8])

% exportgraphics(f,"../../../Figures/busse_mussels.eps", "Resolution", 1000 )