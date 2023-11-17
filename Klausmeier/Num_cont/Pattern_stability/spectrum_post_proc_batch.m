%% Spectrum post processing and save
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023

clear; close all;
no_species = 1;
cont_data = importcontinuationdata_klausmeier('../Pattern_generation/output/b.ptw');
ptw_sol_data = importsolutiondata_klausmeier('../Pattern_generation/output/s.ptw');


N = find(ptw_sol_data(:,1) == 1); % determine no of meshpoints used by AUTO
N = N(1) - 1; 

% find and extract the user-specified ('UZ1') output from AUTO data
cont_data_start_index = find(cont_data(:,1) ~=0);
cont_data_start_index = cont_data_start_index(1);
sol_output_ind = find(cont_data(cont_data_start_index:end,3) ~= 0); 
sol_output_user_defined_ind = find(cont_data(cont_data_start_index:end,3) == -4); %set -4 for user defined output; 6 for BP, -2 for other endpoint of cont
sol_no = find(sol_output_ind == sol_output_user_defined_ind(1));
sol1 = ptw_sol_data((sol_no-1)*(N +4)+1:(sol_no-1)*(N +4)+1+N,:);


%% Parameters
A = cont_data(sol_output_ind(sol_no)+cont_data_start_index-1,8);
L = cont_data(sol_output_ind(sol_no)+cont_data_start_index-1,7);
c = cont_data(sol_output_ind(sol_no)+cont_data_start_index-1,9);
parameters = load('../parameters.dat');
B1 = parameters(1); d = parameters(3); nu = parameters(4); 


f = figure;
hold on
grid on
xlabel('$\Re(\lambda)$', 'interpreter', 'latex')
ylabel('$\Im(\lambda)$', 'interpreter', 'latex')
if no_species > 0
spectrum = importspectrumdata_klausmeier('b.spectrum_full'); spectrum = double(spectrum);
hold on
col = lines;

    diff_data1 = abs(diff(spectrum(:,1)));
    diff_data2 = abs(diff(spectrum(:,2)));
    discont_ind = [0];
    for ii = 2:length(diff_data2)-1
       if (diff_data2(ii) > 3*diff_data2(ii-1) &&  diff_data2(ii) > 3*diff_data1(ii+1)) || (diff_data1(ii) > 3*diff_data1(ii-1) &&  diff_data1(ii) > 3*diff_data1(ii+1))
           discont_ind = [discont_ind, ii];
       end
    end
    discont_ind = [discont_ind, length(spectrum(:,2))];
    spectrum_real = []; spectrum_imag = []; gamma = [];
    for pp = 1:length(discont_ind)-1     
        plot(spectrum(discont_ind(pp)+1:discont_ind(pp+1),1),spectrum(discont_ind(pp)+1:discont_ind(pp+1),2), '.', 'linewidth', 2, 'color', col(1,:));
        spectrum_real = [spectrum_real; spectrum(discont_ind(pp)+1:discont_ind(pp+1),1)];
        spectrum_imag = [spectrum_imag; spectrum(discont_ind(pp)+1:discont_ind(pp+1),2)];
        gamma = [gamma; spectrum(discont_ind(pp)+1:discont_ind(pp+1),3)];
    end

xlim([min(spectrum(:,1)),max(spectrum(:,1))])
Astr = num2str(A);
Astr = strrep(Astr,'.','dot');
cstr = num2str(L);
cstr = strrep(cstr,'.','dot');


    % title(['$A = ', num2str(A), '$, $L= ', num2str(L), '$'], 'interpreter', 'latex')
    figfilename = ['plots/Klausmeier_A_', num2str(Astr),'_and_L_',num2str(cstr)];
    datafilename = ['spectrum_data/Klausmeier_A_', num2str(Astr),'_and_L_',num2str(cstr)];

pbaspect([3 1 1])

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[5.697361111111110 3.033888888888889 15 4])
cd ../spectra_calc_batch/
savefig(figfilename)
saveas(gcf,figfilename,'jpg')
saveas(gcf,figfilename,'epsc')
save(datafilename, "gamma", "spectrum_real", "spectrum_imag", "A", "c", "L", "B1", "nu", "d")

elseif no_species == -1
    col = lines;
    eckhaus_data = importeckhaus_klausmeier('Eckhaus_stab_boundary/b.eckhaus_Ac');
    %pattern_onset_dat1 = load('output/Multi_species_pattern_onset1')
    %load('output/Multi_species_pattern_onset')
    figure
    hold on 
    grid on
    exclude_index = find(eckhaus_data(:,1)==1);
    plot_index = 1:length(eckhaus_data(:,1));
    plot_index = setdiff(plot_index,exclude_index);
    p(1) = plot(eckhaus_data(plot_index,1),eckhaus_data(plot_index,2), '.', 'displayname', 'Eckhaus stability boundary');
    %p(2) = plot(pattern_onset_dat(:,1),pattern_onset_dat(:,2),'.', 'color',col(2,:), 'displayname', 'Pattern onset boundary');
   % plot(pattern_onset_dat1.pattern_onset_dat(:,1),pattern_onset_dat1.pattern_onset_dat(:,2),'.', 'color',col(2,:), 'displayname', 'Pattern onset boundary');


    legend(p, 'location', 'northwest')
    xlabel('$A$', 'interpreter','latex')
    ylabel('$c$', 'interpreter','latex')
    xlim([0,4])
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    pbaspect([1 1 1])

end