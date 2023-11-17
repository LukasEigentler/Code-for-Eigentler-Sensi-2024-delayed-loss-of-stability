%% Spectrum post processing and save

clear; close all;
no_species = 1;
cont_data = importcontinuationdata_mussels('../Pattern_generation/output/b.ptw');
ptw_sol_data = importsolutiondata_mussels('../Pattern_generation/output/s.ptw');


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
delta = cont_data(sol_output_ind(sol_no)+cont_data_start_index-1,8);
L = cont_data(sol_output_ind(sol_no)+cont_data_start_index-1,7);
c = cont_data(sol_output_ind(sol_no)+cont_data_start_index-1,9);
parameters = load('../parameters.dat');
alpha = parameters(1); beta = parameters(2); d = parameters(6); nu = parameters(7); eta = parameters(4); theta = parameters(5);

f = figure;
hold on
grid on
xlabel('$\Re(\lambda)$', 'interpreter', 'latex')
ylabel('$\Im(\lambda)$', 'interpreter', 'latex')

spectrum = importspectrumdata_klausmeier('b.spectrum_full');
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
    spectrum_real = []; spectrum_imag = [];
    for pp = 1:length(discont_ind)-1     
        plot(spectrum(discont_ind(pp)+1:discont_ind(pp+1),1),spectrum(discont_ind(pp)+1:discont_ind(pp+1),2), '.', 'linewidth', 2, 'color', col(1,:));
        spectrum_real = [spectrum_real; spectrum(discont_ind(pp)+1:discont_ind(pp+1),1)];
        spectrum_imag = [spectrum_imag; spectrum(discont_ind(pp)+1:discont_ind(pp+1),2)];
    end

xlim([min(spectrum(:,1)),max(spectrum(:,1))])
deltastr = num2str(delta);
deltastr = strrep(deltastr,'.','dot');
cstr = num2str(L);
cstr = strrep(cstr,'.','dot');


    % title(['$A = ', num2str(A), '$, $L= ', num2str(L), '$'], 'interpreter', 'latex')
    figfilename = ['plots/mussels_delta_', num2str(deltastr),'_and_L_',num2str(cstr)];
    datafilename = ['spectrum_data/mussels_delta_', num2str(deltastr),'_and_L_',num2str(cstr)];

pbaspect([3 1 1])

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[5.697361111111110 3.033888888888889 15 4])
cd ../spectra_calc_batch/
savefig(figfilename)
saveas(gcf,figfilename,'jpg')
saveas(gcf,figfilename,'epsc')
save(datafilename, "spectrum_real", "spectrum_imag", "delta", "c", "L", "alpha", "nu", "d","theta","eta","beta")

