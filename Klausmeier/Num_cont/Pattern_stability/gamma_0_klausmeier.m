%% Stability calculations for /gamma = 0: Matrix eigenvalue problem
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023

clear; close all

no_species = input(['Select a number to choose from the following options: \n', ...
    '-2: trace Hopf stability boundary; \n' ...
    '-1: trace Eckhaus stability boundary; \n' ...
    '1: compute spectrum for a given pattern; \n'])

%% Import data

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

%% interpolation to create solution on mesh with equidistant points
z = transpose(linspace(0,1,N+1)); 
sol = zeros(N+1,5);
sol(:,1) = L*z; % rescale z
sol(:,2:end) = interp1(sol1(:,1),sol1(:,2:end),z);
sol = sol(1:end-1,:); % remove z=L data
z = sol(:,1);
dx = z(2) - z(1);
u1 = sol(:,2); w = sol(:,4);

%% plot solution
f = figure;
subplot(3,1,1)
plot(z,u1)
grid on
xlabel("TW coordinate, $z$", "Interpreter","latex")
ylabel("Plant density, $u(z)$", "Interpreter","latex")

subplot(3,1,2)
plot(z,w)
grid on
xlabel("TW coordinate, $z$", "Interpreter","latex")
ylabel("Water density, $u(z)$", "Interpreter","latex")



%% Derivatives in Jacobian
dfdu1 = w.*2.*u1-B1;
dfdw = u1.*(u1);
dhdu1 = -w.*(u1)-w.*(u1);

dhdw = -1 - (u1).*(u1);

%% Matrix two species model

M = sparse(2*N, 2*N);

% rows affected by BCs
M(1,1) = dfdu1(1) - 2/(dx^2) - c/dx;
M(N+1,N+1) = dhdw(1) - d*2/(dx^2) - (c+nu)/dx;

M(1,2) = 1/(dx^2) + c/dx;
M(N+1,N+2) = d/(dx^2) + (c+nu)/dx;

M(1,N) = 1/(dx^2) ;
M(N+1,2*N) = d/(dx^2);

M(1,N+1) = dfdw(1); 
M(N+1,1) = dhdu1(1); 


M(N,N) = dfdu1(N) - 2/(dx^2) -c/dx;
M(2*N,2*N) = dhdw(N) - d*2/(dx^2)  -(c+nu)/dx;

M(N,N-1) = 1/(dx^2);
M(2*N,2*N-1) = d/(dx^2);

M(N,1) = 1/(dx^2) + c/dx;
M(2*N,N+1) = d/(dx^2) + (c+nu)/dx;

M(N,2*N) = dfdw(N);
M(2*N,N) = dhdu1(N); 
% rows unaffected by BCs

for ii = 2:N-1
    M(ii,ii-1) = 1/(dx^2);
    M(N+ii,N+ii-1) = d/(dx^2);
    
    M(ii,ii) = dfdu1(ii) - 2/(dx^2)- c/dx;
    M(N+ii,N+ii) = dhdw(ii) - d*2/(dx^2)- (c+nu)/dx;
    
    M(ii,ii+1) = 1/(dx^2) + c/dx;
    M(N+ii,N+ii+1) = d/(dx^2) + (c+nu)/dx;
    
    M(ii,N+ii) = dfdw(ii);
    M(N+ii,ii) = dhdu1(ii);
end



%% Eigenvalues
if no_species == 1
    no_eigen = 21;
elseif no_species == -1
    no_eigen = 2;
end
% two species model
[wM, lambdaM] = eigs(M,2*N);
lambdaM = diag(lambdaM);

% split eigenvalues in real and imag parts
lambdaM_real = real(lambdaM);
lambdaM_imag = imag(lambdaM);
[lambdaM_real_sorted,lambdaM_real_sorted_index]  = sort(lambdaM_real); % sort eigenvalues by real part

% plot spectrum
figure
 plot(lambdaM_real,lambdaM_imag,'o')
hold on
xlim([lambdaM_real_sorted(end-no_eigen+1),lambdaM_real_sorted(end)])
grid on
xlabel('$\Re(\lambda)$', 'interpreter', 'latex')
ylabel('$\Im(\lambda)$', 'interpreter', 'latex')

L_A = [L,A, c];
if no_species > 0
    save stab_data/period_rain.dat L_A -ascii
elseif no_species == 0
     save stab_data/period_rain.dat L_A -ascii
elseif no_species == -1
     save Eckhaus_stab_boundary/period_rain.dat L_A -ascii
elseif no_species == -2
     save Hopf_stab_boundary/period_rain.dat L_A -ascii
end

%% Determine eigenvalues with largest real part
largest_index = lambdaM_real_sorted_index(end-no_eigen+1:end);
lambda_start_col = [lambdaM_real(largest_index), lambdaM_imag(largest_index)];
if no_species >0
    save stab_data/eigenvalues_col.dat lambda_start_col -ascii
end


%% Extract data of one eigenvalue as starting value for continuation in gamma (two species model only)
if no_species > 0
for ff = 1:no_eigen

start_index = largest_index(ff);
lambda_real_start = lambdaM_real(start_index); % split into real and imag parts
lambda_imag_start = lambdaM_imag(start_index);
lambda_start = [lambdaM_real(start_index), lambdaM_imag(start_index)];
% parameters_start = [A, B1, b2, c, chi, d, D2, F2, H2, nu];

% derivatives of eigenvector entries
pp = 1:N-1;
w_u1p_start(pp) = (wM(pp+1,start_index) - wM(pp,start_index))/dx;
w_u1p_start(N) = (wM(1,start_index) - wM(N,start_index))/dx;
w_wp_start(pp) = (wM(N+pp+1,start_index) - wM(N+pp,start_index))/dx;
w_wp_start(N) = (wM(N+1,start_index) - wM(2*N,start_index))/dx;



% combine eigenvector data and normalise
w_start = [ wM(1:N,start_index), transpose(w_u1p_start),wM(N+1:2*N,start_index), transpose(w_wp_start)]; 
w_start_real = real(w_start);
w_start_imag = imag(w_start);
toint = 0;
for nn = 1:length(w_start_real(1,:))
    toint = toint + w_start_real(:,nn).^2 + w_start_imag(:,nn).^2;
end
z_norm = z/L;
w_norm = trapz(z_norm,toint);
w_start_real = w_start_real/sqrt(w_norm);
w_start_imag = w_start_imag/sqrt(w_norm);

clear gamma_0_sol_dat
% combine with PTW solution and write to file
gamma_0_sol_dat = [sol(1:N,:), w_start_real, w_start_imag];
gamma_0_sol_dat = [gamma_0_sol_dat; gamma_0_sol_dat(1,:)]; % add z=L data using periodicity
gamma_0_sol_dat(end,1) = L;
gamma_0_sol_dat(:,1) = gamma_0_sol_dat(:,1)/L;

figfilename = ['stab_data/gamma_0_sol', num2str(ff),'.dat'];
save(figfilename, 'gamma_0_sol_dat','-ascii')


end

elseif no_species <= 0
    if no_species == 0
        ff = input('enter index of eigenvalue to be continued to origin \n');
    elseif no_species == -1
        ff = input('enter index of eigenvalue to be checked for Eckhaus stability change \n');
    elseif no_species == -2
        ff = input('enter index of eigenvalue from which Hopf stability change boundary calculation should be started \n');
    end
    start_index = largest_index(ff);
    lambda_real_start = lambdaM_real(start_index); % split into real and imag parts
    lambda_imag_start = lambdaM_imag(start_index);
    lambda_start = [lambdaM_real(start_index), lambdaM_imag(start_index)];
    if no_species == 0
        save stab_data/eigenvalues.dat lambda_start -ascii
    elseif no_species == -1
         save Eckhaus_stab_boundary/eigenvalues.dat lambda_start -ascii
    elseif no_species == -2
         save Hopf_stab_boundary/eigenvalues.dat lambda_start -ascii
    end
        % parameters_start = [A, B1, b2, c, chi, d, D2, F2, H2, nu];

    % derivatives of eigenvector entries
    pp = 1:N-1;
    w_u1p_start(pp) = (wM(pp+1,start_index) - wM(pp,start_index))/dx;
    w_u1p_start(N) = (wM(1,start_index) - wM(N,start_index))/dx;
    
    w_wp_start(pp) = (wM(N+pp+1,start_index) - wM(N+pp,start_index))/dx;
    w_wp_start(N) = (wM(N+1,start_index) - wM(2*N,start_index))/dx;



    % combine eigenvector data and normalise
    w_start = [ wM(1:N,start_index), transpose(w_u1p_start),wM(N+1:2*N,start_index), transpose(w_wp_start)]; 
    w_start_real = real(w_start);
    w_start_imag = imag(w_start);
    toint = 0;
    for nn = 1:length(w_start_real(1,:))
        toint = toint + w_start_real(:,nn).^2 + w_start_imag(:,nn).^2;
    end
    z_norm = z/L;
    w_norm = trapz(z_norm,toint);
    w_start_real = w_start_real/sqrt(w_norm);
    w_start_imag = w_start_imag/sqrt(w_norm);

    clear gamma_0_sol_dat
    % combine with PTW solution and write to file
    gamma_0_sol_dat = [sol(1:N,:), w_start_real, w_start_imag];
    gamma_0_sol_dat = [gamma_0_sol_dat; gamma_0_sol_dat(1,:)]; % add z=L data using periodicity
    gamma_0_sol_dat(end,1) = L;
    gamma_0_sol_dat(:,1) = gamma_0_sol_dat(:,1)/L;
    
    if no_species == 0
        figfilename = 'stab_data/gamma_0_sol.dat';
    elseif no_species == -1
        figfilename = 'Eckhaus_stab_boundary/gamma_0_sol.dat';
           gamma_0_sol_dat = [gamma_0_sol_dat, ones(N+1,12)];
    elseif no_species == -2
        figfilename = 'Hopf_stab_boundary/gamma_0_sol.dat';

    end  
    save(figfilename, 'gamma_0_sol_dat','-ascii')

end

fprintf('Press any key to continue after continuation has been done /n')
pause

%%
figure(f);
subplot(3,1,3)
hold on
grid on
xlabel('$\Re(\lambda)$', 'interpreter', 'latex')
ylabel('$\Im(\lambda)$', 'interpreter', 'latex')
if no_species > 0
spectrum = importspectrumdata_klausmeier('b.spectrum_full');
spectrum = double(spectrum);
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
        % xplot = [xplot, spectrum(discont_ind(pp)+1:discont_ind(pp+1),1)];
        % yplot = [yplot, spectrum(discont_ind(pp)+1:discont_ind(pp+1),2)];
        % plot(spectrum(discont_ind(pp)+1:discont_ind(pp+1),1),spectrum(discont_ind(pp)+1:discont_ind(pp+1),2), '.', 'linewidth', 2, 'color', col(2,:));
        spectrum_real = [spectrum_real; spectrum(discont_ind(pp)+1:discont_ind(pp+1),1)];
        spectrum_imag = [spectrum_imag; spectrum(discont_ind(pp)+1:discont_ind(pp+1),2)];
    end

    plot(spectrum_real,spectrum_imag, '.', 'linewidth', 2, 'color', col(2,:));
xlim([lambdaM_real_sorted(end-floor(no_eigen/2)),max(spectrum(:,1))+1])
ylim([-0.4,0.4])
[p,z] = zoomPlot(spectrum_real,spectrum_imag,[-0.05 0.05],[0.72 0.16 0.16 0.16],[]); 


Astr = num2str(A);
Astr = strrep(Astr,'.','dot');
cstr = num2str(L);
cstr = strrep(cstr,'.','dot');


    % title(['$A = ', num2str(A), '$, $L= ', num2str(L), '$'], 'interpreter', 'latex')
    figfilename = ['plots/Klausmeier_A_', num2str(Astr),'_and_L_',num2str(cstr)];
    datafilename = ['spectrum_data/Klausmeier_A_', num2str(Astr),'_and_L_',num2str(cstr)];

% pbaspect([3 1 1])

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[5.697361111111110 3.033888888888889 12 12])
% exportgraphics(f,"../../../../Mattia_Sensi_Lukas_Eigentler_shared_space/klausmeier_sol_spectra_example.eps", "Resolution", 1000,'ContentType','vector' )


savefig(figfilename)
saveas(gcf,figfilename,'jpg')
saveas(gcf,figfilename,'epsc')
save(datafilename, "spectrum_real", "spectrum_imag", "A", "c", "L", "B1", "nu", "d")

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


%%

