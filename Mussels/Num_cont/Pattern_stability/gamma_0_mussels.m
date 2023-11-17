%% Stability calculations for /gamma = 0: Matrix eigenvalue problem

clear; close all

no_species = input(['Select a number to choose from the following options: \n', ...
    '-2: trace Hopf stability boundary; \n' ...
    '-1: trace Eckhaus stability boundary; \n' ...
    '1: compute spectrum for a given pattern; \n'])

%% Import data

cont_data = importcontinuationdata_klausmeier('../Pattern_generation/output/b.ptw');
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
alpha = parameters(1); beta = parameters(2); eta = parameters(3); theta = parameters(4); D = parameters(6); nu = parameters(7);

%% interpolation to create solution on mesh with equidistant points
z = transpose(linspace(0,1,N+1)); 
sol = zeros(N+1,6);
sol(:,1) = L*z; % rescale z
sol(:,2:end) = interp1(sol1(:,1),sol1(:,2:end),z);
sol = sol(1:end-1,:); % remove z=L data
z = sol(:,1);
dx = z(2) - z(1);
m = sol(:,2); a = sol(:,4); s = sol(:,5);

%% plot solution
f = figure;
subplot(2,2,1)
plot(z,m)
grid on
xlabel("TW coordinate, $z$", "Interpreter","latex")
ylabel("Mussel dens., $m(z)$", "Interpreter","latex")

subplot(2,2,2)
plot(z,s)
grid on
xlabel("TW coordinate, $z$", "Interpreter","latex")
ylabel("Sediment dens., $s(z)$", "Interpreter","latex")

subplot(2,2,3)
plot(z,a)
grid on
xlabel("TW coordinate, $z$", "Interpreter","latex")
ylabel("Algae dens., $a(z)$", "Interpreter","latex")

%% Derivatives in Jacobian
dfdm = (a.*delta.*(eta + s))./(s + 1) - 1;
dfda = (delta*m.*(eta + s))./(s + 1);
dfds = (a.*delta.*m)./(s + 1) - (a.*delta.*m.*(eta + s))./(s + 1).^2;
dgdm = -(a.*beta.*(eta + s))./(s + 1);
dgda = - alpha - (beta*m.*(eta + s))./(s + 1);
dgds = (a.*beta.*m.*(eta + s))./(s + 1).^2 - (a.*beta.*m)./(s + 1);
dhdm = 1.0*ones(1,length(m));
dhda = 0.0*ones(1,length(m));
dhds = -theta*ones(1,length(m));

%% Matrix two species model

M = sparse(3*N, 3*N);

% rows affected by BCs
M(1,1) = dfdm(1) - 2/(dx^2) - c/dx;
M(1,2) = 1/(dx^2) + c/dx;
M(1,N) = 1/(dx^2);
M(1,N+1) = dfda(1); 
M(1,2*N+1) = dfds(1);

M(N+1,N+1) = dgda(1) - (c+nu)/dx;
M(N+1,N+2) = (c+nu)/dx;
M(N+1,1) = dgdm(1);
M(N+1,2*N+1) = dgds(1);


M(2*N+1,2*N+1) = dhds(1) - 2*D/(dx^2) - c/dx;
M(2*N+1,2*N+2) = D/(dx^2) + c/dx;
M(2*N+1,3*N) = D/(dx^2);
M(2*N+1,1) = dhdm(1); 
M(2*N+1,N+1) = dhda(1);



M(N,N) = dfdm(N) - 2/(dx^2) -c/dx;
M(N,N-1) = 1/(dx^2);
M(N,1) = 1/(dx^2) + c/dx;
M(N,2*N) = dfda(N); 
M(N,3*N) = dfds(N);


M(2*N,2*N) = dgda(N)  -(c+nu)/dx;
M(2*N,N+1) = (c+nu)/dx;
M(2*N,N) = dgdm(N); 
M(2*N,3*N) = dgds(N);


M(3*N,3*N) = dhds(N) - 2*D/(dx^2)  -c/dx;
M(3*N,3*N-1) = D/(dx^2);
M(3*N,2*N+1) = D/(dx^2) + c/dx;
M(3*N,N) = dhdm(N);
M(3*N,2*N) = dhda(N);


% rows unaffected by BCs

for ii = 2:N-1

    M(ii,ii) = dfdm(ii) - 2/(dx^2)- c/dx;
    M(ii,ii-1) = 1/(dx^2);
    M(ii,ii+1) = 1/(dx^2) + c/dx;
    M(ii,N+ii) = dfda(ii);
    M(ii,2*N+ii) = dfds(ii);

    M(N+ii,N+ii) = dgda(ii) - (c+nu)/dx;
    M(N+ii,N+ii+1) = (c+nu)/dx;
    M(N+ii,ii) = dgdm(ii);
    M(N+ii,2*N+ii) = dgds(ii);
    
    M(2*N+ii,2*N+ii) = dhds(ii) - 2*D/(dx^2) - c/dx;
    M(2*N+ii,2*N+ii-1) = D/(dx^2);
    M(2*N+ii,2*N+ii+1) = D/(dx^2) + c/dx;
    M(2*N+ii,ii) = dhdm(ii);
    M(2*N+ii,N+ii) = dhda(ii);
   
end




%% Eigenvalues
if no_species == 1
    no_eigen = 21;
elseif no_species == -1
    no_eigen = 3;
end
% two species model
[wM, lambdaM] = eigs(M,3*N);
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

L_delta = [L,delta, c];
if no_species > 0
    save stab_data/period_delta.dat L_delta -ascii
elseif no_species == 0
     save stab_data/period_delta.dat L_delta -ascii
elseif no_species == -1
     save Eckhaus_stab_boundary/period_delta.dat L_delta -ascii
elseif no_species == -2
     save Hopf_stab_boundary/period_delta.dat L_delta -ascii
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
w_mp_start(pp) = (wM(pp+1,start_index) - wM(pp,start_index))/dx;
w_mp_start(N) = (wM(1,start_index) - wM(N,start_index))/dx;
w_sp_start(pp) = (wM(2*N+pp+1,start_index) - wM(2*N+pp,start_index))/dx;
w_sp_start(N) = (wM(2*N+1,start_index) - wM(3*N,start_index))/dx;



% combine eigenvector data and normalise
w_start = [ wM(1:N,start_index), transpose(w_mp_start),wM(N+1:2*N,start_index), wM(2*N+1:3*N,start_index), transpose(w_sp_start)]; 
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
    w_mp_start(pp) = (wM(pp+1,start_index) - wM(pp,start_index))/dx;
    w_mp_start(N) = (wM(1,start_index) - wM(N,start_index))/dx;
    
    w_sp_start(pp) = (wM(2*N+pp+1,start_index) - wM(2*N+pp,start_index))/dx;
    w_sp_start(N) = (wM(2*N+1,start_index) - wM(3*N,start_index))/dx;



    % combine eigenvector data and normalise
    w_start = [ wM(1:N,start_index), transpose(w_mp_start),wM(N+1:2*N,start_index), wM(2*N+1:3*N,start_index), transpose(w_sp_start)]; 
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
           gamma_0_sol_dat = [gamma_0_sol_dat, 0.1*ones(N+1,10)];
    elseif no_species == -2
        figfilename = 'Hopf_stab_boundary/gamma_0_sol.dat';

    end  
    save(figfilename, 'gamma_0_sol_dat','-ascii')

end

fprintf('Press any key to continue after continuation has been done /n')
pause

%%
figure(f);
subplot(2,2,4)
hold on
grid on
xlabel('$\Re(\lambda)$', 'interpreter', 'latex')
ylabel('$\Im(\lambda)$', 'interpreter', 'latex')
if no_species > 0
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
    for pp = 1:length(discont_ind)-1     
        plot(spectrum(discont_ind(pp)+1:discont_ind(pp+1),1),spectrum(discont_ind(pp)+1:discont_ind(pp+1),2), '.', 'linewidth', 2, 'color', col(2,:));
    end

    spectrum_real = []; spectrum_imag = [];
    for pp = 1:length(discont_ind)-1     
        plot(spectrum(discont_ind(pp)+1:discont_ind(pp+1),1),spectrum(discont_ind(pp)+1:discont_ind(pp+1),2), '.', 'linewidth', 2, 'color', col(2,:));
        spectrum_real = [spectrum_real; spectrum(discont_ind(pp)+1:discont_ind(pp+1),1)];
        spectrum_imag = [spectrum_imag; spectrum(discont_ind(pp)+1:discont_ind(pp+1),2)];
    end

%xlim([lambdaM_real_sorted(end-floor(no_eigen/2)),max(spectrum(:,1))])
xlim([-0.5,max(spectrum(:,1))])
deltastr = num2str(delta);
deltastr = strrep(deltastr,'.','dot');
cstr = num2str(c);
cstr = strrep(cstr,'.','dot');
Lstr = num2str(L);
Lstr = strrep(Lstr,'.','dot');


    % title(['$A = ', num2str(A), '$, $L= ', num2str(L), '$'], 'interpreter', 'latex')
    figfilename = ['plots/Mussels_delta_', num2str(deltastr),'_c_',num2str(cstr),'_L_',num2str(Lstr)];
    datafilename = ['spectrum_data/mussels_delta_', num2str(deltastr),'_and_L_',num2str(Lstr)];



set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[5.697361111111110 3.033888888888889 10 10])
% exportgraphics(f,"../../../Mattia_Sensi_Lukas_Eigentler_shared_space/mussels_sol_spectra_example.eps", "Resolution", 1000,'ContentType','vector' )


savefig(figfilename)
saveas(gcf,figfilename,'jpg')
saveas(gcf,figfilename,'epsc')

save(datafilename, "spectrum_real", "spectrum_imag", "alpha", "beta", "L", "eta", "delta", "theta", "nu", "D")

elseif no_species == -1
    col = lines;
    eckhaus_data = importeckhaus_klausmeier('Eckhaus_stab_boundary/b.eckhaus_deltac');
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
    xlabel('$\delta$', 'interpreter','latex')
    ylabel('$c$', 'interpreter','latex')
    xlim([240,350])
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    pbaspect([1 1 1])
    
end