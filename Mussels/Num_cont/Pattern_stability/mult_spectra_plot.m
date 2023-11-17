%% plot multiple spectra at once for the mussels model

clear;
close all;

% L_col = [20,25,30,35,40];
L_col = [15,17.4999,20,24.9999,30,40]; delta = 277.7309; % first transition in numsim
% L_col = [20,21,22,23,24]; delta = 257.4; % transition to extinction in numsim
f = figure(1);
f1 = figure(2);
for ll = 1:length(L_col)
 L = L_col(ll);
deltastr = strrep(num2str(delta),'.','dot'); Lstr = strrep(num2str(L),'.','dot');
datafilename = ['spectrum_data/mussels_delta_', deltastr,'_and_L_',Lstr];
data = load(datafilename);

figure(f)
subplot(length(L_col),1,ll)
plot(data.spectrum_real, data.spectrum_imag, '.')
title("L = " + num2str(L, '%.0f'))
pbaspect([6 1 1])
grid on
xlim([-2,0.1])

figure(f1)
hold on
plot3(data.spectrum_real, data.spectrum_imag, L*ones(1,length(data.spectrum_imag)), '.')
grid on
end

figure(f)
xlabel('$\Re(\lambda)$', 'interpreter', 'latex')
ylabel('$\Im(\lambda)$', 'interpreter', 'latex')
sgtitle("$\delta = " + num2str(delta,'%.0f')+ "$", 'interpreter','latex')
set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[5.697361111111110 3.033888888888889 10 17])
% exportgraphics(f,"../../../Figures/mussels_spectra_col.eps", "Resolution", 500 )

figure(f1)
xlabel('$\Re(\lambda)$', 'interpreter', 'latex')
ylabel('$\Im(\lambda)$', 'interpreter', 'latex')
zlabel('$L$', 'interpreter', 'latex')
xlim([-2,0.1])
sgtitle("$\delta = " + num2str(delta,'%.0f')+"$", 'interpreter', 'latex')
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[5.697361111111110 3.033888888888889 10 10])
%exportgraphics(f1,"../../../Figures/mussels_spectra_col_3d.eps", "Resolution", 500 )