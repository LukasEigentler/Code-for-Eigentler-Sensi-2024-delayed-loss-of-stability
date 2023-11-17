%% plot multiple spectra at once for the Klausmeier model
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023
clear;
close all;

L = 20;
skip = 5; % plot every skip spectrum only
%% find spectra files with correct wavelength
cd ../spectra_calc_batch/spectrum_data/
Files = dir;
correctL = NaN*ones(1,length(Files));
for ff = 1:length(Files)
    correctL(ff) = contains(Files(ff).name, ["L_"+num2str(L),"L_"+num2str(L-1)+"dot99"]);
    filenames(ff) = string(Files(ff).name);
    
end
Files = Files(correctL==1);
filenames = filenames(correctL==1);
for ff = 1:length(Files)
    filename = Files(ff).name;
    startind = strfind(filename,'_A_'); startind = startind+3;
    endind = strfind(filename,'_and_L_'); endind = endind-1;
    Aval(ff) = str2num(strrep(filename(startind:endind),'dot','.')); %determine A value
end


%% loop through spectra
f1 = figure(2);
gammanorm = linspace(0,11*2*pi,1000);
[gammamesh,Amesh] = meshgrid(gammanorm,Aval);
for ll = 1:length(Files)

    data = load(filenames(ll));
    gamma = data.gamma; gammanew = NaN*ones(1,length(gamma));
    zeroind = find(gamma==0);
    tpind(1:length(zeroind)-1) = zeroind(2:end); tpind(length(zeroind)) = length(gamma);
    for zz = 1:length(zeroind)
        gammanew(zeroind(zz):tpind(zz)) = (zz-1)*2*pi + gamma(zeroind(zz):tpind(zz));
    end
    spectrum_real_plot(ll,:) = interp1(gammanew,data.spectrum_real,gammanorm);
    spectrum_imag_plot(ll,:) = interp1(gammanew,data.spectrum_imag,gammanorm);
    A_plot(ll,:) = Aval(ll)*ones(1,length(gammanorm));
    [maxreal(ll),maxrealind] = max(spectrum_real_plot(ll,:));
    imag_at_max_real(ll) = spectrum_imag_plot(ll,maxrealind);
    
    if mod(ll-1,skip) == 0
        figure(f1)
        hold on
        plot3(data.spectrum_real, data.spectrum_imag, Aval(ll)*ones(1,length(data.spectrum_imag)), '.', "Color", [0,0,0], "MarkerSize",1)
        grid on
    end


end


figure(f1)
view(25,35)
[Aval,sortind] = sort(Aval);
plot3(maxreal(sortind),-abs(imag_at_max_real(sortind)),Aval, '-r')
xlabel('$\Re(\lambda)$', 'interpreter', 'latex')
ylabel('$\Im(\lambda)$', 'interpreter', 'latex')
zlabel('Bifurcation parameter, $A$', 'interpreter', 'latex')
xlim([-1,inf])
% sgtitle("$A = " + num2str(A,'%.2f')+"$", 'interpreter', 'latex')
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[5.697361111111110 3.033888888888889 10 10])
%exportgraphics(f1,"../../../../Mattia_Sensi_Lukas_Eigentler_shared_space/klausmeier_spectra_stack.eps", "Resolution", 1000 )

f = figure;
plot(Aval,maxreal(sortind), 'r')
xlabel('Bifurcation parameter, $A$', 'interpreter', 'latex')
ylabel('$\Re(\lambda)$', 'interpreter', 'latex')
grid on
xlim([min(Aval),max(Aval)])
set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[5.697361111111110 3.033888888888889 10 5])
%exportgraphics(f,"../../../../Mattia_Sensi_Lukas_Eigentler_shared_space/klausmeier_max_real_part_vs_A.eps", "Resolution", 1000,'ContentType','vector' )

f2 = figure;
m = 0.005; A_s = 1.7; tval = (A_s-Aval)/m;
plot(tval,maxreal(sortind))
xlim([min(tval),max(tval)])
xlabel('Time, $t$', 'interpreter', 'latex')
ylabel('$\Re(\lambda)$', 'interpreter', 'latex')
grid on
set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',11)
set(f2,'Units','centimeters')
set(f2,'Position',[5.697361111111110 3.033888888888889 10 5])