%% Klausmeier model delay predictions
% This script perfomrs predictions for the time delay that occurs after a
% PTW crosses a stability boundary in the Busse balloon and before it
% undergoes a wavelength change.
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023

clear;
close all;

%% Key input parameters
A_s = 1.7; % start value of bif para
Lstart = 20; % wavelength of PTW
Astab = 1.69688;
kt = linspace(0,1e5,1e4); % time domain
m = 1e-4; % m value for constant rate of change

%% sinosoidal regime
% A_vec = A_s-m*kt.*(2+sin(kt/50)); % sinosoidal regime supplement

%% linear regime
% A_vec = A_s-m*kt; % linear regime

%% change of direction regime
% swapt = 120; tempind = find(kt<swapt); A_vec(tempind) = A_s - m*kt(tempind);
% tempind2 = find(kt>=swapt); A_vec(tempind2) = A_s-swapt*m + m*(kt(tempind2)-swapt);

%% zig zag regime
% swapt1 = 60; tempind = find(kt<swapt1); A_vec(tempind) = A_s - m*kt(tempind);
% tempind2 = find(kt>=swapt1); A_vec(tempind2) = A_s-swapt1*m + m*(kt(tempind2)-swapt1);
% swapt2 = swapt1 + 50; 
% tempind3 = find(kt>=swapt2); A_vec(tempind3) = A_s-swapt1*m + m*(swapt2 - swapt1) - m*(kt(tempind3)-swapt2);

%% zig zag loop
A_s = 1.75; % start value of bif para
swaptloop = 1800;
tcount = 0;
A_vec = A_s;
while tcount < max(kt)
    m = 0.0001;
    tempind = find(kt>tcount); A_vec(tempind) = A_vec(tempind(1)-1) - m*(kt(tempind) - kt(tempind(1)));
    try
    tempind = find(kt>tcount+swaptloop); A_vec(tempind) = A_vec(tempind(1)-1) + m*(kt(tempind) - kt(tempind(1)));
    end
    tcount = tcount + 2*swaptloop;
end

%% linear regime then sudden increase
% increaset = 115;
% A_vec = A_s-m*kt; A_vec(kt>increaset) = A_s+0.3;

%% find tstab
tstabind = find(A_vec<Astab); tstab_ind = tstabind(1);
tstab = kt(tstab_ind);

%% model parameters
B = 0.45; %Plant loss 
nu = 182.5; % advection speed
d = 500; %Water diffusion

%% load spectra info
cd ../Num_cont/spectra_calc_batch/spectrum_data/
Files = dir;
% find spectra data for correct wavelength
correctL = NaN*ones(1,length(Files));
for ff = 1:length(Files)
    correctL(ff) = contains(Files(ff).name, ["L_"+num2str(Lstart),"L_"+num2str(Lstart-1)+"dot99"]);
    filenames(ff) = string(Files(ff).name);
end
Files = Files(correctL==1);
filenames = filenames(correctL==1);
% determine max real part of spectra
maxrespec = NaN*ones(1,length(Files)); Aval_spec = maxrespec;
for aa = 1:length(Files)
    loadind = aa;
    load(Files(loadind).name) %load spectrum data
    [maxrespec(aa),maxrealind] = max(spectrum_real); % extract max real part

    filename = Files(aa).name;
    startind = strfind(filename,'_A_'); startind = startind+3;
    endind = strfind(filename,'_and_L_'); endind = endind-1;
    Aval_spec(aa) = str2num(strrep(filename(startind:endind),'dot','.')); %determine A value
end
%% calculate delay prediction
% intepolation of spectra values
tempind1 = find(A_vec>max(Aval_spec)); 
tempind = find(A_vec>=min(Aval_spec)); 
% tempind = intersect(tempind1,tempind2); %indices for which spectrum data is available
A_vec_stab = A_vec(tempind); %extract A values for which spectra data is available
t_vec_stab = kt(tempind); % and associated times
maxrespec_inter = interp1(Aval_spec,maxrespec,A_vec_stab);
maxrespec_inter(tempind1) = 0;

% calculation of mubar without reset

for ii = 2:length(t_vec_stab)
    mubar_noreset(ii) = trapz(t_vec_stab(1:ii),maxrespec_inter(1:ii));
end

% find delay prediction
delayind = find(mubar_noreset>10); 
if ~isempty(delayind)
    delayind = delayind(1);
    delay_pred_no_reset = t_vec_stab(delayind)-t_vec_stab(1);
    xline(delay_pred_no_reset+t_vec_stab(1),'--r')
else
    delay_pred = NaN;
end

% calculation of mubar with reset
for ii = 2:length(t_vec_stab)
    if maxrespec_inter(ii) == 0
        mubar(ii) = 0;
    else
        zeroind = find(maxrespec_inter==0);
        zeroindsmaller = find(zeroind<ii);
        startind = zeroind(zeroindsmaller(end));
        mubar(ii) = trapz(t_vec_stab(startind:ii),maxrespec_inter(startind:ii));
    end
end


% find delay prediction
delayind = find(mubar>10); 
if ~isempty(delayind)
    delayind = delayind(1);
    delay_pred = t_vec_stab(delayind)-t_vec_stab(1);
    xline(delay_pred+t_vec_stab(1),'--r')
else
    delay_pred = NaN;
end
fprintf("The delay prediction is "+num2str(delay_pred)+"\n")

% plot mubar over time
f = figure;
subplot(2,2,1)
plot(kt,A_vec)
grid on
hold on
yline(Astab, '--k')
% xline(tstab,'--r')
xlim([0,inf])
xlabel("Time, $t$", "Interpreter", "latex")
ylabel("Bif. parameter, $A(t)$", "Interpreter","latex")
% title("Predicted delay = "+num2str(delay_pred,'%.0f'),"Interpreter","latex")
title({"Predicted delay = "+num2str(delay_pred,'%.0f'),"Predicted delay without reset = "+num2str(delay_pred_no_reset,'%.0f')},"Interpreter","latex")

subplot(2,2,3)
semilogy(t_vec_stab,mubar)
hold on
grid on
semilogy(t_vec_stab,mubar_noreset)
yline(10,'--k')
% xline(delay_pred,'--r')
xlabel("Time, $t$", "Interpreter", "latex")
ylabel("$\overline{\mu}$", "Interpreter","latex")
xlim([0,inf])
ylim([1e-3,inf])
legend("reset", "no reset", "Location", "southeast")






%% compare with simulation
cd ../../../Num_sim/
% Space  grid
xmax = 200; % half of the space domain
M = 2^8; % no of space points used in discretisation
x=linspace(-xmax,xmax,M);
dx = x(2)-x(1);

% load IC from numcont file
L = 20; A = A_s;
start_data_op = importsol_klausmeier("../Num_cont/Pattern_generation/output/s.ptw_A"+strrep(num2str(A),'.','dot')+"_L"+strrep(num2str(L),'.','dot'));
start_data_op = start_data_op(1:end-4,:);

ustart = repmat(start_data_op(:,2), 2*xmax/L,1);
wstart = repmat(start_data_op(:,4), 2*xmax/L,1);
xstart = [];
for qq = 1:2*xmax/L
    xstart = [xstart;(qq-1)*L+L*start_data_op(:,1)];
end
xstart = xstart - xmax;

ustart = interp1(xstart,ustart,x); ustart(end) = ustart(1);
wstart = interp1(xstart,wstart,x); wstart(end) = wstart(1);
vin = [ustart,wstart];

% ODE Solver

options = odeset('Stats', 'off'); % max step size needed to capture changing parameter over time!
tic
[t_out,v] = ode15s(@(t,v) klausmeierode(t,v,B,nu,d,M,dx,kt,A_vec,x), [0,kt(end)], vin, options);
toc
v_end = v(end,:);

% Plots
lw = 2;
col = lines;

[~, timeind]  = min(abs(t_out-t_out(end)));
% timeind = 5499;


%% find wavelength
L  = NaN*ones(1,length(t_out));
for tt = 1:length(t_out)
    maxdens = max(v(tt,1:M));
    if maxdens > 1e-2 % if desert, do nothing
        locmaxind = find(islocalmax([v(tt,M),v(tt,1:M),v(tt,1)]));
        locmaxind(locmaxind>M+1) = [];
        locmaxind(locmaxind == 1) = [];
        L(tt) = 2*xmax/length(find(v(tt,locmaxind) > 0));
    end

end

% find delay
tjumpind = find(L>Lstart); 
if ~isempty(tjumpind)
    tjumpind = tjumpind(1); tjump = t_out(tjumpind);
    tdelay_sim = tjump - tstab;
else
    tdelay_sim = NaN;
end
fprintf("The delay detected in the simulation is "+num2str(tdelay_sim)+"\n")


%% contour plot (sols in (x-t) plane)
plot_ind = intersect(find(t_out>000),find(t_out<100000));
A_plot = interp1(kt,A_vec,t_out);
subplot(2,2,2)
contourf(t_out(plot_ind),x,transpose(v(plot_ind,1:M)), 'linestyle', 'none')
colormap(flipud(summer))
shading interp
% c = colorbar;
% c.Ticks = [];
% xlabel('Time, $t$', 'interpreter','latex')
ylabel('Space, $x$', 'interpreter','latex')
% set(gca,'XTick',[])
% set(gca,'YTick',[-xmax,xmax])
% set(gca,"YTickLabel",["0", "L"])
% pbaspect([4 1 1])
title("Observed delay = "+num2str(tdelay_sim,'%.0f'), 'Interpreter','latex')
subplot(2,2,4)
hold on
xlabel('Time, $t$', 'interpreter','latex')
plot(t_out,L)
% xlabel("Time")
ylabel("Wavelength")
% set(gca,'YTick',[20, 40, 100, 200])
% ylim([0,40])

grid on
set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[5.697361111111110 3.033888888888889 16 8])
% set(f,'Position',[5.697361111111110 3.033888888888889 25 8])
% exportgraphics(f,"../../Mattia_Sensi_Lukas_Eigentler_shared_space/klausmeier_delay_predictions_reset.eps", "Resolution", 1000 )

