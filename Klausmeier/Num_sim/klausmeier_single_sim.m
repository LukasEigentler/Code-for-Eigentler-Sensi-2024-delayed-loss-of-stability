%% Single simulation of the Klausmeier model with changing rainfall parameter
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023

set(0,'DefaultFigureWindowStyle','docked')
clear; 
close all
tic
f1 = figure(10);
withpause = 0; % include a pause in the parameter change
withtermination = 0; % include termination of change
withchangem = 0; % include a sudden change in m to a fixed value for all realisations
withchangeA = 0; % include sudden change of bifurcation parameter at the stability boundary
withsuddenchangeA = 1; % sudden change of A beyond stability boundary

%% Solver setup

%% time dependent parameter


m = 0.001; % rate of change of bifurcation parameter
calb_time = 100; %time taken to calibrate solution at start for withsuddenchangeA 
A_s = 1.7; A_e = 1.65;
tspan = [0,calb_time+(A_s-A_e)/m];
kt1 = linspace(1e-5,tspan(end),10000);
A_vec1 =  A_s-m*kt1; %Rainfall
kt0 = linspace(0,calb_time); A_vec0 = A_s*ones(1,length(kt0));
kt=[kt0,calb_time+kt1]; A_vec = [A_vec0,A_vec1];


P = 300; % Time after crossing stab boundary at which pause in parameter change occurs
Plength = 30000; % duration of pause in parameter change
mbar = 5e-4;
% Astab = 1.69688; % L = 20
Astab = 0.8382; % L = 40
if withpause == 1
    Astab = 1.69688;
    [~,boundaryind] = min(abs(A_vec - Astab));
    pauseind = find(kt<kt(boundaryind) + P); pauseind = pauseind(end);
    kt1 = kt(1:pauseind-1); A_vec1 = A_vec(1:pauseind-1);
    kt2 = linspace(kt(pauseind),kt(pauseind)+Plength, 1000); A_vec2 = A_vec(pauseind)*ones(1,1000);
    kt3 = kt(pauseind+1:end)+Plength; A_vec3 = A_vec(pauseind+1:end);
    kt = [kt1,kt2,kt3]; A_vec = [A_vec1,A_vec2, A_vec3];
    tspan(end) = 50000;
elseif withchangem == 1
    tspan = [0,(A_s-A_e)/mbar];
    kt = linspace(0,tspan(end),10000);
    A_vec =  A_s-m*kt; %Rainfall
    [~,boundaryind] = min(abs(A_vec - Astab));
    pauseind = find(kt<kt(boundaryind) + P); pauseind = pauseind(end);
    A_vec1 = A_vec(1:pauseind-1);
    A_vec2 = A_vec(pauseind) + mbar*kt(pauseind) - mbar*kt(pauseind:end);
    A_vec = [A_vec1,A_vec2];
elseif withtermination == 1
        tspan = linspace(0,1e4,1e4);
        [~,boundaryind] = min(abs(A_vec - Astab));
        pauseind = find(kt<kt(boundaryind) + P); pauseind = pauseind(end);
        A_vec1 = A_vec(1:pauseind-1);
        A_vec2 = A_vec1(end)*ones(1,length(A_vec(pauseind:end)));
        A_vec = [A_vec1,A_vec2];
elseif withsuddenchangeA == 1
    tspan = [0,10000];
    kt = linspace(0,tspan(end),10000);
    A_vec =  A_e*ones(1,length(kt)); %Rainfall

    calibrateind = find(kt<calb_time); calibrateind = calibrateind(end); % let solution calibrate for 100 time units, then change A
    kt1 = kt(1:calibrateind-1); A_vec1 = A_s*ones(1,length(A_vec(1:calibrateind-1)));
    kt2 = kt(calibrateind:end); A_vec2 = A_vec(calibrateind:end);
    kt = [kt1,kt2]; A_vec = [A_vec1,A_vec2];
end

%% other parameters
B = 0.45; %Plant loss 
nu = 182.5; % advection speed
d = 500; %Water diffusion


%%  Space and Time grid
xmax = 200; % half of the space domain
M = 2^9; % no of space points used in discretisation

%% Initialise space coordinates - different options, uncomment the desired one
x=linspace(-xmax,xmax,M);
dx = x(2)-x(1);

% perturbations of specific wave number(s)

% kmax = 0.23;
% pert = zeros(length(kmax),length(x));
% for kk = 1:length(kmax)
%     pert(kk,:) = sin(2*kk*(x-x(1))*pi./(x(end)-x(1)));
% end
% for xx = 1:length(x)
%     pertsum(xx) = sum(pert(:,xx));
% end
% % pertsum = pertsum./(10*max(abs(pertsum)));
% 
% 
% pertu = pertsum;
% pertw = pertsum;
%

% random perturbation
% pertu = 0.1*rand(1,M);
% pertw = 0.1*rand(1,M);
% % 
% % 
% vin(1:M) = 2 +pertu; 
% vin(M+1:2*M) = 1+pertw;

% load IC from file

% start_data= load("start_data.mat");
% xstart = start_data.x; vstart = start_data.v;
% vin(1:M) = interp1(xstart,vstart(end,1:end/2),x);
% vin(M+1:2*M) = interp1(xstart,vstart(end,end/2+1:end),x);


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


ustart = interp1(xstart,ustart,x); ustart(end) = ustart(1); ustart(ustart<0) = 0;
wstart = interp1(xstart,wstart,x); wstart(end) = wstart(1);
vin = [ustart,wstart];




%% ODE Solver

options = odeset('Stats', 'on', 'MaxStep',10, 'NonNegative',1:2*M); % max step size needed to capture changing parameter over time!
tic
[t_out,v] = ode15s(@(t,v) klausmeierode(t,v,B,nu,d,M,dx,kt,A_vec,x), tspan, vin, options);
toc
v_end = v(end,:);


%% Plots
lw = 2;
col = lines;

[~, timeind]  = min(abs(t_out-t_out(end)));
% timeind = 5499;

figure
subplot(2,1,1)
hold on
p(1) = plot(x,v(timeind,1:M),'Linewidth', lw, 'displayname', 'Grass density, $u_1(x,t)$', 'color',col(1,:));
ylabel({'Plants', '$u(x,t)$'}, 'Interpreter', 'Latex');
xlabel('Space, $x$', 'Interpreter', 'Latex');
grid on


subplot(2,1,2)
plot(x,v(timeind,M+1:2*M),'Linewidth', lw, 'color',col(4,:));
hold on
grid on
ylabel({'Water', '$w(x,t)$'}, 'Interpreter', 'Latex');
xlabel('Space, $x$', 'Interpreter', 'Latex');
set(findall(gcf,'-property','FontSize'),'FontSize',20)



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
figure(f1)
hold on
plot(t_out,L)
xlabel("Time")
ylabel("Wavelength")
A_plot = interp1(kt,A_vec,t_out);
yyaxis right
plot(t_out,A_plot,'linewidth',2)

ylabel('Resource input, $A$', 'interpreter','latex')

%% contour plot (sols in (x-t) plane)
f = figure;
plot_ind = intersect(find(t_out>000),find(t_out<100000));
A_plot = interp1(kt,A_vec,t_out);
subplot(2,1,1)
contourf(t_out(plot_ind),x,transpose(v(plot_ind,1:M)), 'linestyle', 'none')
colormap(flipud(summer))
shading interp
c = colorbar;
c.Ticks = [];
% xlabel('Time, $t$', 'interpreter','latex')
ylabel('Space, $x$', 'interpreter','latex')
% set(gca,'XTick',[])
% set(gca,'YTick',[-xmax,xmax])
% set(gca,"YTickLabel",["0", "L"])
% pbaspect([4 1 1])
subplot(2,1,2)


plot(t_out(plot_ind),A_plot(plot_ind),'linewidth',2)
hold on
Astab = 1.69688;
yline(Astab,'--')
xlabel('Time, $t$', 'interpreter','latex')
% ylabel('Rainfall, $A$', 'interpreter','latex')
ylabel('Bif. param.', 'interpreter','latex')
% set(gca,'YTick',[])
% set(gca,'YTickLabel',[])
col = lines;
% set(gca,'XTick',tspan)
% set(gca,"XTickLabel",["0", "T"])
set(gca,'YColor',col(1,:))
% pbaspect([4 1 1])
ylim([0,2])
xlim(tspan)

yyaxis right
plot(t_out,L)
% xlabel("Time")
ylabel("Wavelength")
% set(gca,'YTick',[20, 40, 100, 200])
ylim([0,60])
ax = gca(); 
% ax.YTick = [0, 20, 40, 100, 200];
yl = arrayfun(@(y)yline(ax, y,'LineStyle',ax.GridLineStyle,'Color',ax.GridColor,'Alpha',ax.GridAlpha),ax.YTick);


set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
% set(f,'Position',[5.697361111111110 3.033888888888889 21 10])
set(f,'Position',[5.697361111111110 3.033888888888889 8 8])
%set(f, 'Renderer', 'painters')
% exportgraphics(f,"../../Figures/numsim_klausmeier.eps", "Resolution", 500 )
% exportgraphics(f,"../../Figures/numsim_klausmeier_for_grant.eps", "Resolution", 1000 )
% exportgraphics(f,"../../Mattia_Sensi_Lukas_Eigentler_shared_space/numsim_klausmeier_long_delay.eps", "Resolution", 1000 )
