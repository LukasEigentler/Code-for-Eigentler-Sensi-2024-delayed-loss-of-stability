%% Investigation of possible coexistence in the parameter region where both patterned states are stable
set(0,'DefaultFigureWindowStyle','docked')
clear; 
close all
tic
f1 = figure(10);

withsuddenchangedelta = 0; % sudden change of A beyond stability boundary
untilLchange = 0; % set to run simulation until first wavelength change is detected.

%% bif para dependent on time
m = 0.0001; % rate of change of bifurcation parameter
calb_time = 100; %time taken to calibrate solution at start for withsuddenchangeA 
deltastab = 289.1997; % L=15
% deltastab = 261.674; % L=20
delta_s = 290; delta_e = 270;
Lstart = 15;

%% Other parameters
alpha = 50;
beta = 200;
eta = 0.1;
theta = 2.5;
nu = 360;
D = 1; %Water diffusion


%%  Space and Time grid
xmax = 30; %Space domain size of interest
M = 2^8;


%% Initialise space coordinates
x=linspace(-xmax,xmax,M);
dx = x(2)-x(1);

%% load IC from numcont file
delta = delta_s;
start_data_op = importsol_mussels("../Num_cont/Pattern_generation/output/s.ptw_delta"+strrep(num2str(delta),'.','dot')+"_L"+strrep(num2str(Lstart),'.','dot'));
start_data_op = start_data_op(1:end-4,:);

mstart = repmat(start_data_op(:,2), 2*xmax/Lstart,1);
astart = repmat(start_data_op(:,4), 2*xmax/Lstart,1);
sstart = repmat(start_data_op(:,5), 2*xmax/Lstart,1);

xstart = [];
for qq = 1:2*xmax/Lstart
    xstart = [xstart;(qq-1)*Lstart+Lstart*start_data_op(:,1)];
end
xstart = xstart - xmax;


mstart = interp1(xstart,mstart,x); mstart(end) = mstart(1); mstart(mstart<0) = 0;
astart = interp1(xstart,astart,x); astart(end) = astart(1);
sstart = interp1(xstart,sstart,x); sstart(end) = sstart(1);
vin = [mstart,astart,sstart];

%% time dependent parameter


Lchange = 0;
run = 0;
vcol = zeros(1,3*M);
tcol = [];
Lcol = [];
delta_plot = [];

while Lchange == 0 % loop until wavelength change detected
    run = run+1;  
    if withsuddenchangedelta == 1
        tspan = [0,3000];
        if run  > 1
            fprintf("Simulating to t = "+num2str(tcol(end)+tspan(end))+"\n")
        else
            fprintf("Simulating to t = "+num2str(tspan(end))+"\n")
        end
        kt = linspace(0,tspan(end),10000);
        delta_vec =  delta_e*ones(1,length(kt)); %Rainfall
        if run == 1 % calibration in first run
            calibrateind = find(kt<calb_time); calibrateind = calibrateind(end); % let solution calibrate for 100 time units, then change A
            kt1 = kt(1:calibrateind-1); delta_vec1 = delta_s*ones(1,length(delta_vec(1:calibrateind-1)));
            kt2 = kt(calibrateind:end); delta_vec2 = delta_vec(calibrateind:end);
            kt = [kt1,kt2]; delta_vec = [delta_vec1,delta_vec2];
        end   
    else
        tspan = [0,calb_time+(delta_s-delta_e)/m];
        % tspan = [0,1000];
        kt1 = linspace(1e-5,tspan(end),10000);
        delta_vec1 =  delta_s-m*kt1; %Rainfall
        kt0 = linspace(0,calb_time); delta_vec0 = delta_s*ones(1,length(kt0));
        kt=[kt0,calb_time+kt1]; delta_vec = [delta_vec0,delta_vec1];
    end
    
    
    %% ODE Solver
    
    options = odeset('Stats', 'off');
    tic
    [t_out,v] = ode15s(@(t,v) musselode(t,v,nu,M,dx,kt,delta_vec, alpha, beta, theta, D, eta), tspan, vin, options);
    toc
    v_end = v(end,:);
    vcol = [vcol;v];
    tcol = [tcol;(run-1)*tspan(end)+t_out];
    
    
    %% find wavelength
    L  = zeros(1,length(t_out));
    for tt = 1:length(t_out)
        maxdens = max(v(tt,1:M));
        locmaxind = find(islocalmax([v(tt,M),v(tt,1:M),v(tt,1)]));
        locmaxind(locmaxind>M+1) = [];
        locmaxind(locmaxind == 1) = [];
        L(tt) = 2*xmax/length(find(v(tt,locmaxind) > 0));
        if L(tt) <10
            L(tt) = 0;
        end
        % if length(locmaxind)>2
        %     maxdif = diff(x(locmaxind));
        %     L(tt) = mean(maxdif(3:end-1));
        % end
    end
    Lcol = [Lcol,L];
    
    %figure(f1)
    %hold on
    %plot(tcol,Lcol)
    %ylabel("Wavelength")
    %yyaxis right
    delta_plot = [delta_plot; interp1(kt,delta_vec,t_out)];
    %plot(tcol,delta_plot,'linewidth',2)
    %xlabel("Time")
    %ylabel("\delta")
    
    %% find delay
    
    tstab = find(delta_vec<deltastab); tstab = tstab(1); tstab = kt(tstab); 
    tjumpind = find(L>Lstart); 
    if ~isempty(tjumpind)
        tjumpind = tjumpind(1); tjump = t_out(tjumpind);
        tdelay = (run-1)*tspan(end)+tjump - tstab;
        deltajumpind = find(kt>tjump); deltajumpind = deltajumpind(1);
        deltajump = delta_vec(deltajumpind);
        Lchange = 1;
    else
        vin = v_end;
    end

    if untilLchange == 0
        Lchange = 1;
    end
end
%% contour plot (sols in (x-t) plane)
f = figure;
v = vcol(2:end,:);
t_out = tcol;
subplot(2,1,1)
[t_out,tind] = unique(t_out);
v = v(tind,:);
Lcol = Lcol(tind);
plot_ind = intersect(find(t_out>0000),find(t_out<100000));
contourf(t_out(plot_ind),x,transpose(v(plot_ind,1:M)), 'linestyle', 'none')
colormap(flipud(summer))
shading interp

% colorbar
% xlabel('Time, $t$', 'interpreter','latex')
ylabel('Space, $x$', 'interpreter','latex')
% set(gca,'XTick',[])
% set(gca,'YTick',[-xmax,xmax])
% set(gca,"YTickLabel",["0", "L"])
% title('Mussels, $m(x,t)$', 'interpreter','latex')

subplot(2,1,2)


plot(t_out(plot_ind),delta_plot(plot_ind),'linewidth',2)
xlabel('Time, $t$', 'interpreter','latex')
% ylabel('Max. mussel growth rate, $\delta$', 'interpreter','latex')
ylabel('Bif. param.', 'interpreter','latex')
col = lines;
ylim([250,300])
% set(gca,'XTick',tspan)
% set(gca,"XTickLabel",["0", "T"])
set(gca,'YColor',col(1,:))
% set(gca,'YTick',[])


yyaxis right
plot(t_out,Lcol)
ylabel("Wavelength")
ylim([0,50])
ax = gca(); 
ax.YTick = [0,15, 20,30,40,50];
yl = arrayfun(@(y)yline(ax, y,'LineStyle',ax.GridLineStyle,'Color',ax.GridColor,'Alpha',ax.GridAlpha),ax.YTick);

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[5.697361111111110 3.033888888888889 10 9.5])
%set(f, 'Renderer', 'painters')
%exportgraphics(f,"../../Figures/numsim_mussels.eps", "Resolution", 500 )
%exportgraphics(f,"../../Figures/numsim_mussels_for_grant.eps", "Resolution", 1000 )
% exportgraphics(f,"../../Mattia_Sensi_Lukas_Eigentler_shared_space/numsim_mussels_long_delay.eps", "Resolution", 1000 )
