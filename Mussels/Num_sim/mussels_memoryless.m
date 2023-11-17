%% Memorylessness test Mussel model
% This script shows that the delay is independent of the dynamics BEFORE
% hitting the stability boundary

clear; close all;

plotonly = 1;

if plotonly == 0
    %% other parameters
alpha = 50;
beta = 200;
eta = 0.1;
theta = 2.5;
nu = 360;
D = 1; %Water diffusion
calb_time = 100;
    
%%  Space and Time grid
xmax = 60; % half of the space domain
M = 2^8; % no of space points used in discretisation

%% Initialise space coordinates - different options, uncomment the desired one
x=linspace(-xmax,xmax,M);
dx = x(2)-x(1);



%% different start values of delta with sudden drop
mcol = logspace(0,1,10);
plotind = [1,5,10];
delta_e = 280;
Lstart = 15; L=Lstart; deltastab = 289.211;
for aa = 1:length(mcol)
    fprintf("Step "+num2str(aa) + " of " + num2str(length(mcol)))
    delta_s = 310; delta = delta_s;
    m = mcol(aa);
    

    tspan = [0,1000+ (delta_s-deltastab)/m];
    kt = linspace(0,tspan(end),10000);
    delta_vec =  delta_e*ones(1,length(kt)); %Rainfall
    delta_vec0 = delta_s - m*kt; delta_vec0(delta_vec0<deltastab) = [];
    kt0 = kt(1:length(delta_vec0));

    calibrateind = find(kt<calb_time); calibrateind = calibrateind(end); % let solution calibrate for 100 time units, then change A
    kt2 = kt(length(kt0)+1:end); delta_vec2 = delta_vec(length(kt0)+1:end);
    kt = [kt0,kt2]; delta_vec = [delta_vec0,delta_vec2];

% load IC from numcont file
    
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
    
    %% ODE Solver
    
    options = odeset('Stats', 'on', 'MaxStep',10, 'NonNegative',1:3*M); % max step size needed to capture changing parameter over time!
    tic
    [t_out,v] = ode15s(@(t,v) musselode(t,v,nu,M,dx,kt,delta_vec, alpha, beta, theta, D, eta), tspan, vin, options);
    toc
    v_end = v(end,:);

      % find wavelength
    L  = NaN*ones(1,length(t_out));
    for tt = 1:length(t_out)
        maxdens = max(v(tt,1:M));
        locmaxind = find(islocalmax([v(tt,M),v(tt,1:M),v(tt,1)]));
        locmaxind(locmaxind>M+1) = [];
        locmaxind(locmaxind == 1) = [];
        L(tt) = 2*xmax/length(find(v(tt,locmaxind) > 0));
    end
     % find delay  
    tstab = find(delta_vec<deltastab); tstab = tstab(1); tstab = kt(tstab); 
    tjumpind = find(L>Lstart); tjumpind = tjumpind(1); tjump = t_out(tjumpind);
    tdelay(aa) = tjump - tstab;

    % plot
    if ~isempty(intersect(aa,plotind))
        f(aa) = figure;
        plot_ind = intersect(find(t_out>000),find(t_out<100000));
        delta_plot = interp1(kt,delta_vec,t_out);
        subplot(2,1,1)
        contourf(t_out(plot_ind),x,transpose(v(plot_ind,1:M)), 'linestyle', 'none')
        colormap(flipud(summer))
        shading interp
        c = colorbar;
        c.Ticks = [];
        ylabel('Space, $x$', 'interpreter','latex')
        % set(gca,'XTick',[])
        % set(gca,'YTick',[-xmax,xmax])
        % set(gca,"YTickLabel",["0", "L"])
        % pbaspect([4 1 1])
        subplot(2,1,2)
        
        
        plot(t_out(plot_ind),delta_plot(plot_ind),'linewidth',2)
        hold on
        
        yline(deltastab,'--')
        xlabel('Time, $t$', 'interpreter','latex')
        ylabel('Bif. param.', 'interpreter','latex')
        % set(gca,'YTick',[])
        % set(gca,'YTickLabel',[])
        col = lines;
        % set(gca,'XTick',tspan)
        % set(gca,"XTickLabel",["0", "T"])
        set(gca,'YColor',col(1,:))
        % pbaspect([4 1 1])
        ylim([delta_e-10,delta_s+10])
        xlim(tspan)
        
        yyaxis right
        plot(t_out,L)
        % xlabel("Time")
        ylabel("Wavelength")
        % set(gca,'YTick',[20, 40, 100, 200])
        ylim([0,50])
        ax = gca(); 
        % ax.YTick = [0, 20, 40, 100, 200];
        yl = arrayfun(@(y)yline(ax, y,'LineStyle',ax.GridLineStyle,'Color',ax.GridColor,'Alpha',ax.GridAlpha),ax.YTick);
        
        %sgtitle("$A = "+num2str(A)+"$", "Interpreter","latex")
        set(f(aa),'Windowstyle','normal')
        set(findall(f(aa),'-property','FontSize'),'FontSize',11)
        set(f(aa),'Units','centimeters')
        % set(f,'Position',[5.697361111111110 3.033888888888889 21 10])
        set(f(aa),'Position',[5.697361111111110 3.033888888888889 8 8])
    end

  


end


try
    old_data = load("memoryless_data");
    mcol = [old_data.mcol,mcol];
    [mcol, sortind] = sort(mcol);
    tdelay = [old_data.tdelay, tdelay];
    tdelay = tdelay(sortind);
catch
    disp("First run detected. Creating new file...")
end

save("memoryless_data","mcol","tdelay","alpha","beta","eta","theta", "nu", "D", "delta_e")
else
    load("memoryless_data");
end


%% delay plot
f1 = figure;
loglog(mcol,tdelay,'o')
hold on
grid on
xlabel({"bif. parameter rate of change", "before stability boundary, $-m$"}, "Interpreter", "latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter", "latex")
ylim([10,1e4])
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[5 3 8 8])

exportgraphics(f1,"../../Mattia_Sensi_Lukas_Eigentler_shared_space/mussels_memoryless.eps", "Resolution", 1000,'ContentType','vector' )
