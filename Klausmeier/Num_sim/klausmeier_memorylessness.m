%% Memorylessness test
% This script shows that the delay is independent of the dynamics BEFORE
% hitting the stability boundary
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023

clear; close all;

plotonly = 0;

if plotonly == 0
    %% other parameters
    B = 0.45; %Plant loss 
    nu = 182.5; % advection speed
    d = 500; %Water diffusion
    calb_time = 100;
    
    %%  Space and Time grid
    xmax = 200; % half of the space domain
    M = 2^8; % no of space points used in discretisation
    
    %% Initialise space coordinates - different options, uncomment the desired one
    x=linspace(-xmax,xmax,M);
    dx = x(2)-x(1);
    
    
    
    %% different start values of A with sudden drop
    mcol = [1e-5:1e-5:9e-5];
    plotind = [];
    A_e = 1.4;
    Lstart = 20; L=Lstart; Astab = 1.69688;
    for aa = 1:length(mcol)
        fprintf("Step "+num2str(aa) + " of " + num2str(length(mcol)))
        A_s = 2; A = A_s;
        m = mcol(aa);
        
    
        tspan = [0,1000+ (A_s-A_e)/m];
        kt = linspace(0,tspan(end),10000);
        A_vec =  A_e*ones(1,length(kt)); %Rainfall
        A_vec0 = A_s - m*kt; A_vec0(A_vec0<Astab) = [];
        kt0 = kt(1:length(A_vec0));
    
        calibrateind = find(kt<calb_time); calibrateind = calibrateind(end); % let solution calibrate for 100 time units, then change A
        kt2 = kt(length(kt0)+1:end); A_vec2 = A_vec(length(kt0)+1:end);
        kt = [kt0,kt2]; A_vec = [A_vec0,A_vec2];
    
    % load IC from numcont file
        
        start_data_op = importsol_klausmeier("../Num_cont/Pattern_generation/output/s.ptw_A"+strrep(num2str(A),'.','dot')+"_L"+strrep(num2str(Lstart),'.','dot'));
        start_data_op = start_data_op(1:end-4,:);
        
        ustart = repmat(start_data_op(:,2), 2*xmax/Lstart,1);
        wstart = repmat(start_data_op(:,4), 2*xmax/Lstart,1);
        xstart = [];
        for qq = 1:2*xmax/Lstart
            xstart = [xstart;(qq-1)*Lstart+Lstart*start_data_op(:,1)];
        end
        xstart = xstart - xmax;
        
        
        ustart = interp1(xstart,ustart,x); ustart(end) = ustart(1); ustart(ustart<0) = 0;
        wstart = interp1(xstart,wstart,x); wstart(end) = wstart(1);
        vin = [ustart,wstart];
        
        %% ODE Solver
        
        options = odeset('Stats', 'on', 'MaxStep',0.1, 'NonNegative',1:2*M); % max step size needed to capture changing parameter over time!
        tic
        [t_out,v] = ode15s(@(t,v) klausmeierode(t,v,B,nu,d,M,dx,kt,A_vec,x), tspan, vin, options);
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
        tstab = find(A_vec<Astab); tstab = tstab(1); tstab = kt(tstab); 
        tjumpind = find(L>Lstart); tjumpind = tjumpind(1); tjump = t_out(tjumpind);
        tdelay(aa) = tjump - tstab;
    
        % plot
        if ~isempty(intersect(aa,plotind))
            f(aa) = figure;
            plot_ind = intersect(find(t_out>000),find(t_out<100000));
            A_plot = interp1(kt,A_vec,t_out);
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
            
            
            plot(t_out(plot_ind),A_plot(plot_ind),'linewidth',2)
            hold on
            
            yline(Astab,'--')
            xlabel('Time, $t$', 'interpreter','latex')
            ylabel('Bif. param.', 'interpreter','latex')
            % set(gca,'YTick',[])
            % set(gca,'YTickLabel',[])
            col = lines;
            % set(gca,'XTick',tspan)
            % set(gca,"XTickLabel",["0", "T"])
            set(gca,'YColor',col(1,:))
            % pbaspect([4 1 1])
            ylim([0,A_s])
            xlim(tspan)
            
            yyaxis right
            plot(t_out,L)
            % xlabel("Time")
            ylabel("Wavelength")
            % set(gca,'YTick',[20, 40, 100, 200])
            ylim([0,40])
            ax = gca(); 
            % ax.YTick = [0, 20, 40, 100, 200];
            yl = arrayfun(@(y)yline(ax, y,'LineStyle',ax.GridLineStyle,'Color',ax.GridColor,'Alpha',ax.GridAlpha),ax.YTick);
            
            %sgtitle("$A = "+num2str(A)+"$", "Interpreter","latex")
            set(f(aa),'Windowstyle','normal')
            set(findall(f(aa),'-property','FontSize'),'FontSize',11)
            set(f(aa),'Units','centimeters')
            % set(f,'Position',[5.697361111111110 3.033888888888889 21 10])
            set(f(aa),'Position',[5.697361111111110 3.033888888888889 8 8])
            % exportgraphics(f(1),"../../Mattia_Sensi_Lukas_Eigentler_shared_space/klausmeier_memoryless_example_1e-4.eps", "Resolution", 1000)
            % exportgraphics(f(2),"../../Mattia_Sensi_Lukas_Eigentler_shared_space/klausmeier_memoryless_example_1e-2.eps", "Resolution", 1000)
            % exportgraphics(f(3),"../../Mattia_Sensi_Lukas_Eigentler_shared_space/klausmeier_memoryless_example_1e0.eps", "Resolution", 1000)
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
    
    save("memoryless_data","mcol","tdelay","B","nu","d","A_e")
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

exportgraphics(f1,"../../Mattia_Sensi_Lukas_Eigentler_shared_space/klausmeier_memoryless.eps", "Resolution", 1000,'ContentType','vector' )
