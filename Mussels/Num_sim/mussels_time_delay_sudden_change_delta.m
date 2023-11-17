%% Mussels: Time delay simulations
% This script numerically calculates the time delay that occurs after
% crossing a stability boundary to the change in wavelength. For this, a 
% stable pattern is taken close to the stability boundary and then the 
% bifurcation parameter is changed suddenly.

clear; 
close all
plotonly = 1;

%% Parameters
Lstart = 15; 
if Lstart == 15
    deltastab = 289.1997; % L=15
    delta_s = 290; % initial rainfall parameter at which sol is loaded from AUTO
elseif Lstart == 20
    deltastab = 261.674; % L=20
    delta_s = 262; % initial rainfall parameter at which sol is loaded from AUTO
end
deltadist = 21:1:23;
deltatarget = deltastab - deltadist; % target A values 

L = Lstart; delta = delta_s; %wavelength

alpha = 50;
beta = 200;
eta = 0.1;
theta = 2.5;
nu = 360;
D = 1; 
calb_time = 100; %time taken to calibrate solution at start
filename = "sudden_change_delta_L"+num2str(Lstart);
%%  Space and Time grid
xmax = 30; % half of the space domain
M = 2^8; % no of space points used in discretisation
x=linspace(-xmax,xmax,M);
dx = x(2)-x(1);
options = odeset('Stats', 'off'); % max step size needed to capture changing parameter over time! - REVIEW as no gradual change in this case

%% load IC from numcont file



if plotonly ~=1
    %% loop over deltatarget
    tdelay = NaN*ones(1,length(deltatarget)); 
    v_jump = NaN*ones(length(deltatarget),3*M);
    for mm = 1:length(deltatarget)

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

        Lchange = 0;
        vcol = zeros(1,3*M);
        tcol = [];
        Lcol = [];
        delta_plot = [];
        disp("Step "+num2str(mm)+" of "+num2str(length(deltatarget)))
        run = 0;
        while Lchange == 0
            run = run+1;
            tspan = [0,500];
            if run  > 1
                fprintf("Simulating to t = "+num2str(tcol(end)+tspan(end))+"\n")
            else
                fprintf("Simulating to t = "+num2str(tspan(end))+"\n")
            end
            kt = linspace(0,tspan(end),10000);
            delta_vec =  deltatarget(mm)*ones(1,length(kt)); %Rainfall
            if run == 1
                calibrateind = find(kt<calb_time); calibrateind = calibrateind(end); % let solution calibrate for 100 time units, then change A
                kt1 = kt(1:calibrateind-1); delta_vec1 = delta_s*ones(1,length(delta_vec(1:calibrateind-1)));
                kt2 = kt(calibrateind:end); delta_vec2 = delta_vec(calibrateind:end);
                kt = [kt1,kt2]; delta_vec = [delta_vec1,delta_vec2];
            end
        
            [t_out,v] = ode15s(@(t,v) musselode(t,v,nu,M,dx,kt,delta_vec, alpha, beta, theta, D, eta), tspan, vin, options);
            v_end = v(end,:);
            vcol = [vcol;v];
            tcol = [tcol;(run-1)*tspan(end)+t_out];
        
        % find wavelength
            L  = NaN*ones(1,length(t_out));
            for tt = 1:length(t_out)
                maxdens = max(v(tt,1:M));
                if maxdens > 1e-2 % if desert do nothing
                    locmaxind = find(islocalmax([v(tt,M),v(tt,1:M),v(tt,1)]));
                    locmaxind(locmaxind>M+1) = [];
                    locmaxind(locmaxind == 1) = [];
                    L(tt) = 2*xmax/length(find(v(tt,locmaxind) > 0));
                end
            end
            Lcol = [Lcol,L];
            delta_plot = [delta_plot; interp1(kt,delta_vec,t_out)];
         % find delay
            tstab = find(delta_vec<deltastab); tstab = tstab(1); tstab = kt(tstab); 
            tjumpind = find(L>Lstart);
            if ~isempty(tjumpind)
                tjumpind = tjumpind(1); tjump = t_out(tjumpind);
                tdelay(mm) = (run-1)*tspan(end)+tjump - tstab;
                deltajumpind = find(kt>tjump); deltajumpind = deltajumpind(1);
                  % find solution at jump
                v_jump(mm,:) = v(tjumpind-1,:);
                Lchange = 1;
            else
                vin = v_end;
            end
        end

   
        % catch
        %     warning("Error occurred")
        % end
    end
    %%
    
    

    try
        old_data = load(filename);
        deltatarget = [old_data.deltatarget,deltatarget];
        [deltatarget, sortind] = sort(deltatarget);
        tdelay = [old_data.tdelay, tdelay];
        tdelay = tdelay(sortind);
        v_jump = [old_data.v_jump;v_jump];
        v_jump = v_jump(sortind,:);
    catch
        disp("First run detected. Creating new file...")
    end
else
    load(filename);
end
deltatarget = deltatarget(~isnan(tdelay));
tdelay(isnan(tdelay)) = [];

f = figure;
loglog(deltastab - deltatarget, tdelay, '--o')
grid on
hold on
loglog(deltastab - deltatarget,10*(deltastab - deltatarget).^(-8/4), '--')

xlabel({"bifurcation parameter: distance to", "stability boundary, $\delta_{\textrm{stab}} - \delta_{\textrm{target}}$"}, "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
title("Instantaneous change")


set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[9 9 8 8])
%exportgraphics(f,"../../Figures/mussels_time_delay_L"+num2str(Lstart)+".eps", "Resolution", 500 )
deltamdata = [deltatarget;tdelay]; % save to dat file if used in batch calculation of spectra




%% compare with spectrum data

cd ../Num_cont/spectra_calc_batch/spectrum_data/
Files = dir;

% extract files corresponding to correct wavelength
correctL = NaN*ones(1,length(Files));
for ff = 1:length(Files)
    correctL(ff) = contains(Files(ff).name, ["L_"+num2str(Lstart),"L_"+num2str(L-1)+"dot99"]);
    filenames(ff) = string(Files(ff).name);
end
Files = Files(correctL==1);
filenames = filenames(correctL==1);
maxrespec = NaN*ones(1,length(deltatarget)); lambda_area = maxrespec;
for aa = 1:length(deltatarget)
    loadind = find(contains(filenames,"delta_"+strrep(num2str(deltatarget(aa)),".","dot")+"_and")); % find correct spectrum file
    try % some Atarget have no spectrum because no PTW with Atarget exists
        load(Files(loadind).name) %load spectrum data
        [maxrespec(aa),maxrealind] = max(spectrum_real); % extract max real part
        tempind1 = find(spectrum_real>0); tempind2 = find(spectrum_imag>0); areaind = intersect(tempind1,tempind2);
        lambda_area(aa) = polyarea(spectrum_real(areaind),spectrum_imag(areaind));
    end
end
cd ../../../Num_sim/
save(filename, "deltatarget", "tdelay", "v_jump", "maxrespec")

f1 = figure;
loglog(maxrespec, tdelay, '--o')
grid on
hold on
loglog(maxrespec,6*maxrespec.^(-1), '--')
xlabel("Max real part of spectrum", "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
title("Instantaneous change")
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 9 8 8])

f2 = figure;
semilogx(maxrespec, deltatarget, '--o')
grid on
hold on
xlabel("Max real part of spectrum", "Interpreter","latex")
ylabel("bifurcation parameter, $\delta_{\textrm{target}}$", "Interpreter","latex")
title("Instantaneous change")
set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',11)
set(f2,'Units','centimeters')
set(f2,'Position',[27 9 8 8])

%% delay against A integral

deltaint = tdelay.*deltatarget;


f3 = figure;
loglog(deltaint, tdelay, '--o')
grid on
hold on

xlabel("$\overline{\delta}$", "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
title("Instantaneous change")


set(f3,'Windowstyle','normal')
set(findall(f3,'-property','FontSize'),'FontSize',11)
set(f3,'Units','centimeters')
set(f3,'Position',[0 9 8 8])

%% delay against spectrum integral

muint = maxrespec.*tdelay;
lambdaint = lambda_area.*tdelay;

f4 = figure;
loglog(muint,tdelay,'--o')
grid on
xlabel("$\overline{\mu}$", "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
title("Instantaneous change")


set(f4,'Windowstyle','normal')
set(findall(f4,'-property','FontSize'),'FontSize',11)
set(f4,'Units','centimeters')
set(f4,'Position',[0 0 8 8])

f5 = figure;
loglog(lambdaint,tdelay,'--o')
grid on
xlabel("$\overline{\lambda}$", "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
title("Instantaneous change")


set(f5,'Windowstyle','normal')
set(findall(f5,'-property','FontSize'),'FontSize',11)
set(f5,'Units','centimeters')
set(f5,'Position',[0 0 8 8])