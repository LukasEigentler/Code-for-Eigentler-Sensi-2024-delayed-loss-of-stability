%% Klausmeier: Time delay simulations
% This script numerically calculates the time delay that occurs after
% crossing a stability boundary to the change in wavelength. For this, a 
% stable pattern is taken close to the stability boundary and then the 
% bifurcation parameter is changed suddenly.
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023

clear; 
% close all
plotonly = 1;

%% Parameters
Atarget = 1.65:-0.025:0.5; % target A values 
A_s = 1.7; % initial rainfall parameter at which sol is loaded from AUTO
Lstart = 20; L = Lstart; A = A_s; %wavelength
Astab = 1.69688; % bif para at crossing of boundary: not automated!! 
B = 0.45; %Plant loss 
nu = 182.5; % advection speed
d = 500; %Water diffusion
calb_time = 100; %time taken to calibrate solution at start
filename = "sudden_change_A_L"+num2str(Lstart);
%%  Space and Time grid
xmax = 200; % half of the space domain
M = 2^9; % no of space points used in discretisation
x=linspace(-xmax,xmax,M);
dx = x(2)-x(1);
options = odeset('Stats', 'on'); % max step size needed to capture changing parameter over time! - REVIEW as no gradual change in this case

%% load IC from numcont file

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

if plotonly ~=1
    %% loop over Atarget
    tdelay = NaN*ones(1,length(Atarget)); Ajump = tdelay;
    v_jump = NaN*ones(length(Atarget),2*M);
    for mm = 1:length(Atarget)
        try
        disp("Step "+num2str(mm)+" of "+num2str(length(Atarget)))
        if mm == 1
        tspan = [0,10000];
        else
            tspan = [0,max([200,calb_time+tdelay(mm-1)+100])];
        end
        kt = linspace(0,tspan(end),10000);
        A_vec =  Atarget(mm)*ones(1,length(kt)); %Rainfall

        calibrateind = find(kt<calb_time); calibrateind = calibrateind(end); % let solution calibrate for 100 time units, then change A
        kt1 = kt(1:calibrateind-1); A_vec1 = A_s*ones(1,length(A_vec(1:calibrateind-1)));
        kt2 = kt(calibrateind:end); A_vec2 = A_vec(calibrateind:end);
        kt = [kt1,kt2]; A_vec = [A_vec1,A_vec2];
        
    
        [t_out,v] = ode15s(@(t,v) klausmeierode(t,v,B,nu,d,M,dx,kt,A_vec), tspan, vin, options);
    
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
     % find delay
        tstab = find(A_vec<Astab); tstab = tstab(1); tstab = kt(tstab); 
        tjumpind = find(L>Lstart);
        if ~isempty(tjumpind)
            tjumpind = tjumpind(1); tjump = t_out(tjumpind);
            tdelay(mm) = tjump - tstab;
            Ajumpind = find(kt>tjump); Ajumpind = Ajumpind(1);
              % find solution at jump
            v_jump(mm,:) = v(tjumpind-1,:);
        end


   
        catch
            warning("Error occurred")
        end
    end
    %%
    
    

    try
        old_data = load(filename);
        Atarget = [old_data.Atarget,Atarget];
        [Atarget, sortind] = sort(Atarget);
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
Atarget = Atarget(~isnan(tdelay));
tdelay(isnan(tdelay)) = [];

f = figure;
loglog(Astab - Atarget, tdelay, '--o')
grid on
hold on
loglog(Astab - Atarget,10*(Astab - Atarget).^(-8/4), '--')

xlabel({"bifurcation parameter: distance to", "stability boundary, $A_{stab} - A_{target}$"}, "Interpreter","latex")
ylabel("Time delay, $t_{delay}$", "Interpreter","latex")
title("Instantaneous change")


set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[9 9 8 8])
%exportgraphics(f,"../../Figures/klausmeier_time_delay_L"+num2str(Lstart)+".eps", "Resolution", 500 )
Amdata = [Atarget;tdelay]; % save to dat file if used in batch calculation of spectra


save(filename, "Atarget", "tdelay", "v_jump", "maxrespec")

%% compare with spectrum data

cd ../Num_cont/spectra_calc_batch/spectrum_data/
Files = dir;

% extract files corresponding to correct wavelength
correctL = NaN*ones(1,length(Files));
for ff = 1:length(Files)
    correctL(ff) = contains(Files(ff).name, ["L_"+num2str(20),"L_"+num2str(L-1)+"dot99"]);
    filenames(ff) = string(Files(ff).name);
end
Files = Files(correctL==1);
filenames = filenames(correctL==1);
maxrespec = NaN*ones(1,length(Atarget)); lambda_area = maxrespec;
for aa = 1:length(Atarget)
    loadind = find(contains(filenames,"A_"+strrep(num2str(Atarget(aa)),".","dot")+"_and")); % find correct spectrum file
    try % some Atarget have no spectrum because no PTW with Atarget exists
        load(Files(loadind).name) %load spectrum data
        [maxrespec(aa),maxrealind] = max(spectrum_real); % extract max real part
        tempind1 = find(spectrum_real>0); tempind2 = find(spectrum_imag>0); areaind = intersect(tempind1,tempind2);
        lambda_area(aa) = polyarea(spectrum_real(areaind),spectrum_imag(areaind));
    end
end

f1 = figure;
loglog(maxrespec, tdelay, '--o')
grid on
hold on
loglog(maxrespec,6*maxrespec.^(-1), '--')
xlabel("Max real part of spectrum", "Interpreter","latex")
ylabel("Time delay, $t_{delay}$", "Interpreter","latex")
title("Instantaneous change")
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 9 8 8])

f2 = figure;
semilogx(maxrespec, Atarget, '--o')
grid on
hold on
xlabel("Max real part of spectrum", "Interpreter","latex")
ylabel("bifurcation parameter, $A_{target}$", "Interpreter","latex")
title("Instantaneous change")
set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',11)
set(f2,'Units','centimeters')
set(f2,'Position',[27 9 8 8])

%% delay against A integral

Aint = tdelay.*Atarget;


f3 = figure;
loglog(Aint, tdelay, '--o')
grid on
hold on

xlabel("$\overline{A}$", "Interpreter","latex")
ylabel("Time delay, $t_{delay}$", "Interpreter","latex")
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
ylabel("Time delay, $t_{delay}$", "Interpreter","latex")
title("Instantaneous change")


set(f4,'Windowstyle','normal')
set(findall(f4,'-property','FontSize'),'FontSize',11)
set(f4,'Units','centimeters')
set(f4,'Position',[0 0 8 8])

f5 = figure;
loglog(lambdaint,tdelay,'--o')
grid on
xlabel("$\overline{\lambda}$", "Interpreter","latex")
ylabel("Time delay, $t_{delay}$", "Interpreter","latex")
title("Instantaneous change")


set(f5,'Windowstyle','normal')
set(findall(f5,'-property','FontSize'),'FontSize',11)
set(f5,'Units','centimeters')
set(f5,'Position',[0 0 8 8])