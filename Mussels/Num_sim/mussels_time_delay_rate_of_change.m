%% Mussels: Time delay simulations
% This script numerically calculates the time delay that occurs after
% crossing a stability boundary to the change in wavelength for different
% constant rates of change of the bifurcation parameter.

% keep f4
clear; 
close all
plotonly = 1;

%% Parameters
Lstart = 15; L = Lstart;  %wavelength

mcol = 0.0001:0.0001:0.0009; % rate of change A = As - mt



alpha = 50;
beta = 200;
eta = 0.1;
theta = 2.5;
nu = 360;
D = 1; 
calb_time = 100;



if Lstart == 15
    deltastab = 289.1997; % L=15
    delta_s = 290; delta_e = 270; delta = delta_s;
    xmax = 30; % half of the space domain
elseif Lstart == 20
    deltastab = 261.674; % L=20
    delta_s = 262; delta_e = 240; delta = delta_s;
    xmax = 50; % half of the space domain
end



%%  Space and Time grid

M = 2^8; % no of space points used in discretisation
x=linspace(-xmax,xmax,M);
dx = x(2)-x(1);
options = odeset('Stats', 'off'); % max step size needed to capture changing parameter over time!

%% load IC from numcont file

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

if plotonly ~=1
    %% loop over m
    tdelay = NaN*ones(1,length(mcol)); deltajump = tdelay;
    v_jump = NaN*ones(length(mcol),3*M);
    for mm = 1:length(mcol)
        % try
        disp("Step "+num2str(mm)+" of "+num2str(length(mcol)))
        m = mcol(mm);
        % tspan = [0,calb_time+(A_s-A_e)/m];
        if mm == 1
            tspan = [0,1e5];
        else
            tspan = [0,calb_time + (delta_s-deltastab)/m + 1.1*tdelay(mm-1)];
        end
        kt1 = linspace(1e-5,tspan(end),10000);
        delta_vec1 =  delta_s-m*kt1; %Rainfall
        kt0 = linspace(0,calb_time); delta_vec0 = delta_s*ones(1,length(kt0));
        kt=[kt0,calb_time+kt1]; delta_vec = [delta_vec0,delta_vec1];
    
        [t_out,v] = ode15s(@(t,v) musselode(t,v,nu,M,dx,kt,delta_vec, alpha, beta, theta, D, eta), tspan, vin, options);
    
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
        tdelay(mm) = tjump - tstab;
        deltajumpind = find(kt>tjump); deltajumpind = deltajumpind(1);
        deltajump(mm) = delta_vec(deltajumpind);

     % find solution at jump
        v_jump(mm,:) = v(tjumpind-1,:);
        % catch
        %     warning("Error occurred")
        % end
    end
    %%
    filename = "jumpdata_L"+num2str(Lstart);
    try
        old_data = load(filename);
        mcol = [old_data.mcol,mcol];
        [mcol, sortind] = sort(mcol);
        deltajump = [old_data.deltajump, deltajump];
        deltajump = deltajump(sortind);
        tdelay = [old_data.tdelay, tdelay];
        tdelay = tdelay(sortind);
        v_jump = [old_data.v_jump;v_jump];
        v_jump = v_jump(sortind,:);
    catch
        disp("First run detected. Creating new file...")
    end
else
    
    filename = "jumpdata_L"+num2str(Lstart);
    load(filename);
end
deltajump = deltajump(~isnan(tdelay));
mcol = mcol(~isnan(tdelay));
tdelay(isnan(tdelay)) = [];
f100 = figure(100);
subplot(1,3,1)
loglog(mcol, tdelay, '--o')
grid on
hold on
plot(mcol, 5*mcol.^(-2/3), '--')
xlabel("bif para rate of change, $m$", "Interpreter","latex")

subplot(1,3,2)
semilogx(mcol, deltajump, '--o')
grid on
hold on
plot(mcol, deltastab - 5*mcol.^(1/3), '--')
xlabel("bif para rate of change, $m$", "Interpreter","latex")
ylabel("bif para at jump")

f10 = figure(10);
loglog(deltastab - deltajump, tdelay, '--o')
grid on
hold on

if Lstart == 15
    tdelaycurrent = tdelay;
    load('sudden_change_delta_L15.mat')
    loglog(deltastab - deltatarget, tdelay, '--o')
    tdelay = tdelaycurrent;
    leg = legend(["Constant rate", "Instantaneous"], "Location", "southwest");
    leg.AutoUpdate = "off";
end
refvec = [min([deltastab-deltajump,deltastab - deltatarget]),max([deltastab-deltajump,deltastab - deltatarget])];
loglog(refvec,50000*refvec.^(-2), '--')
xlabel({"distance to", "stability boundary, $\delta_{\textrm{stab}} - \delta_{\textrm{change}}$"}, "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
% title("Constant rate of change")
set(f10,'Windowstyle','normal')
set(findall(f10,'-property','FontSize'),'FontSize',11)
set(f10,'Units','centimeters')
set(f10,'Position',[9 9 8 8])


set(f100,'Windowstyle','normal')
set(findall(f100,'-property','FontSize'),'FontSize',11)
set(f100,'Units','centimeters')
set(f100,'Position',[0 0 24 8])
%exportgraphics(f,"../../Figures/klausmeier_time_delay_L"+num2str(Lstart)+".eps", "Resolution", 500 )
Amdata = [mcol;deltajump;tdelay]; % save to dat file if used in batch calculation of spectra


save(filename, "mcol", "tdelay", "deltajump", "v_jump")



%% compare with spectrum data

cd ../Num_cont/spectra_calc_batch/spectrum_data/
Files = dir;

% extract files corresponding to correct wavelength
correctL = NaN*ones(1,length(Files));
for ff = 1:length(Files)
    correctL(ff) = contains(Files(ff).name, ["L_"+num2str(Lstart),"L_"+num2str(Lstart-1)+"dot99"]);
    filenames(ff) = string(Files(ff).name);
end
Files = Files(correctL==1);
filenames = filenames(correctL==1);
maxrespec = NaN*ones(1,length(deltajump));
faildelta = [];
for aa = 1:length(deltajump)
    loadind = find(contains(filenames,"delta_"+strrep(num2str(deltajump(aa)),".","dot")+"_and")); % find correct spectrum file
    try % some Atarget have no spectrum because no PTW with Atarget exists
        load(Files(loadind).name) %load spectrum data
        [maxrespec(aa),maxrealind] = max(spectrum_real); % extract max real part
    catch
        faildelta = [faildelta,deltajump(aa)];
    end
end

f11 = figure(11);
loglog(maxrespec, tdelay, '--o')
grid on
hold on

cd ../../../Num_sim/
if Lstart == 15
    tdelaycurrent = tdelay;
    maxrespeccurrent =  maxrespec;
    load('sudden_change_delta_L15.mat')
    loglog(maxrespec, tdelay, '--o')
    refvec = [min([maxrespeccurrent,maxrespec]),max([maxrespeccurrent,maxrespec])];
    loglog(refvec,5*refvec.^(-1), '--')
    tdelay = tdelaycurrent;
    maxrespec = maxrespeccurrent;
    leg = legend(["Constant rate", "Instantaneous"], "Location", "northeast");
    leg.AutoUpdate = "off";
end
xlabel("Max real part of spectrum", "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
% title("Constant rate of change")
set(f11,'Windowstyle','normal')
set(findall(f11,'-property','FontSize'),'FontSize',11)
set(f11,'Units','centimeters')
set(f11,'Position',[18 9 8 8])

f2 = figure(2);
semilogx(maxrespec, deltajump, '--o')
grid on
hold on
xlabel("Max real part of spectrum", "Interpreter","latex")
ylabel("bifurcation parameter, $\delta_{\textrm{target}}$", "Interpreter","latex")
title("Constant rate of change")
set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',11)
set(f2,'Units','centimeters')
set(f2,'Position',[27 9 8 8])

%% delay against A integral and compare with sudden change data

deltaint = tdelay.*(deltastab - mcol.*tdelay/2);
f3 = figure(3);
loglog(deltaint, tdelay, '--o')
grid on
hold on
if Lstart == 15
    % load data for sudden rate of change
    tdelaycurrent = tdelay;
    load('sudden_change_delta_L15.mat')
    deltaint_sudden = tdelay.*deltatarget;
    
    loglog(deltaint_sudden, tdelay, '--o')
    leg = legend(["Constant rate", "Instantaneous"], "Location", "southeast");
    leg.AutoUpdate = "off";
    tdelay = tdelaycurrent;
end
refvec = [min([deltaint,deltaint_sudden]), max([deltaint,deltaint_sudden])];
loglog(refvec, 0.006*refvec.^(1), "--")

xlabel({"Accumulated distance from stab.", "boundary $\overline{\delta}(t_{\textrm{stab}} + t_{\textrm{delay}})$"}, "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")



set(f3,'Windowstyle','normal')
set(findall(f3,'-property','FontSize'),'FontSize',11)
set(f3,'Units','centimeters')
set(f3,'Position',[0 9 8 8])

%% delay against spectrum integral

% extract all data held on spectra
maxrespec = NaN*ones(1,length(Files)); lambda_area = maxrespec;
cd ../Num_cont/spectra_calc_batch/spectrum_data/
try
for aa = 1:length(Files)
    loadind = aa; % find correct spectrum file
    filename = Files(aa).name;
    startind = strfind(filename,'_delta_'); startind = startind+7;
    endind = strfind(filename,'_and_L_'); endind = endind-1;
    deltaval_spec(aa) = str2num(strrep(filename(startind:endind),'dot','.')); %determine A value
    try % some Atarget have no spectrum because no PTW with Atarget exists
        load(Files(loadind).name) %load spectrum data
        [maxrespec(aa),maxrealind] = max(spectrum_real); % extract max real part
        tempind1 = find(spectrum_real>0); tempind2 = find(spectrum_imag>0); areaind = intersect(tempind1,tempind2);
        lambda_area(aa) = polyarea(spectrum_real(areaind),spectrum_imag(areaind));
    end
end
end
try
    figure(f2);
    hold on
    semilogx(maxrespec, deltaval_spec, '--o')
    
    [deltaval_spec,uniqueind] = unique(deltaval_spec);
    maxrespec = maxrespec(uniqueind);
    for ss = 1:length(tdelay)
        tvec = linspace(1,tdelay(ss),100);
        deltavec = deltastab - mcol(ss)*tvec;
        maxrespec_inter = interp1(deltaval_spec,maxrespec,deltavec);
        lambda_area_inter = interp1(deltaval_spec,lambda_area,deltavec);
        muint(ss) = trapz(tvec,maxrespec_inter);
        if isnan(muint(ss))
            1+1;
        end
        lambdaint(ss) = trapz(tvec,lambda_area_inter);
    end
end

try
    f4 = figure(4);
    scatter(muint,tdelay,'o')
    set(gca,'yscale','log')
    type = "constant rate";
    type = repmat(type,length(muint),1);
    datacomb = table(type,muint');
    datacomb.Properties.VariableNames = ["type","muint"];
    data1 = muint;
    grid on
    if Lstart == 15
        hold on
        load("muint_data_sudden_change")
        scatter(muint,tdelay,'o')
        legend(["Constant rate", "Instantaneous"], "Location", "southeast")

        type = "instantaneous";
        type = repmat(type,length(muint),1);
        datacomb1 = table(type,muint');
        datacomb1.Properties.VariableNames = ["type","muint"];

        datacomb = [datacomb;datacomb1];
        data2 = muint;
    end
    
      xlabel("Accumulated maximal instability, $\overline{\mu}$", "Interpreter","latex")
    ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
    % title("$L = "+num2str(Lstart)+"$", "Interpreter","latex")
    
    set(f4,'Windowstyle','normal')
    set(findall(f4,'-property','FontSize'),'FontSize',11)
    set(f4,'Units','centimeters')
    set(f4,'Position',[0 0 8 8])
end


