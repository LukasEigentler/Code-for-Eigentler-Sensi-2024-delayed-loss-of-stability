%% Klausmeier: Time delay simulations
% This script numerically calculates the time delay that occurs after
% crossing a stability boundary to the change in wavelength for different
% constant rates of change of the bifurcation parameter.
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023

% keep f4
clear; 
close all
plotonly = 1;
withpause = 0; % include a pause in the parameter change
withtermination = 0; % include termination of change
withchangem = 0; % include a sudden change in m to a fixed value for all realisations
withchangeA = 0; % include sudden change of bifurcation parameter at the stability boundary
if withpause*withchangem ~= 0 
    error("Wrong parameter selection: one of withpause or withchangem must be zero")
end
%% Parameters
Lstart = 20; L = Lstart;  %wavelength
B = 0.45; %Plant loss 
nu = 182.5; % advection speed
d = 500; %Water diffusion
P = 300; % Time after crossing stab boundary at which pause or change of m/A in parameter change occurs
Plength = 300; % duration of pause in parameter change
mbar = 1e-2;
calb_time = 100;
mcol = logspace(-6,-5,10); % rate of change A = As - mt
mcol = flip(mcol);

if Lstart == 20
    A_s = 1.7; A_e = 0.5; A = A_s;
    Astab = 1.69688; % bif para at crossing of boundary: not automated!!
    xmax = 200; % half of the space domain
elseif Lstart == 40
    A_s = 0.85; A_e = 0.4; A = A_s;
    Astab = 0.8382; % bif para at crossing of boundary: not automated!! 
    xmax = 200; % half of the space domain
elseif Lstart == 60
    A_s = 0.6; A_e = 0.35; A = A_s;
    Astab = 0.5787; % bif para at crossing of boundary: not automated!! 
    xmax = 240; % half of the space domain
end

%%  Space and Time grid

M = 2^9; % no of space points used in discretisation
x=linspace(-xmax,xmax,M);
dx = x(2)-x(1);
options = odeset('Stats', 'on', 'MaxStep',10); % max step size needed to capture changing parameter over time!

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
    %% loop over m
    tdelay = NaN*ones(1,length(mcol)); Ajump = tdelay;
    v_jump = NaN*ones(length(mcol),2*M);
    for mm = 1:length(mcol)
        try
        disp("Step "+num2str(mm)+" of "+num2str(length(mcol)))
        m = mcol(mm);
        % tspan = [0,calb_time+(A_s-A_e)/m];
        tspan = [0,2e4];
        kt1 = linspace(1e-5,tspan(end),10000);
        A_vec1 =  A_s-m*kt1; %Rainfall
        kt0 = linspace(0,calb_time); A_vec0 = A_s*ones(1,length(kt0));
        kt=[kt0,calb_time+kt1]; A_vec = [A_vec0,A_vec1];

        if withpause == 1
            [~,boundaryind] = min(abs(A_vec - Astab));
            pauseind = find(kt<kt(boundaryind) + P); pauseind = pauseind(end);
            kt1 = kt(1:pauseind-1); A_vec1 = A_vec(1:pauseind-1);
            kt2 = linspace(kt(pauseind),kt(pauseind)+Plength, 1000); A_vec2 = A_vec(pauseind)*ones(1,1000);
            kt3 = kt(pauseind+1:end)+Plength; A_vec3 = A_vec(pauseind+1:end);
            kt = [kt1,kt2,kt3]; A_vec = [A_vec1,A_vec2, A_vec3];
        elseif withtermination == 1
            tspan = [0,1e4];
            [~,boundaryind] = min(abs(A_vec - Astab));
            pauseind = find(kt<kt(boundaryind) + P); pauseind = pauseind(end);
            A_vec1 = A_vec(1:pauseind-1);
            A_vec2 = A_vec1(end)*ones(1,length(A_vec(pauseind:end)));
            A_vec = [A_vec1,A_vec2];
        elseif withchangem == 1
            tspan = [0,(A_s-A_e)/min(m,mbar)];
            kt = linspace(0,tspan(end),10000);
            A_vec =  A_s-m*kt; %Rainfall
            [~,boundaryind] = min(abs(A_vec - Astab));
            pauseind = find(kt<kt(boundaryind) + P); pauseind = pauseind(end);
            A_vec1 = A_vec(1:pauseind-1);
            A_vec2 = A_vec(pauseind) + mbar*kt(pauseind) - mbar*kt(pauseind:end);
            A_vec = [A_vec1,A_vec2];
        elseif withchangeA == 1
            [~,boundaryind] = min(abs(A_vec - Astab));
            pauseind = find(kt<kt(boundaryind) + P); pauseind = pauseind(end);
            A_vec1 = A_vec(1:pauseind-1);
            A_vec2 = 0.8*ones(1,length(A_vec(pauseind:end)));
            A_vec = [A_vec1,A_vec2];
        end
    
    
        [t_out,v] = ode15s(@(t,v) klausmeierode(t,v,B,nu,d,M,dx,kt,A_vec), tspan, vin, options);
    
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
        tdelay(mm) = tjump - tstab;
        Ajumpind = find(kt>tjump); Ajumpind = Ajumpind(1);
        Ajump(mm) = A_vec(Ajumpind);

     % find solution at jump
        v_jump(mm,:) = v(tjumpind-1,:);
        catch
            warning("Error occurred")
        end
    end
    %%
    if withpause == 1
        filename = "jumpdata_withpause_L"+num2str(Lstart);
    elseif withtermination == 1
        filename = "jumpdata_withtermination_L"+num2str(Lstart);
    elseif withchangem == 1
        filename = "jumpdata_withchangem_L"+num2str(Lstart);
    elseif withchangeA == 1
        filename = "jumpdata_withchangeA_L"+num2str(Lstart);
    else
        filename = "jumpdata_L"+num2str(Lstart);
    end
    try
        old_data = load(filename);
        mcol = [old_data.mcol,mcol];
        [mcol, sortind] = sort(mcol);
        Ajump = [old_data.Ajump, Ajump];
        Ajump = Ajump(sortind);
        tdelay = [old_data.tdelay, tdelay];
        tdelay = tdelay(sortind);
        v_jump = [old_data.v_jump;v_jump];
        v_jump = v_jump(sortind,:);
    catch
        disp("First run detected. Creating new file...")
    end
else
    if withpause == 1
        filename = "jumpdata_withpause_L"+num2str(Lstart);
    elseif withtermination == 1
        filename = "jumpdata_withtermination_L"+num2str(Lstart);
    elseif withchangem == 1
        filename = "jumpdata_withchangem_L"+num2str(Lstart);
    elseif withchangeA == 1
        filename = "jumpdata_withchangeA_L"+num2str(Lstart);
    else
        filename = "jumpdata_L"+num2str(Lstart);
    end
    load(filename);
end
Ajump = Ajump(~isnan(tdelay));
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
plot(mcol, Ajump, '--o')
grid on
hold on
plot(mcol, Astab - 5*mcol.^(1/3), '--')
xlabel("bif para rate of change, $m$", "Interpreter","latex")
ylabel("bif para at jump")

f10 = figure(10);
loglog(Astab - Ajump, tdelay, '--o')
grid on
hold on


if Lstart == 20
    tdelaycurrent = tdelay;
    load('sudden_change_A_L20.mat')
    loglog(Astab - Atarget, tdelay, '--o')
    tdelay = tdelaycurrent;
    leg = legend(["Constant rate", "Instantaneous"], "Location", "southwest");
    leg.AutoUpdate = "off";
end
refvec = [min([Astab - Atarget,Astab - Ajump]),max([Astab - Atarget,Astab - Ajump])];
loglog(refvec,30*refvec.^(-2), '--')
xlabel({"distance to", "stability boundary, $A_{\textrm{stab}} - A_{\textrm{change}}$"}, "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
%title("Constant rate of change")
set(f10,'Windowstyle','normal')
set(findall(f10,'-property','FontSize'),'FontSize',11)
set(f10,'Units','centimeters')
set(f10,'Position',[9 9 8 8])


set(f100,'Windowstyle','normal')
set(findall(f100,'-property','FontSize'),'FontSize',11)
set(f100,'Units','centimeters')
set(f100,'Position',[0 0 24 8])
%exportgraphics(f,"../../Figures/klausmeier_time_delay_L"+num2str(Lstart)+".eps", "Resolution", 500 )
Amdata = [mcol;Ajump;tdelay]; % save to dat file if used in batch calculation of spectra


save(filename, "mcol", "tdelay", "Ajump", "v_jump")

%% plot both case with pause and case with no pause on one plot
try
    
    filename1 = "jumpdata_L"+num2str(Lstart);
    filename2 = "jumpdata_withpause_L"+num2str(Lstart);
    filename3 = "jumpdata_withtermination_L"+num2str(Lstart);
    filename4 = "jumpdata_withchangem_L"+num2str(Lstart);
    filename5 = "jumpdata_withchangeA_L"+num2str(Lstart);

    load(filename1);
    Ajump1 = Ajump(~isnan(tdelay));
    mcol1 = mcol(~isnan(tdelay));
    tdelay(isnan(tdelay)) = []; tdelay1 = tdelay;

    load(filename2);
    Ajump2 = Ajump(~isnan(tdelay));
    mcol2 = mcol(~isnan(tdelay));
    tdelay(isnan(tdelay)) = []; tdelay2 = tdelay;

    load(filename3);
    Ajump3 = Ajump(~isnan(tdelay));
    mcol3 = mcol(~isnan(tdelay));
    tdelay(isnan(tdelay)) = []; tdelay3 = tdelay;

    load(filename4);
    Ajump4 = Ajump(~isnan(tdelay));
    mcol4 = mcol(~isnan(tdelay));
    tdelay(isnan(tdelay)) = []; tdelay4 = tdelay;

    load(filename5);
    Ajump5 = Ajump(~isnan(tdelay));
    mcol5 = mcol(~isnan(tdelay));
    tdelay(isnan(tdelay)) = []; tdelay5 = tdelay;

    f100 = figure;
    subplot(1,2,1)
    p(1) = loglog(mcol1, tdelay1, '--o', "DisplayName", "Case 1");
    grid on
    hold on
    p(2) = loglog(mcol2, tdelay2, '--o', "DisplayName", "Case 2");
    p(3) = loglog(mcol3, tdelay3, '--o', "DisplayName", "Case 3");
    p(4) = loglog(mcol4, tdelay4, '--o', "DisplayName", "Case 4");
    p(5) = loglog(mcol5, tdelay5, '--o', "DisplayName", "Case 5");
    loglog(mcol5,P*ones(1,length(mcol5)),'--')
    xlabel("bif para rate of change, $m$", "Interpreter","latex")
    ylabel("Time delay after crossing stab boundary")
    legend(p, "location", "southwest")
    xlim([-inf,1e-2])

    subplot(1,2,2)
    loglog(mcol1, Ajump1, '--o')
    grid on
    hold on
    semilogx(mcol2, Ajump2, '--o')
    semilogx(mcol3, Ajump3, '--o')
    semilogx(mcol4, Ajump4, '--o')
    semilogx(mcol5, Ajump5, '--o')
    xlabel("bif para rate of change, $m$", "Interpreter","latex")
    ylabel("bif para at jump")
    xlim([-inf,1e-2])
    set(f100,'Windowstyle','normal')
    set(findall(f100,'-property','FontSize'),'FontSize',11)
    set(f100,'Units','centimeters')
    set(f100,'Position',[0 0 20 9])
end


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
maxrespec = NaN*ones(1,length(Ajump1));
Afail = [];
for aa = 1:length(Ajump1)
    loadind = find(contains(filenames,"A_"+strrep(num2str(Ajump1(aa)),".","dot")+"_and")); % find correct spectrum file
    try % some Atarget have no spectrum because no PTW with Atarget exists
        load(Files(loadind).name) %load spectrum data
        [maxrespec(aa),maxrealind] = max(spectrum_real); % extract max real part
    catch
        Afail = [Afail,Ajump1(aa)];
    end
end

f11 = figure(11);
loglog(maxrespec, tdelay1, '--o')
grid on
hold on
cd ../../../Num_sim/
if Lstart == 20
    tdelaycurrent = tdelay;
    maxrespeccurrent =  maxrespec;
    load('sudden_change_A_L20.mat')
    loglog(maxrespec, tdelay, '--o')
    refvec = [min([maxrespeccurrent,maxrespec]),max([maxrespeccurrent,maxrespec])];
    loglog(refvec,10*refvec.^(-1), '--')
    tdelay = tdelaycurrent;
    maxrespec = maxrespeccurrent;
    leg = legend(["Constant rate", "Instantaneous"], "Location", "southwest");
    leg.AutoUpdate = "off";
end

xlabel("Max real part of spectrum", "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")
%title("Constant rate of change")
set(f11,'Windowstyle','normal')
set(findall(f11,'-property','FontSize'),'FontSize',11)
set(f11,'Units','centimeters')
set(f11,'Position',[18 9 8 8])

f2 = figure(2);
semilogx(maxrespec, Ajump1, '--o')
grid on
hold on
xlabel("Max real part of spectrum", "Interpreter","latex")
ylabel("bifurcation parameter, $A_{target}$", "Interpreter","latex")
title("Constant rate of change")
set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',11)
set(f2,'Units','centimeters')
set(f2,'Position',[27 9 8 8])

%% delay against A integral and compare with sudden change data

Aint = tdelay1.*(Astab - mcol1.*tdelay1/2);
f3 = figure(3);
loglog(Aint, tdelay1, '--o')
grid on
hold on

if Lstart == 20
    % load data for sudden rate of change
    
    load('sudden_change_A_L20.mat')
    Aint_sudden = tdelay.*Atarget;
    
    loglog(Aint_sudden, tdelay, '--o')
    leg = legend(["Constant rate", "Instantaneous"], "Location", "southeast");
    leg.AutoUpdate = "off";
end
refvec = [min([Aint_sudden,Aint]),max([Aint_sudden,Aint])];
loglog(refvec, 1*refvec.^(1), "--")

xlabel({"Accumulated distance from stab.", "boundary, $\overline{A}(t_{\textrm{stab}} + t_{\textrm{delay}})$"}, "Interpreter","latex")
ylabel("Time delay, $t_{\textrm{delay}}$", "Interpreter","latex")



set(f3,'Windowstyle','normal')
set(findall(f3,'-property','FontSize'),'FontSize',11)
set(f3,'Units','centimeters')
set(f3,'Position',[0 9 8 8])

%% delay against spectrum integral

% extract all data held on spectra
maxrespec = NaN*ones(1,length(Files)); lambda_area = maxrespec;
cd ../Num_cont/spectra_calc_batch/spectrum_data/
for aa = 1:length(Files)
    loadind = aa; % find correct spectrum file
    filename = Files(aa).name;
    startind = strfind(filename,'_A_'); startind = startind+3;
    endind = strfind(filename,'_and_L_'); endind = endind-1;
    Aval_spec(aa) = str2num(strrep(filename(startind:endind),'dot','.')); %determine A value
    try % some Atarget have no spectrum because no PTW with Atarget exists
        load(Files(loadind).name) %load spectrum data
        [maxrespec(aa),maxrealind] = max(spectrum_real); % extract max real part
        tempind1 = find(spectrum_real>0); tempind2 = find(spectrum_imag>0); areaind = intersect(tempind1,tempind2);
        lambda_area(aa) = polyarea(spectrum_real(areaind),spectrum_imag(areaind));
    end
end
try
    figure(f2);
    semilogx(maxrespec, Aval_spec, '--o')
    
    [Aval_spec,uniqueind] = unique(Aval_spec);
    maxrespec = maxrespec(uniqueind);
    for ss = length(tdelay1):-1:1
        tvec = linspace(1e-4,tdelay1(ss),100);
        Avec = Astab - mcol1(ss)*tvec;
        maxrespec_inter = interp1(Aval_spec,maxrespec,Avec);
        lambda_area_inter = interp1(Aval_spec,lambda_area,Avec);
        muint(ss) = trapz(tvec,maxrespec_inter);
        lambdaint(ss) = trapz(tvec,lambda_area_inter);
    end
end

try
    f4 = figure(4);
    scatter(muint,tdelay1,'o')
    set(gca,'yscale','log')
    type = "constant rate";
    type = repmat(type,length(muint),1);
    datacomb = table(type,muint');
    datacomb.Properties.VariableNames = ["type","muint"];
    data1 = muint;
    grid on
    if Lstart == 20
        hold on
        load("muint_data_sudden_change")
        scatter(muint,tdelay,'o')
        legend("Constant rate", "Instantaneous")
        
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

% f5 = figure(5);
% semilogy(lambdaint,tdelay1,'o')
% type = "constant rate";
% type = repmat(type,length(lambdaint),1);
% datacomb = table(type,lambdaint');
% datacomb.Properties.VariableNames = ["type","lambdaint"];
% data1 = lambdaint;
% grid on
% if Lstart == 20
%     hold on
%     load("lambdaint_data_sudden_change")
%     semilogy(lambdaint,tdelay,'o')
%     legend("Constant rate", "Instantaneous")
% 
%     type = "instantaneous";
%     type = repmat(type,length(lambdaint),1);
%     datacomb1 = table(type,lambdaint');
%     datacomb1.Properties.VariableNames = ["type","lambdaint"];
% 
%     datacomb = [datacomb;datacomb1];
%     data2 = lambdaint;
% end
% 
% xlabel("$\overline{\lambda}$", "Interpreter","latex")
% ylabel("Time delay, $t_{delay}$", "Interpreter","latex")
% 
% 
% set(f5,'Windowstyle','normal')
% set(findall(f5,'-property','FontSize'),'FontSize',11)
% set(f5,'Units','centimeters')
% set(f5,'Position',[0 0 8 8])

% %% mubar scatter plot
% f5 = figure;
% xdata = categorical(datacomb.type);
% scatter(xdata,datacomb.muint, 'o')
% grid on
% xlabel("Parameter change regime", "Interpreter", "latex")
% ylabel("$\overline{\mu}$", "Interpreter","latex")
% 
% set(f5,'Windowstyle','normal')
% set(findall(f5,'-property','FontSize'),'FontSize',11)
% set(f5,'Units','centimeters')
% set(f5,'Position',[0 0 8 8])
