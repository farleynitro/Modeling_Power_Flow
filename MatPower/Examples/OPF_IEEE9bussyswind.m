%==========================================================================
%         Power flow calculation for the IEEE 9-bus test system
%==========================================================================
%By: Dr.ir. J.L. Rueda Torres
%    Assistant Professor - IEPG Section/ESE Department/EWI Faculty/TU Delft
%Date: 09-10-2017

close all
clear all
clc


%% Step 1: Load the data of the case
mpc_ieee9busWind = loadcase('caseIEEE9bussysWind');
Pwpp_act=180; %Define dispatch (MW) of the wind power plant (WPP)
wpcoshphi=1; %Operation of WPP with power factor between 0.95 overexcited and 0.95 underexcited
%WPP data
%Max WPP P: 180 MW
%Min WPP P: 24 MW
%WPP Q (MVar): 25.6486 (cosphi=0.99) 36.5506 (cosphi=0.98) 45.1123 (cosphi=0.97) 59.1631 (cosphi=0.95)  87.1780 (cosphi=0.9)

%Update the generated random values of WPP active power output
mpc_ieee9busWind.gen(4,2)=Pwpp_act; %Active power output of WPP
mpc_ieee9busWind.gen(4,4)=Pwpp_act*tan(acos(wpcoshphi)); %Qmax of WPP
mpc_ieee9busWind.gen(4,5)=Pwpp_act*tan(acos(wpcoshphi)); %Qmin of WPP
mpc_ieee9busWind.gen(4,9)=Pwpp_act; %Pmax - Fixed P dispatch of WPP
mpc_ieee9busWind.gen(4,10)=Pwpp_act; %Pmin -Fixed P dispatch of WPP

%% Step 2: Set the options for optimal power flow calculation
% mpopt = mpoption('pf.alg','NR', 'verbose', 3, 'pf.tol', 1e-8, 'pf.enforce_q_lims', 2);
opt = mpoption('opf.ac.solver','MIPS', 'opf.flow_lim','S','verbose', 3,'out.lim.all',2);% 'out.sys_sum', 0, 'out.area_sum', 0,'out.bus', 0,'out.branch', 0,'out.gen', 0,'out.lim.all',0);

%% Step 3: Run power flow calculation
% PFresults = runpf(mpc_ieee9busWind, mpopt);
% OPFresults =runopf(mpc_ieee9busWind);
OPFresults =opf(mpc_ieee9busWind,opt);

%Max WPP P: 90  MW
%Min WPP P: 23.8893 MW
%WPPQ: 25.6486 (cosphi=0.99) 36.5506 (cosphi=0.98) 45.1123 (cosphi=0.97) 59.1631 (cosphi=0.95)  87.1780 (cosphi=0.9)

%% Plots
figure
bar(OPFresults.bus(:,8));  %Bus voltage magnitudes
xlabel('Bus number')
ylabel('Voltage magnitude (p.u.)')

figure
bar(abs(OPFresults.branch(:,14)+OPFresults.branch(:,16)));  %Active power branch losses
xlabel('Branch number')
ylabel('Ploss (MW)')

figure
stem(OPFresults.branch(:,6), 'rs')
hold on
Sact_from=sqrt(OPFresults.branch(:,14).^2+OPFresults.branch(:,15).^2);
Sact_to=sqrt(OPFresults.branch(:,16).^2+OPFresults.branch(:,17).^2);
stem(max(Sact_from,Sact_to)); %Distance to thermal limit per branch
xlabel('Branch number (from bus - to bus)')
ylabel('Apparent power (MVA)')
hold off
% legend('MVA rating', 'Actual MVA')

% load('ConvInfo.mat')
% Vm=abs(Vconv);
% Vang=angle(Vconv)*180/pi;
 
% figure
% plot(Vang(2,:),Fconv(1,:),'bs')
% plot(Vang(3,:),Fconv(2,:),'bo')
% plot(Vm(end,:),Fconv(end,:),'bo')
 
 
 
 
