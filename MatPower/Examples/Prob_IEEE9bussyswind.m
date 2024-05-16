%==========================================================================
%    Probabilistic Power flow calculation for the IEEE 9-bus test system
%==========================================================================
%By: Dr.ir. J.L. Rueda Torres
%    Assistant Professor - IEPG Section/ESE Department/EWI Faculty/TU Delft
%Date: 09-10-2017

close all
clear all
clc


%% Step 1: Load the data of the case
mpc_ieee9busWind = loadcase('caseIEEE9bussysWind');
pload_idx=find(mpc_ieee9busWind.bus(:,3)~=0);%Location/index of loads in bus matrix
pq_phi_ratio_load= mpc_ieee9busWind.bus(pload_idx,4)./mpc_ieee9busWind.bus(pload_idx,3); %Q/P of loads 

%Parameters of Weibull distribution of wind speed
k=2.02; %Shape parameter in p.u.
lambda=11; %Scale parameter in m/s
%Parameters of wind power plant
Pwpp=180; %Max. MW of the wind power plant
% Pwt=2; %Rated MW of individual wind generator.
wsin=3; % Cut-in wind speed in (m/s)
wsr=12; %Rated wind speed in (m/s)
wsout=20; % Cut-off wind speed in (m/s)
wpcoshphi=1; %Operation at unit power factor

%WPP data
%Max WPP P: 180  MW
%Min WPP P: 24 MW
%WPP Q (MVar): 25.6486 (cosphi=0.99) 36.5506 (cosphi=0.98) 45.1123 (cosphi=0.97) 59.1631 (cosphi=0.95)  87.1780 (cosphi=0.9)


%%Step2: Generate random variation of load and wind power
N=300; %Amount of random samples to be generated

rand('twister',5489); % Fix seeds of the random number generator

%Wind power
ws_rand = wblrnd(lambda,k,N,1); %Random wind speed following Weibull distribution
Pwpp_act=zeros(size(ws_rand));
for i=1:length(ws_rand) %Consideration of power curve
    if ws_rand(i)<wsin
        Pwpp_act(i)=0;
    elseif ws_rand(i)>wsin & ws_rand(i)<wsr
        Pwpp_act(i)=Pwpp*(ws_rand(i)^3-wsin^3)/(wsr^3-wsin^3);
    elseif ws_rand(i)>wsr & ws_rand(i)<wsout
        Pwpp_act(i)=Pwpp;
    elseif ws_rand(i)>wsout
        Pwpp_act(i)=0;
    end
end

% % histogram(Pwpp_act)

%Define mean values of nodal loads
load_mean=mpc_ieee9busWind.bus(pload_idx,3);

%Define standard deviations of nodal loads
load_std=zeros(length(load_mean),1);

for i=1:length(load_mean)
    load_std(i)=unifrnd(1,30);
end

%Create random load variations with normal distribution
samples_load=zeros(N,length(load_mean));
for i=1:length(load_mean)
    samples_load(:,i)=normrnd(load_mean(i),load_std(i),N,1);
end

% %Create random load variations with multivariate normal distribution
% cov_load=cov(samples_load);
% mmuu=mean(samples_load);
% total_loads=mvnrnd(mmuu,cov_load,N);
% 
% %To save random generation and load
% save dat_incert_ieee9busWind samples_load total_loads Pwpp_act

% load dat_incert_ieee9busWind

%% Step 3: Run probabilistic power flow calculation
count_flag=0;
printff=100;
feas_scn_idx=zeros(N,1);
timeOPF=zeros(N,1);

%Counters
iss=0;
issa=0;
issb=0;

% Set the options for optimal power flow calculation
% mpopt = mpoption('pf.alg','NR', 'verbose', 3, 'pf.tol', 1e-8, 'pf.enforce_q_lims', 2);
opt = mpoption('opf.ac.solver','MIPS', 'opf.flow_lim','S','opf.ignore_angle_lim',1,'verbose', 0, 'out.sys_sum', 0, 'out.area_sum', 0,'out.bus', 0,'out.branch', 0,'out.gen', 0,'out.lim.all',0);

while (count_flag==0)
    iss=iss+1;
    
    %Update the generated random values of nodal load (P &Q) and wind
    %generation (P)
    mpc_ieee9busWind.bus(pload_idx,3)=samples_load(iss,:)'; %Load active power
    mpc_ieee9busWind.bus(pload_idx,4)=samples_load(iss,:)'.*pq_phi_ratio_load; %Load reactive power

    %Update the generated random values of WPP active power output
    mpc_ieee9busWind.gen(4,2)=Pwpp_act(iss); %Active power output of WPP
    mpc_ieee9busWind.gen(4,4)=Pwpp_act(iss)*tan(acos(wpcoshphi)); %Qmax of WPP
    mpc_ieee9busWind.gen(4,5)=-Pwpp_act(iss)*tan(acos(wpcoshphi)); %Qmin of WPP
    mpc_ieee9busWind.gen(4,9)=Pwpp_act(iss); %Pmax - Fixed P dispatch of WPP
    mpc_ieee9busWind.gen(4,10)=Pwpp_act(iss); %Pmin -Fixed P dispatch of WPP
     
    %Run OPF (minimization of costs)
    tStart = cputime;
    OPFresults =opf(mpc_ieee9busWind,opt);
    timeOPF(iss)=cputime-tStart;
    
    if OPFresults.success==1 %if it converges
        issa=issa+1;    
        feas_scn_idx(iss)=OPFresults.success;
        
        %Store data of OPF output variables
        Pload(:,issa)=OPFresults.bus(:,3);
        Qload(:,issa)=OPFresults.bus(:,4);
        Pgen(:,issa)=OPFresults.gen(:,2);
        Qgen(:,issa)=OPFresults.gen(:,3);
        Vmag(:,issa)=OPFresults.bus(:,8);
        Vang(:,issa)=OPFresults.bus(:,9);
        Sact_from=sqrt(OPFresults.branch(:,14).^2+OPFresults.branch(:,15).^2);
        Sact_to=sqrt(OPFresults.branch(:,16).^2+OPFresults.branch(:,17).^2);
        Sfrom_act(:,issa)=Sact_from;
        Sto_act(:,issa)=Sact_to;
        Smax(:,issa)=max(Sact_from,Sact_to);
        Pbranch_loss(:,issa)=abs(OPFresults.branch(:,14)+OPFresults.branch(:,16));
        
        if issa>1 %Lines 135-142 are intented to illustrate how to compute statistical metrics of active power losses.
            issb=issb+1;
            Ploss_mean=mean(Pbranch_loss,2);
            Ploss_std=std(Pbranch_loss,0,2);
            [maxPlossmean,cc]=max(Ploss_mean);
            Conv_Ploss_mean(:,issb)=maxPlossmean;
            Conv_Ploss_error(:,issb)=Ploss_std(cc)/maxPlossmean;
        end
  
    
        if ((issa == 1) || (mod(issa,printff) == 0))
            fprintf('%7s    %7d  %7d\n','Iter:  ',issa, iss); 
        end
       
               
    end
    
    if (iss==N)
        count_flag=1;
    end
     
    
    
end

figure
plot(Vmag,'bo')
xlabel('Bus number')
ylabel('Voltage magnitude (p.u.)')


% boxplot(Vmag')
% boxplot(Smax','plotstyle','compact')

figure
plot(Smax,'bo')
hold on
plot(OPFresults.branch(:,6), 'rs')
xlabel('Branch number (from bus - to bus)')
ylabel('Apparent power (MVA)')
hold off

 
 
 
 
