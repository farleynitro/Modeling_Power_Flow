%Authors: Farley Rimon & Marouane Mastouri
%Last updated: 18/05/2020 (dd-mm-yyyy)

%Modeling PF for a Substation with a Wind-and Solar Park
%A case study is loaded with 13 strings of Wind Turbine Generators and 3
%strings from a PV-module.
%Step 2 involves the modeling of the power flow in situations. These are
%useful to determine limiting constraints in the substation's power flow

clear all;
close all;
clc;

%% Step 1: Load data of the case
mpc_substation_system13 = loadcase('system_13');
P_idx = find(mpc_substation_system13.gen(:,2)~=0); %Find indexes of generating strings in generator data
P_wt_idx = P_idx(5:end);
P_pv_idx = P_idx(1:4);

for(i=1:length(P_pv_idx))
    
mpc_substation_system13.gen(P_wt_idx(i),8) = 0; %disconnect WTG to run only PV test

end 

mpc_substation_system13.bus(8,2) = 4;
mpc_substation_system13.bus(9,2) = 4;
mpc_substation_system13.bus(10,2) = 4;

mpc_substation_system13.bus(13,2) = 4;
mpc_substation_system13.bus(14,2) = 4;
mpc_substation_system13.bus(15,2) = 4;

mpc_substation_system13.bus(18,2) = 4;
mpc_substation_system13.bus(19,2) = 4;
mpc_substation_system13.bus(20,2) = 4;
mpc_substation_system13.bus(21,2) = 4;

mpc_substation_system13.bus(24,2) = 4;
mpc_substation_system13.bus(25,2) = 4;
mpc_substation_system13.bus(26,2) = 4;

mpc_substation_system13.branch(7,11) = 0;
mpc_substation_system13.branch(8,11) = 0;
mpc_substation_system13.branch(9,11) = 0;

mpc_substation_system13.branch(13,11) = 0;
mpc_substation_system13.branch(14,11) = 0;
mpc_substation_system13.branch(15,11) = 0;

mpc_substation_system13.branch(19,11) = 0;
mpc_substation_system13.branch(20,11) = 0;
mpc_substation_system13.branch(21,11) = 0;
mpc_substation_system13.branch(22,11) = 0;

mpc_substation_system13.branch(26,11) = 0;
mpc_substation_system13.branch(27,11) = 0;
mpc_substation_system13.branch(28,11) = 0;


%branch_status_idx = find(mpc_substation_system13.branch(:,11)==0); %Find indexes of generating strings in generator data

%define all constants
  N_strings = 13;
  
%PV farm

  N_daily_dispatch = 16; %total available wind profiles
  N_strings_pv = 4; %total strings
  P_pv_max = [30.3   30.3  30.3  30.3];  %Max active power in MW generated per PV string 1-4 
 
  A_pv = 181196.1722; %Area in m2 of each PV system , a total of 4 PV
  t_sun_h = 11; % Sun hours per day
  t_increase_h = 6; % 15minutes increase per day needed to go from 0 to max solar irradiance in W/m2
  t_steps = [1:0.25:6];
 
  %Insert seasons
  solar_irr_max= 800; %max solar irradiance in W/m2
  n_pv = 0.22; %efficiency of pv panel
  n_inv =0.95; %efficiency of inverter
  rc_g_solar_irr = (solar_irr_max)/((10-4)); %rc in W/m2 * 15min from 0 solar irradiance to max in, linear increase
  
 
  solar_irr = rc_g_solar_irr * t_steps;
  
  %% Step 2: Generate random variables

% %Each time, one of the substeps is uncommented and tested.
% %Random variables are created with PDFs
 Q_pv_max = (1/3)*P_pv_max; %Max reactive power in MVAr per PV string 1-4
 P_pv = zeros(N_daily_dispatch, N_strings_pv); %Make matrix for generated nominal solar irradiance\
 
% %Calculate power generated at each PV String out of dispatch Solar
% %Irradiance (W/m2)
% 
% %Summer model for solar
   for(i=1: length(P_pv_max)) %makes a 1x4
       k=1;
         for(j=1: N_daily_dispatch) %makes a 96x1
            
            if (j < (4/24)*N_daily_dispatch ) %Assume before 4 am , no PV energy generated
                P_pv(j,i) = 0;
        
            elseif( (j >= ((4/24)*N_daily_dispatch)) && (j < ((10/24)*N_daily_dispatch)) )%Assume 4 - 10 am , PV power increasing linearly
             
                 P_pv(j,i) = ((solar_irr(k) * A_pv)*n_inv*n_pv)*10^-6; %Output active power in MW
                 if (k<length(t_steps))
                    k=k+1;
                 end 
            elseif( (j >= ((10/24)*N_daily_dispatch)) && (j < ((21/24)*N_daily_dispatch)) ) %Assume 10 am - 9pm , PV max energy generated
                   P_pv(j,i) = P_pv_max(1,i);  
         
            elseif( (j >= ((21/24)*N_daily_dispatch))  ) %Assume 9 pm -12 am, no PV energy generated
                 P_pv(j,i) = 0;
         end
    end
end

% Step 2b: Generate random Qs for fixed solar irradiance (and fixed tap
%position, neglect wind farm)

%1 Assume +-Qmax is 1/3 Pmax

%2 Assume Pmax, Q=0 & P=0, Qmax, extra task)

samples_Q_pv = zeros(N_daily_dispatch,N_strings_pv); %matrix for the random Q_demands

Q_mean_pv = zeros(1, N_strings_pv); %matrix with mu for each string 1-4

samples_Q_pv = zeros(length(Q_mean_pv),1); %initialize matrix for Q random samples

for i=1:length(Q_mean_pv)
    Q_std_pv(i)=unifrnd(1,Q_pv_max(i)); %make a uniform distribution for the standard deviation string 1-4
end


for(i=1:N_strings_pv) %makes a 1x4
  for(j=1: N_daily_dispatch) %makes a 96x1
    samples_Q_pv(j,i)= normrnd(Q_mean_pv(i), 1*Q_std_pv(i), [1, 1]); %Creates random samples for Q in 96x4, some values violating the constraints
  end
end

constraints_neg_pv = zeros(N_daily_dispatch,N_strings_pv);
constraints_pos_pv = zeros(N_daily_dispatch,N_strings_pv);
   for(k=1: N_strings_pv)
   constraints_neg_pv = samples_Q_pv(samples_Q_pv < -Q_pv_max(k)); %Note values greater than positive max reactive 
   constraints_pos_pv = samples_Q_pv(samples_Q_pv > Q_pv_max(k));  %Note values greater than negative max reactive
     
    idx_constraints_neg_pv(:,k) = (samples_Q_pv(:,k) < -Q_pv_max(k)); %Note values greater than max reactive 
    idx_constraints_pos_pv(:,k) = (samples_Q_pv(:,k) >  Q_pv_max(k));
      
   end 
   
%% Step 3: Update case file and run Power Flow
stop =0;
printff = 4; %Prints values every 4 iterations

while(stop == 0)
   for(i=1:N_daily_dispatch) 
   for (j=1:N_strings_pv)
   %for (k=1:N_branches_shunt)
       
          if(i < N_daily_dispatch)
        mpc_substation_system13.gen(P_pv_idx(j),3) = samples_Q_pv(i,j);
        %mpc_substation_system13.bus(28,2) =   upd_status_bus(k,i) ; %Find 
        %mpc_substation_system13.branch(branch_r_status_idx(k),11) = upd_status(k,i); %Find indexes of generating strings in generator data
        
        %run power flow for each dispatch
        tStart = cputime;
        results_pf = runpf(mpc_substation_system13); %run power flow once Q is loaded per dispatch
        results_loss = get_losses(results_pf); %losses 
        
        %Make matrix to specify branch and its loss
        results_loss_matrix = zeros(31, 4); %first row is bus "from", second row bus "to", third row is loss

        bus_from = mpc_substation_system13.branch(:,1);
        bus_to = mpc_substation_system13.branch(:,2);
    
        results_loss_real = real(results_loss);
        results_loss_im = imag(results_loss);
        results_loss_matrix= [bus_from bus_to results_loss_real results_loss_im];
        
        
        %Perform operation only if system converges
          if (results_pf.success == 1)
              
        %Save data for each power flow  
        Vmag(:,i) = results_pf.bus(:,8); %Saves voltage magnitude for every bus
        Vangle(:,i) = results_pf.bus(:,9); %Saves angle magnitude for every bus

          end 

          elseif(i == N_daily_dispatch)
            mpc_substation_system13.gen(P_pv_idx(j),3) = samples_Q_pv(i,j);
            
         if (results_pf.success == 1)
        %Save data for each power flow 
        Vmag(:,i) = results_pf.bus(:,8);
        Vangle(:,i) = results_pf.bus(:,9);    
        
        constraints_v_low = Vmag(Vmag < 0.85); %Note values greater than minimum p.u.
        constraints_v_high = Vmag(Vmag > 1.05);  %Note values greater than maximum p.u.
%     

idx_constraints_v_low = zeros(28, N_daily_dispatch); %28 is amount of busses in the circuit
for (a=1:N_daily_dispatch)
        idx_constraints_v_low(:,a) = (Vmag(:,a) < 0.85); %Note values greater than max reactive 
        idx_constraints_v_high(:,a) = (Vmag(:,a) > 1.05);
end
        
        %Print  voltage data
        figure(1)
        subplot(2,1,1)
        plot(Vmag,'o');
        xlabel('bus number');
        ylabel('V magnitude (p.u.)');
        hold on;
        
        subplot(2,1,2)
        plot(Vangle,'s');
        xlabel('bus number');
        ylabel('V angle (degrees)');
        hold on;
        
        %Print losses data
        figure(2)
        subplot(2,1,1)
        plot(results_loss_real,'o');
        xlabel('Branch');
        ylabel('Real power loss(MW)');

        
        subplot(2,1,2)
        plot(results_loss_im,'o');
        xlabel('Branch');
        ylabel('Reactive power loss(MVAr)');
        hold on;
   
                     end
                 end             
             end
              
      stop = 1; %stop while loop
            end
        end
%end 
