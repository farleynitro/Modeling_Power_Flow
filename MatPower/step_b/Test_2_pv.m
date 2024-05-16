
close all;
clear all;
clc

%PV farm
 N_daily_dispatch = 96; %total available wind profiles
 N_strings_pv = 4; %total strings
  P_pv_max = [30.296   30.296  30.296  30.296];  %Max active power in MW generated per PV string 1-4 
  A_pv = 181196.1722; %Area in m2 of each PV system , a total of 4 PV
  t_sun_h = 11; % Sun hours per day
  t_increase_h = 6; % 15minutes increase per day needed to go from 0 to max solar irradiance in W/m2
  t_steps = [1:0.25:6];
  Q_pv_max = (1/3)*P_pv_max; %Max reactive power in MVAr per PV string 1-4
  %Insert seasons
  solar_irr_max= 800; %max solar irradiance in W/m2
  n_pv = 0.22; %efficiency of pv panel
  n_inv =0.95; %efficiency of inverter
  rc_g_solar_irr = (solar_irr_max)/((10-4)); %rc in W/m2 * 15min from 0 solar irradiance to max in, linear increase
  
  P_pv = zeros(N_daily_dispatch, N_strings_pv); %Make matrix for generated nominal solar irradiance\

  solar_irr = rc_g_solar_irr * t_steps;
 
  
% 
 

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


%1 Assume +-Qmax is 1/3 Pmax

%(2 Assume Pmax, Q=0 & P=0, Qmax, extra task)

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