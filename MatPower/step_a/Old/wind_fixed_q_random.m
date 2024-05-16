%Authors: Farley Rimon & Marouane Mastouri
%Last updated: 13/05/2020 (dd-mm-yyyy)

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
%branch_status_idx = find(mpc_substation_system13.branch(:,11)==0); %Find indexes of generating strings in generator data

%define all constants
%WTG Farm
% type_wtg= 4; %total types of WTG
  v_c_in = 3 ; %general cut in speed in m/s
  v_r = 14 ; %general rated speed in m/s, fixed for step 2a
  v_w = 14 ; %nominal wind speed in m/s
  v_c_off = 25 ; %general cut off speed in m/s
  global P_wt_max;
  P_wt_max = [33.6 29.4 29.4 16.8 32 29.4 16 28 33 29.15 28 28.8 28.8]; %Rated power in MW for each string
  global N_strings;
  N_strings = 13;

  %% Step 2: Generate random variables

%Each time, one of the substeps is uncommented and tested.
%Random variables are created with PDFs
N_daily_dispatch = 96; %total available wind dispatch profiles
N_strings = 13 %total Wind Turbine Generators strings

%Shunt random status, connected or not.


% %update status of branches for reactor
% % 
% % N_branches_shunt = 2; %Numbers of branches connected to the shunt
% % upd_status = zeros(N_branches_shunt,N_daily_dispatch); %initialize matrix that updates status
% % 
% % for(i=1:N_daily_dispatch)
% %     for(j= N_branches_shunt)
% %         
% %         upd_status(i,j) = randi([0,1]); %assign random 0 or 1 to branches for every dispatch
% %     end
% % end
% 

%Convert generated power at rated wind speed

P_wt = zeros(N_daily_dispatch, N_strings); %Make matrix for generated nominal wind speed


for(i=1: length(P_wt_max)) %makes a 1x13
    for(j=1: N_daily_dispatch) %makes a 96x1
    %%apply boundary conditions to determine which equation is valid
     if(v_w(j) <= v_c_in) 
         P_wt(j,i) =0; 

    elseif ( (v_w(j) > v_c_in) &&  (v_w(j) <= v_r) ) 
    
         P_wt(j,i)=P_wt_max(i)*(v_w(j)^3-v_c_in^3)/(v_r^3-v_c_in^3);
    
    elseif ((v_w(j) > v_r) && (v_w(j) <= v_c_off))
         P_wt(j,i) = P_wt_max(i); %assigns every power to every dispatch 
    
     elseif (v_w(j) > v_c_off)
         P_wt(j,i) = 0;
     end 
    end
end

%% Step 2a: Generate random Qs for fixed nominal wind speed (fixed tap position, neglect PV farm)

%Random Q_demands of Normal distribution for each string/generator

samples_Q = zeros(N_daily_dispatch,N_strings); %matrix for the random Q_demands
Q_wt_max = [21.2   18.55   19.15   10.9    22.4    19.15   12.2    19.6    22.4    13.248  19.6    19.6    19.6] ; %Max positive/negative reactive power generated for eacht string 1-13

Q_mean = zeros(1, N_strings); %matrix with mu for each string 1-13

samples_Q = zeros(length(Q_mean),1); %initialize matrix for Q random samples

for i=1:length(Q_mean)
    Q_std(i)=unifrnd(1,Q_wt_max(i)); %make a uniform distribution for the standard deviation string 1-13
end


for(i=1:N_strings) %makes a 1x13
  for(j=1: N_daily_dispatch) %makes a 96x1
    samples_Q(j,i)= normrnd(Q_mean(i), 0.75*Q_std(i), [1, 1]); %Creates random samples for Q in 96x13, some values violating the constraints
  end
end

constraints_neg = zeros(N_daily_dispatch,N_strings);
constraints_pos= zeros(N_daily_dispatch,N_strings);
   for(k=1: N_strings)
   constraints_neg = samples_Q(samples_Q < -Q_wt_max(k)); %Note values greater than positive max reactive 
   constraints_pos = samples_Q(samples_Q > Q_wt_max(k));  %Note values greater than negative max reactive
%     
    idx_constraints_neg(:,k) = (samples_Q(:,k) < -Q_wt_max(k)); %Note values greater than max reactive 
    idx_constraints_pos(:,k) = (samples_Q(:,k) >  Q_wt_max(k));
      
   end 
   
%% Step 3: Update case file and run Power Flow

stop =0;
%x=0;
while(stop == 0)
   for(i=1:N_daily_dispatch) 
   for (j=1:N_strings) 
   %for (i=1:N_daily_dispatch)
        %mpc_substation_system13.status(var_branch,11) = upd_branch(k) %%update status of the branch
       
        if(i < N_daily_dispatch)
        mpc_substation_system13.gen(P_wt_idx(j),3) = samples_Q(i,j);
        results_pf = runpf(mpc_substation_system13); %run power flow once Q is loaded per dispatch
        elseif(i == N_daily_dispatch)
            mpc_substation_system13.gen(P_wt_idx(j),3) = samples_Q(i,j);
            stop = 1;
% if(x == ((N_daily_dispatch * N_strings)-1000))
%     stop = 1; %will stop loading the Qs
 
   end
%%Put some graphs here to show behaviour for busses/branches of interest
    %end
end
   end
end


