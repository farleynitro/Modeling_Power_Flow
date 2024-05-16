%Authors: Farley Rimon & Marouane Mastouri
%Last updated: 04/06/2020 (dd-mm-yyyy)

%Modeling PF for a Substation with a Wind-and Solar Park
%A case study is loaded with 13 strings of Wind Turbine Generators and 4
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

for(i=1:length(P_wt_idx))
    
mpc_substation_system13.gen(P_wt_idx(i),8) = 0; %disconnect WTG to run only PV test

end 

%initialize normal active power generation
mpc_substation_system13.bus(1,3) = 6.31584; 
mpc_substation_system13.gen(2:5,2) = 1.57896;

B_limit_1_2 = 400;
B_limit_4_7 = 185.19;
B_limit_4_23 = 185.19;
B_limit_6_12 = 185.19;
B_limit_6_17 = 185.19;


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
  
%PV farm

  N_daily_dispatch = 300; %total available solar profiles
  N_strings_pv = 4; %total strings
  P_pv_max = [12.63168  12.63168  12.63168  12.63168];  %Avg active power in MW generated per PV string 1-4
  P_pv_min = [1.57896   1.57896   1.57896   1.57896];
  
  A_pv = 71478; %Area in m2 of each PV system , a total of 4 PV
  t_sun_h = 11; % Sun hours per day
  t_increase_h = 24/N_daily_dispatch; % 15minutes increase per day needed to go from 0 to max solar irradiance in W/m2
  t_steps = [1:t_increase_h:(6+t_increase_h)];
 
  %Insert season
  solar_irr_max= 800;                         %max solar irradiance in W/m2
  P_module_max = 288;                         %output per module for 800 W/m2
  rc_g_solar_power = (P_module_max)/((10-4)); %rc in W/m2 * 15min from 0 solar irradiance to max in, linear increase
  N_modules = 43860;                          %number of modules
  
  P_module_solar = rc_g_solar_power * t_steps;
  
 %% Step 2: Generate random variables

% %Each time, one of the substeps is uncommented and tested.
% %Random variables are created with PDFs
 Q_pv_max = (1/3)*P_pv_max; %Max reactive power in MVAr per PV string 1-4
%  P_pv = zeros(N_daily_dispatch, N_strings_pv); %Make matrix for generated nominal solar irradiance\
 
% %Calculate power generated at each PV String out of dispatch Solar
% %Irradiance (W/m2)
% 
% % %Summer model for solar
%    for(j=1: length(P_pv_max)) %makes a 1x4
%        k=1;
%          for(i=1: N_daily_dispatch) %makes a 1000x1
%             
%             if (i < (4/24)*N_daily_dispatch ) %Assume before 4 am , no PV energy generated
%                 P_pv(i,j) = 0;
%         
%             elseif( (i >= ((4/24)*N_daily_dispatch)) && (i < ((10/24)*N_daily_dispatch)) )%Assume 4 - 10 am , PV power increasing linearly
%                  P_pv(i,j) = (P_module_solar(k)*N_modules)*10^-6; %Output active power in MW
%                  if (P_pv(i,j) > P_pv_max(:,j))
%                      P_pv(i,j) = P_pv_max(:,j);
%                  end
%                  if (k<length(t_steps))
%                     k=k+1;
%                  end 
%             elseif( (i >= ((10/24)*N_daily_dispatch)) && (i < ((21/24)*N_daily_dispatch)) ) %Assume 10 am - 9pm , PV max energy generated
%                    P_pv(i,j) = P_pv_max(1,j);  
%          
%             elseif( (i >= ((21/24)*N_daily_dispatch))  ) %Assume 9 pm -12 am, no PV energy generated
%                  P_pv(i,j) = 0;
%          end
%     end
% end
%% Step 2b: Generate random Qs for fixed solar irradiance (and fixed tap
%position, neglect wind farm)

%1 Assume +-Qmax is 1/3 Pmax

%2 Assume Pmax, Q=0 & P=0, Qmax, extra task)
a= 0.08;

for(j=1:N_strings_pv)
for(i=1:(N_daily_dispatch))
    if(i <= (8*N_daily_dispatch/12))
        Q_mean1(i,j) = -4.2106;
    elseif ((i > (8*N_daily_dispatch/12)) && (i <=(12*N_daily_dispatch/12)))
    for(i=1:(4*N_daily_dispatch/12))
        Q_mean2(i,j) = 0;  
    end
    end
end
end

for(i=1:N_daily_dispatch)
    Q_pv_max_matrix(i,:) = Q_pv_max;
end

for(i=1:N_daily_dispatch)
    Q_std_pv(i,:) = unifrnd(1,Q_pv_max_matrix(i,:));
end

for(i=1:N_strings_pv) %create 4x1
for(j=1:(8*N_daily_dispatch/12)) %makes a 96x1
    samples_Q_pv(j,i)= normrnd(Q_mean1(i), a*Q_std_pv(i), [1, 1]); %Creates random samples for Q in 96x13, some values violating the constraints
  end
  for(j=(8*N_daily_dispatch/12):(12*N_daily_dispatch/12))
       samples_Q_pv(j,i)= normrnd(Q_mean2(i), a*Q_std_pv(i), [1, 1]); %Creates random samples for Q 
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

%% Step 3: Update case and run powerflow

stop = 0;
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
        [results_loss results_inject] = get_losses(results_pf); %losses 
        
        %Make matrix to specify branch losses and reactive power injection
        results_loss_matrix = zeros(31, 50); %first row is bus "from", second row bus "to", third row is loss

        bus_from = mpc_substation_system13.branch(:,1);
        bus_to = mpc_substation_system13.branch(:,2);
    
        %Perform operation only if system converges
          if (results_pf.success == 1)
              
        %Save data for each power flow  
        Vmag(:,i) = results_pf.bus(:,8); %Saves voltage magnitude for every bus
        Vangle(:,i) = results_pf.bus(:,9); %Saves angle magnitude for every bus
        P_flow_1_2(:,i) = results_pf.branch(1,16); %Saves Power flow in MW for Branch to PCC
        P_flow_4_23(:,i) = results_pf.branch(25,16); %Saves Power flow in MW for Branch to Input transformer 1
        P_flow_4_7(:,i) = results_pf.branch(6,16); %Saves Power flow in MW for Branch to Input transformer 1
        P_flow_6_12(:,i) = results_pf.branch(12,16); %Saves Power flow in MW for Branch to Input transformer 2
        P_flow_6_17(:,i) = results_pf.branch(18,16); %Saves Power flow in MW for Branch to Input transformer 2
        
        results_loss_real(:,i) = real(results_loss);
        results_loss_im(:,i) = imag(results_loss);
        results_inject_t(:,i) = results_inject;
        results_loss_matrix = [bus_from bus_to results_loss_real results_loss_im    results_inject_t];
        
          end 

          elseif(i == N_daily_dispatch)
            mpc_substation_system13.gen(P_pv_idx(j),3) = samples_Q_pv(i,j);
            %mpc_substation_system13.gen(P_pv_idx(j),2) = P_pv_max(j);
            mpc_substation_system13.gen(P_pv_idx(j),2) = P_pv_min(j);
            
         if (results_pf.success == 1)
             
        %Save data for each power flow 
        Vmag(:,i) = results_pf.bus(:,8);
        Vangle(:,i) = results_pf.bus(:,9);
        P_flow_1_2(:,i) = results_pf.branch(1,16); %Saves Power flow in MW for Branch to PCC
        P_flow_4_23(:,i) = results_pf.branch(25,16); %Saves Power flow in MW for Branch to Input transformer 1
        P_flow_4_7(:,i) = results_pf.branch(6,16); %Saves Power flow in MW for Branch to Input transformer 1
        P_flow_6_12(:,i) = results_pf.branch(12,16); %Saves Power flow in MW for Branch to Input transformer 2
        P_flow_6_17(:,i) = results_pf.branch(18,16); %Saves Power flow in MW for Branch to Input transformer 2
        
        %Save all the losses generated in each bus/branch for all
        %dispatches
        results_loss_real(:,i) = real(results_loss);
        results_loss_im(:,i) = imag(results_loss);
        results_inject_t(:,i) = results_inject;
        results_loss_matrix = [bus_from bus_to results_loss_real results_loss_im    results_inject_t];
        
        %Computes location and value of maximum loss
        %[p_max_loss p_max_loss_idx] = max(results_loss_real,[],2);
        %[q_max_loss q_max_loss_idx]=  max(results_loss_im,[],2);
        constraints_v_low = Vmag(Vmag < 0.8567); %Note values greater than minimum p.u.
        constraints_v_high = Vmag(Vmag > 1.14);  %Note values greater than maximum p.u.
     
        idx_constraints_v_low = zeros(28, N_daily_dispatch); %28 is amount of busses in the circuit
        
for (a=1:N_daily_dispatch)
        %Vlow is ..... for shunt connected
        %Vlow is ..... for shunt disconnected
        %Vhigh is ..... for shunt connected
        %Vhigh is ..... for shunt disconnected
    
        idx_constraints_v_low(:,a) = (Vmag(:,a) < 0.8567); 
        idx_constraints_v_high(:,a) = (Vmag(:,a) > 1.14);
        [idx_r_Q_v_low idx_c_Q_v_low] = find(idx_constraints_v_low~=0); %find location in matrix where there is v low
        [idx_r_Q_v_high idx_c_Q_v_high] = find(idx_constraints_v_high~=0);
end
        
        Q_dem_constraints_v_low = samples_Q_pv(idx_c_Q_v_low,:); %Find the Q demands for that specific dispatch
        Q_dem_constraints_v_high = samples_Q_pv(idx_c_Q_v_high,:); %Find the Q demands for that specific dispatch
        
        P_loss_v_low = results_loss_real(:, idx_c_Q_v_low);
        P_loss_v_high = results_loss_real(:, idx_c_Q_v_high);
        
        Q_loss_v_low = results_loss_im(:, idx_c_Q_v_low);
        Q_loss_v_high = results_loss_im(:, idx_c_Q_v_high);
        
        Q_inj_v_low = results_inject_t(:, idx_c_Q_v_low);
        Q_inj_v_high = results_inject_t(:, idx_c_Q_v_high);
        
         end
 end             
 end        
      stop = 1; %stop while loop
 end
 end
%end 


        %% Step 4: Plot results in a visable manner
        
        %Print  voltage data
        figure(1)
        subplot(2,1,1)
        plot(Vmag,'o');
        yline(0.8567,'r','LineWidth',1);
        yline(1.1424,'r','LineWidth',1);
        xlabel('Bus number','FontSize',12);
        ylabel('V magnitude (p.u.)','FontSize',12);
        title('Voltage magnitude in p.u. for bus numbers','FontSize',12)
        hold on;
        
        subplot(2,1,2)
        plot(Vangle,'s');
        xlabel('Bus number','FontSize',12);
        ylabel('V angle (degrees)','FontSize',12);
        title('Voltage angles in degrees for bus numbers','FontSize',12)
        hold on;
        
        %Print losses data
        figure(2)
        subplot(3,1,1)
        plot(results_loss_real,'o');
        xlabel('Branch #','FontSize',12);
        ylabel('P loss(MW)','FontSize',12);
        title('Active power loss in MW for branch','FontSize',12)
        hold on;

        
        subplot(3,1,2)
        plot(results_loss_im,'s');
        xlabel('Branch #','FontSize',12);
        ylabel('Q injected(MVAr)','FontSize',12);
        title('Reactive power inject in MVAr for branch','FontSize',12)
        hold on;
        
        subplot(3,1,3)
        plot(results_inject_t,'o');
        xlabel('Branch #','FontSize',12);
        ylabel('Q injected(MVAr)','FontSize',12);
        title('Reactive power inject in MVAr for branch','FontSize',12)
        hold on;
        
        figure(6)

        plot(P_flow_1_2,'o');
       xlabel('Dispatch','FontSize',12);
        ylabel('Power flow (MW)','FontSize',12);
        yline(B_limit_4_7 ,'r','LineWidth',1);
        title('Power flow to input Transformer 1 Branch 4-7','FontSize',12);
        hold on;
        
        figure(7)
        subplot(2,1,1)
     
        plot(P_flow_4_7,'o');
       xlabel('Dispatch','FontSize',12);
        ylabel('Power flow (MW)','FontSize',12);
        yline(B_limit_4_7 ,'r','LineWidth',1);
        title('Power flow to input Transformer 1 Branch 4-7','FontSize',12);
        hold on;
        
        subplot(2,1,2)
     
        plot(P_flow_4_23,'o');
        xlabel('Dispatch','FontSize',12);
        ylabel('Power flow (MW)','FontSize',12);
        yline(B_limit_4_23 ,'r','LineWidth',1);
        title('Power flow to input Transformer 1 Branch 4-23','FontSize',12)
        hold on;
        
        figure(8)
        subplot(2,1,1)
     
        plot(P_flow_6_12,'o');
        xlabel('Dispatch','FontSize',12);
        ylabel('Power flow (MW)','FontSize',12);
        yline(B_limit_6_12 ,'r','LineWidth',1);
        title('Power flow to input Transformer 2 Branch 6-12','FontSize',12)
        hold on;
        
        subplot(2,1,2)
     
        plot(P_flow_6_17,'o');
        xlabel('Dispatch','FontSize',12);
        ylabel('Power flow (MW)','FontSize',12);
        yline(B_limit_6_17 ,'r','LineWidth',1);
        title('Power flow to input Transformer 2 Branch 6-17','FontSize',12)
        hold on;
        
        %Print relation between demand, total system losses and voltage
        
        %Calculate average voltage for specific Q_demand
        v_avg = mean(Vmag);
        
        %Calculate average minimum voltage for specific Q_demand
        %this is equal to constraints_v_low
        
        %Calculate total losses for average voltage of each dispatch
        P_loss_tot = sum(results_loss_real);
        Q_loss_tot = sum(results_loss_im);
        Q_inj_tot = sum(results_inject_t);
        
        %Calculate total losses for voltage constraints
        P_loss_tot_v_low = sum(P_loss_v_low);
        Q_loss_tot_v_low = sum(Q_loss_v_low);
        Q_inj_tot_v_low = sum(Q_inj_v_low);
        
        P_loss_tot_v_high = sum(P_loss_v_high);
        Q_loss_tot_v_high = sum(Q_loss_v_high);
        Q_inj_tot_v_high = sum(Q_inj_v_high);
        
        %Calculate total Q_demand
        
        Q_tot = sum(samples_Q_pv,2);
        
        %Calculate total Q_demand for dispatches causing violation
       
        Q_tot_constraints_v_low = sum(Q_dem_constraints_v_low,2);
        Q_tot_constraints_v_high = sum(Q_dem_constraints_v_high,2);
         
        %Plot 3D plot Q demand, P loss and (minimum and maximum) average voltage
        figure(3)
       
        plot3(Q_tot, P_loss_tot, v_avg,'o','Color', 'r', 'MarkerSize',5); %plot for avg
        hold on;
        plot3(Q_tot_constraints_v_low, P_loss_tot_v_low,constraints_v_low,'o','Color', 'b', 'MarkerSize',5); %plot for v constraints low
        hold on;
        %plot3(Q_tot_constraints_v_high, P_loss_tot_v_high, constraints_v_high,'o','Color', 'k', 'MarkerSize',5); %plot for v constraints low
        hold on;
        
        grid on;
        
        
        xlabel('Q demand (MVAr)','FontSize',12);
        ylabel('P loss total (MW)','FontSize',12);
        zlabel('V average magnitude (p.u.)','FontSize',12);
        
        %Plot 3D plot Q loss for minimum and average voltage
        figure(4)
      
        plot3(Q_tot, Q_loss_tot, v_avg,'o','Color', 'r', 'MarkerSize',5); %plot for avg
        hold on;
        plot3(Q_tot_constraints_v_low, Q_loss_tot_v_low, constraints_v_low,'o','Color', 'b', 'MarkerSize',5); %plot for v constraints low
        hold on;
        %plot3(Q_tot_constraints_v_high, Q_loss_tot_v_high', constraints_v_high,'o','Color', 'k', 'MarkerSize',5); %plot for v constraints low
        hold on;
        
        grid on;
        xlabel('Q demand (MVAr)','FontSize',12);
        ylabel('Q loss total (MVAr)','FontSize',12);
        zlabel('V average magnitude (p.u.)','FontSize',12);
       
        %Plot 3D plot for minimum and average
        figure(5)
        
        plot3(Q_tot, Q_inj_tot,v_avg,'o','Color', 'r', 'MarkerSize',5); %plot for avg
        hold on;
        plot3(Q_tot_constraints_v_low,Q_inj_tot_v_low, constraints_v_low,'o','Color', 'b', 'MarkerSize',5); %plot for v constraints low
        hold on;
        %plot3(Q_tot_constraints_v_high, Q_inj_tot_v_high', constraints_v_high,'o','Color', 'k', 'MarkerSize',5); %plot for v constraints low
        hold on;
        
        grid on;
        
        xlabel('Q demand (MVAr)','FontSize',12);
        ylabel('Q INJ total (MVAr)','FontSize',12);
        zlabel('V average magnitude (p.u.)','FontSize',12);
