%Authors: Farley Rimon & Marouane Mastouri
%Last updated: 09/06/2020 (dd-mm-yyyy)

%Modeling PF for a Substation with a Wind Farm
%A case study is loaded with 13 strings of Wind Turbine Generators 
%Step 2 involves the modeling of the power flow in situations. These are
%useful to determine limiting constraints in the substation's power flow

%Q demand at PCC is dependent of the total generation for the specific wind
%speed
%Q demand is divided equally among all string by means of a normal
%distribution with means going from -Q demand to +Q demand


%If the behaviour around certain Q setpoints want to be tested, the
%Q_dem_tot and Q_mean must be changed. 

%3 tests are ran:
% v= 4m/s, this is low generation wind speed
% v= 7.5 m/s, this is nominal generation wind speed
% v= 14 m/s, this is high generation wind speed
%This value has to be updated for each test run

clear;

%% Step 1: Load data of the case
mpc_substation_system13 = loadcase('system_13');
P_idx = find(mpc_substation_system13.gen(:,2)~=0); %Find indexes of generating strings in generator data
P_wt_idx = P_idx(5:end);
P_pv_idx = P_idx(1:4);

%Indicate branch limits for components that cause the constraints
B_limit_1_2 = 400;
B_limit_4_7 = 185.19;
B_limit_4_23 = 185.19;
B_limit_6_12 = 185.19;
B_limit_6_17 = 185.19;

%disconnect PV to run only WTG test

for(i=1:length(P_pv_idx))
    
mpc_substation_system13.gen(P_pv_idx(i),8) = 0; %disconnect PV to run only WTG test

end 

%Initialize active power generation

mpc_substation_system13.bus(11,2) = 4;
mpc_substation_system13.bus(16,2) = 4;
mpc_substation_system13.bus(22,2) = 4;
mpc_substation_system13.bus(27,2) = 4;

mpc_substation_system13.branch(10,11) = 0;
mpc_substation_system13.branch(16,11) = 0;
mpc_substation_system13.branch(23,11) = 0;
mpc_substation_system13.branch(29,11) = 0;

%define all constants
%WTG Farm
  N_daily_dispatch = 300; %total available wind dispatch profiles
  N_strings = 13; %total Wind Turbine Generators strings
  v_c_in = 3 ; %general cut in speed in m/s
  v_r = 14 ; %general rated speed in m/s, fixed for step 2a
  v_c_off = 25 ; %general cut off speed in m/s
  v_w = zeros(N_daily_dispatch,1);
    
  %Initialize wind speed in m/s
  for (i=1:N_daily_dispatch) 
    v_w(i,:) = 7.5; 
  end  

  %Maximum P(MW) and Q(MVAr) for each string 1-13 
  P_wt_max = [33.6 29.4 29.4 16.8 32 29.4 16 28 33 29.15 28 28.8 28.8];
  Q_wt_max = [21.2 18.55 19.15 10.9 22.4 19.15 12.2 19.6 22.4 13.248 19.6 19.6 19.6] ; 
  

   %% Step 2: Generate random variables
  
%%Generate a random wind distribution with shape parameters provided for
%%case study
% 
%  a= 2.54; %shape parameter in p.u., the lower the more uniformly spread, 1.2 is best
%  b= 7.86;%scale parameter in m/s, the higher the more the probability is spread, 12 is best
%  rand('twister',5489); % Fix seeds of the random number generator
% % rng( 'default'); %specifies seed for the random number generator, seed = 0
%  v_w = wblrnd(b,a, [N_daily_dispatch, 1]); %Initialise matrix for random wind speeds
%  v_low =zeros(1,length(v_w));
%  v_high=zeros(1,length(v_w));
%  for(q=1:length(v_w))
%  if(v_w(q,:) <= 4)
%      v_low(1,q) = v_w(q,1);
%      amount_low = find(v_low(1,:)~=0)
%  elseif(v_w(q,:) >= 13)
%      v_high(1,q) = v_w(q,1);
%      amount_high = find(v_high(1,:)~=0)
%      
%  end
%  end
%  
  
%Calculate generated power at given wind speed

P_wt = zeros(N_daily_dispatch, N_strings);


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

%Initialize active power generation in system topology
P_wt_total = sum(P_wt, 2);
mpc_substation_system13.gen(6:18,2) = P_wt(1,:); 
mpc_substation_system13.bus(1,3) = P_wt_total(1,:); 
mpc_substation_system13.gen(2:5,2) = 0;

%% %% Random Qs for WTG strings

%Approximate non-ideal P-Q relation of the WTGs 
rc_string_in = [6.625  6.625   11.969  13.625  6.22    11.969  6.22    6.22    11.2    7.36    6.22    12.25 10.64]; %specifies slope MVAr/MW at beginning for each WTG string
P_reg_in = [0.1 0.1 0.0544  0.047619    0.1 0.0544  0.1 0.1 0.0606  0.061   0.1 0.0556  0.0556]; %Percentage of total power to reach Q max in capability curve
rc_string_end = [1.5    1.5     0.507   3.143   1.48   0.507   1.8519  1.48    0.4408  0.4047   1.48    0.531   0.531]; %specifies slope MVAr/MW at end for each WTG string
P_reg_end = [0.88   0.88    0.372   0.917   0.83125   0.372   0.83125   0.83125   0.3023 0.225  0.83125   0.346   0.346]; %Percentage of total power to reach final Q in capability curve
 
Q_wt = zeros (N_daily_dispatch, N_strings);


    for i=1:N_strings
        for j=1:N_daily_dispatch  
            if (P_wt(j,i) < (P_reg_in(i)*P_wt_max(i)))
                Q_wt(j,i) = rc_string_in(i)*P_wt(j,i);
                if (Q_wt(j,i) >= Q_wt_max(i))
                    Q_wt(j,i) = Q_wt_max(i);
                end 
            elseif ((P_wt(j,i) >= (P_reg_in(i)*P_wt_max(i))) && (P_wt(j,i) < (P_reg_end(i)*P_wt_max(i))))
                Q_wt(j,i) = Q_wt_max(i);
            elseif ((P_wt(j,i) >= (P_reg_end(i)*P_wt_max(i))))
                Q_wt(j,i) = Q_wt_max(i) - rc_string_end(i).*(P_wt(j,i)-P_reg_end(i)*P_wt_max(i));
                if (Q_wt(j,i) < 0)
                    Q_wt(j,i) = 0;
                end 
            end 
        end
    end


Q_dem_tot = sum(Q_wt(1,:)); %Calculate the total Q generated by the strings
Q_mean_wt = zeros(1, N_strings);%matrix with mu for each string 

%Initialize normal random distribution
a = 0.08;%Percentage of standard deviation, derived from papers to be 2% to 8% of total Q

for(i=1:N_daily_dispatch)
    Q_wt_max_matrix(i,:) = Q_wt_max;
end

for(i=1:N_daily_dispatch)
    Q_std_wt(i,:) = unifrnd(1,Q_wt_max_matrix(i,:));
end

%Generate different means for the random distributions

for(j=1:N_strings)
for(i=1:(N_daily_dispatch))
    if(i <= (N_daily_dispatch/20))
        Q_mean_wt1(i,j) = -20*(Q_dem_tot/ N_strings)/20 ;
    elseif ((i > (1*N_daily_dispatch/20)) && (i <=(2*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt2(i,j) = -18*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (2*N_daily_dispatch/20)) && (i <=(3*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt3(i,j) = -16*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (3*N_daily_dispatch/20)) && (i <=(4*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt4(i,j) = -14*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (4*N_daily_dispatch/20)) && (i <=(5*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt5(i,j) = -12*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (5*N_daily_dispatch/20)) && (i <=(6*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt6(i,j) = -10*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (6*N_daily_dispatch/20)) && (i <=(7*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt7(i,j) = -8*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (7*N_daily_dispatch/20)) && (i <=(8*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt8(i,j) = -6*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (8*N_daily_dispatch/20)) && (i <=(9*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt9(i,j) = -4*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (9*N_daily_dispatch/20)) && (i <=(10*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt10(i,j) = -2*(Q_dem_tot/ N_strings)/20;  
    end
    elseif ((i > (10*N_daily_dispatch/20)) && (i <=(11*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt11(i,j) = 0;  
    end
    elseif ((i > (11*N_daily_dispatch/20)) && (i <=(12*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt12(i,j) = 2*(Q_dem_tot/ N_strings)/20;  
    end
        elseif ((i > (12*N_daily_dispatch/20)) && (i <=(13*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt13(i,j) = 4*(Q_dem_tot/ N_strings)/20;  
    end
        elseif ((i > (13*N_daily_dispatch/20)) && (i <=(14*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt14(i,j) = 6*(Q_dem_tot/ N_strings)/20;  
    end
        elseif ((i > (14*N_daily_dispatch/20)) && (i <=(15*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt15(i,j) = 8*(Q_dem_tot/ N_strings)/20;  
    end
        elseif ((i > (15*N_daily_dispatch/20)) && (i <=(16*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt16(i,j) = 10*(Q_dem_tot/ N_strings)/20;  
    end
    
        elseif ((i > (16*N_daily_dispatch/20)) && (i <=(17*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt17(i,j) = 12*(Q_dem_tot/ N_strings)/20;  
    end
    
        elseif ((i > (17*N_daily_dispatch/20)) && (i <=(18*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt18(i,j) = 14*(Q_dem_tot/ N_strings)/20;  
    end
    
        elseif ((i > (18*N_daily_dispatch/20)) && (i <=(19*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt19(i,j) = 16*(Q_dem_tot/ N_strings)/20;  
    end
    
        elseif ((i > (19*N_daily_dispatch/20)) && (i <=(20*N_daily_dispatch/20)))
    for(i=1:(N_daily_dispatch/20))
        Q_mean_wt20(i,j) = 20*(Q_dem_tot/ N_strings)/20;  
    end
end
end
end

for(i=1:N_strings) %makes a 1x13
  for(j=1:(N_daily_dispatch/20)) %makes a 96x1
    samples_Q_wt(j,i)= normrnd(Q_mean_wt1(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q in 96x13, some values violating the constraints
  end
  for(j=(1*N_daily_dispatch/20):(2*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt2(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(2*N_daily_dispatch/20):(3*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt3(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(3*N_daily_dispatch/20):(4*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt4(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(4*N_daily_dispatch/20):(5*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt5(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(5*N_daily_dispatch/20):(6*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt6(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(6*N_daily_dispatch/20):(7*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt7(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(7*N_daily_dispatch/20):(8*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt8(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(8*N_daily_dispatch/20):(9*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt9(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(9*N_daily_dispatch/20):(10*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt10(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
  for(j=(10*N_daily_dispatch/20):(11*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt11(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
  for(j=(11*N_daily_dispatch/20):(12*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt12(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
  end
    for(j=(12*N_daily_dispatch/20):(13*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt13(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
    for(j=(13*N_daily_dispatch/20):(14*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt14(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
    for(j=(14*N_daily_dispatch/20):(15*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt15(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
    for(j=(15*N_daily_dispatch/20):(16*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt16(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
    for(j=(16*N_daily_dispatch/20):(17*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt17(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
    for(j=(17*N_daily_dispatch/20):(18*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt18(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
    for(j=(18*N_daily_dispatch/20):(19*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt19(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
    for(j=(19*N_daily_dispatch/20):(20*N_daily_dispatch/20))
       samples_Q_wt(j,i)= normrnd(Q_mean_wt20(i), a*Q_std_wt(i), [1, 1]); %Creates random samples for Q 
    end
end

%% Step 3: Update case file and run Power Flow

stop =0;

while(stop == 0)
   for(i=1:N_daily_dispatch) 
   for (j=1:N_strings)

          if(i < N_daily_dispatch)
        mpc_substation_system13.gen(P_wt_idx(j),3) = samples_Q_wt(i,j);
        
        %Initialize active power generation in system topology
        mpc_substation_system13.bus(1,4) = sum(samples_Q_wt(i,j));

  
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
              
        %Save relevant data for each power flow  
        Vmag(:,i) = results_pf.bus(:,8); %Saves voltage magnitude for every bus
        Vangle(:,i) = results_pf.bus(:,9); %Saves angle magnitude for every bus
        
        %Save power flows in MW from different branches a_b with a bus, b bus 2
        P_flow_1_2(:,i) = results_pf.branch(1,16); %Branch to PCC
        P_flow_4_23(:,i) = results_pf.branch(25,16); %Branch 23 to Input transformer 1
        P_flow_4_7(:,i) = results_pf.branch(6,16); %Branch 7 to Input transformer 1
        P_flow_6_12(:,i) = results_pf.branch(12,16); %Branch 12 to Input transformer 2
        P_flow_6_17(:,i) = results_pf.branch(18,16); %Branch 17 to Input transformer 2
        
        %Store loss results with the specific branch in which it is
        %generated
        results_loss_real(:,i) = real(results_loss);
        results_loss_im(:,i) = imag(results_loss);
        results_inject_t(:,i) = results_inject;
        results_loss_matrix = [bus_from bus_to results_loss_real results_loss_im    results_inject_t];
 
          end 

          elseif(i == N_daily_dispatch)
          
            
         if (results_pf.success == 1)
        mpc_substation_system13.gen(P_wt_idx(j),3) = samples_Q_wt(i,j);
        mpc_substation_system13.bus(1,4) = sum(samples_Q_wt(i,j));

        %Save relevant data for each power flow  
        Vmag(:,i) = results_pf.bus(:,8); %Saves voltage magnitude for every bus
        Vangle(:,i) = results_pf.bus(:,9); %Saves angle magnitude for every bus
        
        %Save power flows in MW from different branches a_b with a bus, b bus 2
        P_flow_1_2(:,i) = results_pf.branch(1,16); %Branch to PCC
        P_flow_4_23(:,i) = results_pf.branch(25,16); %Branch 23 to Input transformer 1
        P_flow_4_7(:,i) = results_pf.branch(6,16); %Branch 7 to Input transformer 1
        P_flow_6_12(:,i) = results_pf.branch(12,16); %Branch 12 to Input transformer 2
        P_flow_6_17(:,i) = results_pf.branch(18,16); %Branch 17 to Input transformer 2
        
        %Store loss results with the specific branch in which it is
        %generated
        results_loss_real(:,i) = real(results_loss);
        results_loss_im(:,i) = imag(results_loss);
        results_inject_t(:,i) = results_inject;
        results_loss_matrix = [bus_from bus_to results_loss_real results_loss_im    results_inject_t];
 
        %Store data about the matrix location and value of the dispatches
        %and busses violating system constraints

        constraints_v_low = Vmag(Vmag < 0.8567); %Min voltage in p.u.
        constraints_v_high = Vmag(Vmag > 1.1424);  %Max voltage in p.u.
     
        idx_constraints_v_low = zeros(28, N_daily_dispatch); %28 is amount of busses in the circuit
for (a=1:N_daily_dispatch)

    
        idx_constraints_v_low(:,a) = (Vmag(:,a) < 0.8567); 
        idx_constraints_v_high(:,a) = (Vmag(:,a) > 1.1424);
        [idx_r_Q_v_low idx_c_Q_v_low] = find(idx_constraints_v_low~=0); %find location in matrix where there is v low
        [idx_r_Q_v_high idx_c_Q_v_high] = find(idx_constraints_v_high~=0);
end
        %Calculate lowest voltage magnitude value on the dispatches where a
        %violation has occured
        unique_idx_c_Q_v_low = unique(idx_c_Q_v_low);
        unique_idx_r_Q_v_low = unique(idx_r_Q_v_low);
            V_lowest = Vmag(:,unique_idx_c_Q_v_low);
            
        for(b=1:length(unique_idx_c_Q_v_low))
             constraints_v_lowest(:,b) =  min(V_lowest(:,b));
        end    
        
        %Calculate highest voltage magnitude value on the dispatches where a
        %violation has occured
        
        unique_idx_c_Q_v_high = unique(idx_c_Q_v_high);
        unique_idx_r_Q_v_high = unique(idx_r_Q_v_high);
            V_highest = Vmag(:,unique_idx_c_Q_v_high);
            
          for(c=1:length(unique_idx_c_Q_v_high))
             constraints_v_highest(:,c) =  max(V_highest(:,c));
          end    
  
 %Find Q demands, P loss, Q loss, Q injection corresponding to the
 %dispatches in which voltage violations occur
 
        Q_dem_constraints_v_low = samples_Q_wt(idx_c_Q_v_low,:); 
        Q_dem_constraints_v_high = samples_Q_wt(idx_c_Q_v_high,:); 
        
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
        
        %Print losses data per branch
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
        ylabel('Q loss(MVAr)','FontSize',12);
        title('Reactive power loss in MVAr for branch','FontSize',12)
        hold on;
        
        subplot(3,1,3)
        plot(results_inject_t,'o');
        xlabel('Branch #','FontSize',12);
        ylabel('Q injected(MVAr)','FontSize',12);
        title('Reactive power inject in MVAr for branch','FontSize',12)
        hold on;
        
        %Plot power flows through branches of importance
        figure(6)

        plot(P_flow_1_2,'o');
        xlabel('Dispatch','FontSize',12);
        ylabel('Power flow (MW)','FontSize',12);
        yline(B_limit_1_2 ,'r','LineWidth',1);
        title('Power flow through PCC Branch 1-2','FontSize',12);
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
        
        P_loss_tot_v_low = unique(P_loss_tot_v_low);
        Q_loss_tot_v_low = unique(Q_loss_tot_v_low);
        Q_inj_tot_v_low = unique(Q_inj_tot_v_low); 
        
        P_loss_tot_v_high = sum(P_loss_v_high);
        Q_loss_tot_v_high = sum(Q_loss_v_high);
        Q_inj_tot_v_high = sum(Q_inj_v_high);
        
                
        P_loss_tot_v_high = unique(P_loss_tot_v_high);
        Q_loss_tot_v_high = unique(Q_loss_tot_v_high);
        Q_inj_tot_v_high = unique(Q_inj_tot_v_high);

        %Calculate total Q_demand
    
        Q_tot = sum(samples_Q_wt,2);
        
        %Calculate total Q_demand for dispatches causing violation

        Q_tot_constraints_v_low = sum(Q_dem_constraints_v_low,2);
        Q_tot_constraints_v_high = sum(Q_dem_constraints_v_high,2);
         
        Q_tot_constraints_v_low = unique(Q_tot_constraints_v_low);
        Q_tot_constraints_v_high = unique(Q_tot_constraints_v_high);
        
        %add zeros in order to plot with equal amount of sample points
%         
%         P_loss_tot_v_low(numel(Q_tot_constraints_v_low)) = 0;
%         Q_loss_tot_v_low(numel(Q_tot_constraints_v_low)) = 0;
%         Q_inj_tot_v_low(numel(Q_tot_constraints_v_low)) = 0;
%         
%         P_loss_tot_v_high(numel(Q_tot_constraints_v_high)) = 0;
%         Q_loss_tot_v_high(numel(Q_tot_constraints_v_high)) = 0;
%         Q_inj_tot_v_high(numel(Q_tot_constraints_v_high)) = 0;
        
        %Plot 3D plot Q demand, P loss and (minimum and maximum) average voltage
        figure(3)
       
        plot3(Q_tot, P_loss_tot, v_avg,'o','Color', 'r', 'MarkerSize',5); %plot for avg
        hold on;
        %plot3(Q_tot_constraints_v_low, P_loss_tot_v_low,constraints_v_lowest,'o','Color', 'b', 'MarkerSize',5); %plot for v constraints low
        hold on;
        %plot3(Q_tot_constraints_v_high, P_loss_tot_v_high, constraints_v_highest,'o','Color', 'k', 'MarkerSize',5); %plot for v constraints low
        hold on;
        
        grid on;
        
        xlabel('Q demand (MVAr)','FontSize',12);
        ylabel('P loss total (MW)','FontSize',12);
        zlabel('V average magnitude (p.u.)','FontSize',12);
        
        legend('7.5 m/s');
        %Plot 3D plot Q loss for minimum and average voltage
        figure(4)
      
        plot3(Q_tot, Q_loss_tot, v_avg,'o','Color', 'r', 'MarkerSize',5); %plot for avg
        hold on;
        %plot3(Q_tot_constraints_v_low, Q_loss_tot_v_low, constraints_v_lowest,'o','Color', 'b', 'MarkerSize',5); %plot for v constraints low
        hold on;
        %plot3(Q_tot_constraints_v_high, Q_loss_tot_v_high', constraints_v_highest,'o','Color', 'k', 'MarkerSize',5); %plot for v constraints low
        hold on;
        
        grid on;
        
        xlabel('Q demand (MVAr)','FontSize',12);
        ylabel('Q loss total (MVAr)','FontSize',12);
        zlabel('V average magnitude (p.u.)','FontSize',12);
       
        legend('7.5 m/s');
        %Plot 3D plot for minimum and average
        figure(5)
        
        plot3(Q_tot, Q_inj_tot,v_avg,'o','Color', 'r', 'MarkerSize',5); %plot for avg
        hold on;
        %plot3(Q_tot_constraints_v_low,Q_inj_tot_v_low, constraints_v_lowest,'o','Color', 'b', 'MarkerSize',5); %plot for v constraints low
        hold on;
        %plot3(Q_tot_constraints_v_high, Q_inj_tot_v_high', constraints_v_highest,'o','Color', 'b', 'MarkerSize',5); %plot for v constraints low
        hold on;
        
        grid on;
        
        xlabel('Q demand (MVAr)','FontSize',12);
        ylabel('Q INJ total (MVAr)','FontSize',12);
        zlabel('V average magnitude (p.u.)','FontSize',12);
        
        legend('7.5 m/s');
        