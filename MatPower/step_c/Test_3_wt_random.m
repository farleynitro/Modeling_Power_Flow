close all;
clear all;
clc;

  v_c_in = 3 ; %general cut in speed in m/s
  v_r = 14 ; %general rated speed in m/s, fixed for step 2a
  v_c_off = 23 ; %general cut off speed in m/s
  
N_daily_dispatch = 96;
N_strings = 13;
P_wt_max = [33.6 29.4 29.4 28.8 32 29.4 16 28 33 29.5 28 28.8 28.8]  %Rated power in MW for each string 1-13
Q_wt_max = [21.2 18.55 19.15 10.9 22.4 19.15 12.2 19.6 22.4 13.248 19.6 19.6 19.6] ; %Max positive/negative reactive power in MVAr generated for eacht string 1-13

%%Generate a random wind distribution

a= 1.2; %shape parameter in p.u., the lower the more uniformly spread
b= 12;%scale parameter in m/s, the higher the more the probability is spread
%rand('twister',5489); % Fix seeds of the random number generator
rng( 'default'); %specifies seed for the random number generator, seed = 0
v_w = wblrnd(b,a, [N_daily_dispatch, 1]); %Initialise matrix for random wind speeds

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

%%Random Qs for WTG module


%Q_wt_max = [21.2   18.55   19.15   10.9    22.4    19.15   12.2    19.6    22.4    13.248  19.6    19.6    19.6] ; %Max positive/negative reactive power generated for eacht string 1-13
rc_string_in = [6.625  6.625   11.969  13.625  6.22    11.969  6.22    6.22    11.2    7.36    6.22    12.25 10.64]; %specifies slope MVAr/MW at beginning for each WTG string
P_reg_in = [0.1 0.1 0.0544  0.047619    0.1 0.0544  0.1 0.1 0.0606  0.061   0.1 0.0556  0.0556]; %Percentage of total power to reach Q max in capability curve
rc_string_end = -[1.5    1.5     0.507   3.143   1.48   0.507   1.8519  1.48    0.4408  0.4047   1.48    0.531   0.531]; %specifies slope MVAr/MW at end for each WTG string
%new_Q_max = [0.72   0.72    0.5117   0.596   0.65    0.5117   0.59    0.643   0.547   0.302   0.643 0.4898   0.4898]' * Q_wt_max; %Calculates new available MVAr at P_wt_max
P_reg_end = [0.88   0.88    0.372   0.917   0.83125   0.372   0.83125   0.83125   0.3023 0.225  0.83125   0.346   0.346]; %Percentage of total power to reach final Q in capability curve

samples_Q = zeros(N_daily_dispatch,N_strings); %matrix for the random Q_demands
Q_wt = zeros (N_daily_dispatch, N_strings);


for(j=1:N_daily_dispatch)
    for(i=1:N_strings)
        %check on boundary values 
%           k = 1;
        if (P_wt(j,i) < (P_reg_in(i) .*P_wt_max(i)))
              Q_wt(j,i) = rc_string_in(i).*P_wt(j,i);
                if (Q_wt(j,i) >= Q_wt_max(i))
                    Q_wt(j,i) = Q_wt_max(i);
                end 
        elseif ((P_wt(j,i) >= (P_reg_in(i).*P_wt_max(i))) && (P_wt(j,i) < (P_reg_end(i).*P_wt_max(i))))
               
               Q_wt(j,i) = Q_wt_max(i);
%                if (k<=N_strings)
%                k = k + 1;
%             end
        elseif ((P_wt(j,i) >= (P_reg_end(i).*P_wt_max(i))) && (P_wt(j,i) <= (P_wt_max(i)))) %%this works
              Q_wt(j,i) = Q_wt_max(i) + rc_string_end(i).*(P_wt(j,i)-P_reg_end(i).*P_wt_max(i));
              
        end 
    end
end

Q_mean = zeros(1, N_strings); %matrix with mu for each string 1-13

%samples_Q = zeros(length(Q_mean),1); %initialize matrix for Q random samples

for i=1:length(Q_mean)
    Q_std(i)=unifrnd(1,Q_wt(i)); %make a uniform distribution for the standard deviation string 1-13
end


for(i=1:N_strings) %makes a 1x13
  for(j=1: N_daily_dispatch) %makes a 96x1
    samples_Q(j,i)= normrnd(Q_mean(i), 0.75*Q_std(i),[1,1]); %Creates random samples for Q in 96x13, some values violating the constraints
  end
end


constraints_neg = zeros(N_daily_dispatch,N_strings);
constraints_pos = zeros(N_daily_dispatch,N_strings);
for(i=1:N_daily_dispatch)   
for(k=1: N_strings)
    idx_constraints_neg(:,k) = (samples_Q(:,k) < -Q_wt(:,k)); %Note values greater than max reactive 
    idx_constraints_pos(:,k) = (samples_Q(:,k) >  Q_wt(:,k));
    
    constraints_neg(i,k) = idx_constraints_neg(i,k).*-Q_wt(i,k);
    constraints_pos(i,k) = idx_constraints_pos(i,k).*-Q_wt(i,k);


    %constraints_neg(:,k)= samples_Q(samples_Q < -Q_wt(:,k)); %Note values greater than positive max reactive 
    %constraints_pos(:,k) = samples_Q(samples_Q > Q_wt(:,k));  %Note values greater than negative max reactive
% %     
   end
   end
 %constraints_neg(:,k) = idx_constraints_neg(:,k).*-Q_wt(:,k);
 %constraints_pos(:,k)=  idx_constraints_pos(:,k).*Q_wt(:,k);

%constraints_neg(:,k)= samples_Q(samples_Q < -Q_wt(:,k)); %Note values greater than positive max reactive
%constraints_pos(:,k) = samples_Q(samples_Q > Q_wt(:,k));  %Note values greater than negative max reactive


