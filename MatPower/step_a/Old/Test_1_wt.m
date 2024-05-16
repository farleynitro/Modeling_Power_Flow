% P_gen_max = [33.6 29.4 29.4 28.8 32 29.4 16 28 33 29.5 28 28.8 28.8]  %Rated power for each string
% 
clear all;
close all;

 N_daily_dispatch = 96; %total available wind profiles
 N_strings = 13; %total strings
% 
% P_gen = zeros(N_daily_dispatch, N_strings); %Make matrix for generated nominal wind speed
% 
% for(i=1: length(P_gen_max))
%     for(j=1: N_daily_dispatch)
%         P_gen(j,i) = P_gen_max(i); %assigns 
%     end
% end

%Random Q_demands of Normal distribution for each string/generator



samples_Q = zeros(N_daily_dispatch,N_strings); %matrix for the random Q_demands
Q_gen_max = [21.2   18.55   19.15   10.9    22.4    19.15   12.2    19.6    22.4    13.248  19.6    19.6    19.6] ; %Max positive/negative reactive power generated for eacht string 1-13

Q_mean = zeros(1, N_strings); %mu for each string 1-13
% 
% samples_Q = zeros(length(Q_mean),1);

for i=1:length(Q_mean)
    Q_std(i)=unifrnd(1,Q_gen_max(i)); %make a uniform distribution for the standard deviation string 1-13
end


for(i=1:N_strings) %makes a 1x13
  for(j=1: N_daily_dispatch) %makes a 96x1
    samples_Q(j,i)= normrnd(Q_mean(i), 0.75*Q_std(i), [1, 1]); %Creates random samples for Q in 96x13
  end
end

%it only searches for the first value in every column, also searches per
%column and not per row
constraints_neg = zeros(N_daily_dispatch,N_strings);
constraints_pos= zeros(N_daily_dispatch,N_strings);
   for(k=1: N_strings)
%    constraints_neg = samples_Q(samples_Q < -Q_gen_max(k)); %Note values greater than positive max reactive 
%    constraints_pos = samples_Q(samples_Q > Q_gen_max(k));  %Note values greater than negative max reactive
% %     
    idx_constraints_neg(:,k) = (samples_Q(:,k) < -Q_gen_max(k)); %Note values greater than max reactive 
    idx_constraints_pos(:,k) = (samples_Q(:,k) >  Q_gen_max(k));
       
% 
%     [a,b] = (find(idx_constraints_pos)==1);
%     [c,d] = (find(idx_constraints_neg)==1);
   end 

% N_daily_dispatch = 96
% 

