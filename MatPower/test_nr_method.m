%%Authors: F. Rimon and M. Mastouri
%Last edited: 05/06/2020

%This is a method that simulates the hand calculations done to solve a
%Newton Raphson Power flow

%The case study used is a branch of the case study of the Bachelor Thesis
%Project of Group 09.01
close all;
clear all;
clc
%% Power flow for Node 1(Bus 24) of Branch 23-24 with Q from samples_Q data #Dispatch 1000
% %Workspace needs to be included

done = 0;
iter = 1; %count iterations needed
while(done == 0 )
    
mag = 1;
angle = 0;
P_1 = 0.294;
Q_1 = 0.0169;

P1  =(mag*7.1621*cosd(angle-118.6)+7.326*cosd(-80)*(mag^2));
Q1 = (mag*7.1621*sind(angle-118.6)-7.326*sind(-80)*(mag^2));

Calculate derivatives
diff_P1_angle = -7.1621*mag*sind(angle-118.6); 
diff_P1_mag = 7.1621*cosd(angle-118.6)+7.326*cosd(-80)*2*mag;

diff_Q1_angle = 7.1621*mag*cosd(angle-118.6);
diff_Q1_mag = 7.1621*sind(angle-118.6)-7.326*sind(-80)*2*mag;

equation to be solved
new_P = P_1-P1;
new_Q = Q_1-Q1;
delta_flow = [new_P;new_Q];

Put derivatives in Jacobian
jacob = [diff_P1_angle diff_P1_mag; diff_Q1_angle diff_Q1_mag];

Solve linear equation
x = linsolve(jacob,delta_flow);
delta_angle_degree = (x(1,:)./(2*pi))*360; %convert radian angle to degree angle

specify new angle and magnitude
new_angle = angle - delta_angle_degree;
new_mag = mag - x(2,:);

determine if result is below threshold

if(new_P <= 0.1)
    done = 1;
    iter = iter;
  
else 
    done = 0;
    mag = new_mag;
    angle = new_angle;
    iter = iter + 1;
end
end 

%After the final answers are determined, insert them manually into these
%equations to determine the losses in the system.
new_angle = 19.4;
new_mag = 1.0;
Z_l = 0.02356859504 + 0.134446281*i;
new_mag_phasor = new_mag*(cosd(new_angle)+sind(new_angle)*i);
old_mag_phasor = 0.9776*(cosd(18.5814)+sind(18.5814)*i);
I_23_24 = (new_mag_phasor - old_mag_phasor)/Z_l;
P_loss = (abs(I_23_24))^2*real(Z_l);
Q_loss = (abs(I_23_24))^2*imag(Z_l);

Base = 100;

P_loss_real = P_loss *Base;
Q_loss_real = Q_loss *Base;
%% Power flow for Node 2 (Bus 23) of Branch 23-24 with Q from samples_Q data #Dispatch 1000
%Workspace needs to be included
% 
% done_2 = 0;
% iter_2 = 1; %count iterations needed
% while(done_2 == 0 )
%     
% mag = 1;
% angle = 18;
% P_2 = -0.294;
% Q_2 = -0.0169;
% % solution = nr_method(0.9776,18.5814,0.294,0.000169)
% 
% P2  =(mag*7.18*cosd(angle-119.25)+7.326*cosd(-80)*(mag^2));
% Q2 = (mag*7.18*sind(angle-119.25)-7.326*sind(-80)*(mag^2));
% 
% %Calculate derivatives
% diff_P2_angle = -7.18*mag*sind(angle-119.25); 
% diff_P2_mag = 7.18*cosd(angle-119.25)+7.326*cosd(-80)*2*mag;
% 
% diff_Q2_angle = 7.18*mag*cosd(angle-119.25);
% diff_Q2_mag = 7.18*sind(angle-119.25)-7.326*sind(-80)*2*mag;
% 
% %equation to be solved
% new_P_2 = P_2-P2;
% new_Q_2 = Q_2-Q2;
% delta_flow = [new_P_2;new_Q_2];
% 
% %Put derivatives in Jacobian
% jacob = [diff_P2_angle diff_P2_mag; diff_Q2_angle diff_Q2_mag];
% 
% %Solve linear equation
% x = linsolve(jacob,delta_flow);
% delta_angle_degree = (x(1,:)./(2*pi))*360; %convert radian angle to degree angle
% 
% %specify new angle and magnitude
% new_angle = angle - delta_angle_degree;
% new_mag = mag - x(2,:);
% 
% %determine if result is below threshold
% 
% if(new_P_2 <= 0.00001)
%     done_2 = 1;
%     iter_2 = iter_2;
%   
% else 
%     done_2 = 0;
%     mag = new_mag;
%     angle = new_angle;
%     iter_2 = iter_2 + 1;
% end
% end 
%% Power flow for Node 2 of Branch 23-24 with Q from generator pf data #Dispatch 1000
% % %Workspace needs to be included
% 
% % done = 0;
% % iter = 1; %count iterations needed
% % while(done == 0 )
% %     
% % mag = 1.01;
% % angle = 20.7;
% % P_1 = 0.294;
% % Q_1 = 0.000256;
% % 
% % P1  =(mag*7.1621*cosd(angle-118.6)+7.326*cosd(-80)*(mag^2));
% % Q1 = (mag*7.1621*sind(angle-118.6)-7.326*sind(-80)*(mag^2));
% % 
% % Calculate derivatives
% % diff_P1_angle = -7.1621*mag*sind(angle-118.6); 
% % diff_P1_mag = 7.1621*cosd(angle-118.6)+7.326*cosd(-80)*2*mag;
% % 
% % diff_Q1_angle = 7.1621*mag*cosd(angle-118.6);
% % diff_Q1_mag = 7.1621*sind(angle-118.6)-7.326*sind(-80)*2*mag;
% % 
% % equation to be solved
% % new_P = P_1-P1;
% % new_Q = Q_1-Q1;
% % delta_flow = [new_P;new_Q];
% % 
% % Put derivatives in Jacobian
% % jacob = [diff_P1_angle diff_P1_mag; diff_Q1_angle diff_Q1_mag];
% % 
% % Solve linear equation
% % x = linsolve(jacob,delta_flow);
% % delta_angle_degree = (x(1,:)./(2*pi))*360; %convert radian angle to degree angle
% % 
% % specify new angle and magnitude
% % new_angle = angle - delta_angle_degree;
% % new_mag = mag - x(2,:);
% % 
% % determine if result is below threshold
% % 
% % if(new_P <= 0.001)
% %     done = 1;
% %     iter = iter;
% %   
% % else 
% %     done = 0;
% %     mag = new_mag;
% %     angle = new_angle;
% %     iter = iter + 1;
% % end
% % end 

% new_angle = 19.7;
% new_mag = 0.9801;
% Z_l = 0.02356859504 + 0.134446281*i;
% new_mag_phasor = new_mag*(cosd(new_angle)+sind(new_angle)*i);
% old_mag_phasor = 0.9776*(cosd(18.5814)+sind(18.5814)*i);
% I_23_24 = (new_mag_phasor - old_mag_phasor)/Z_l;
% P_loss = (abs(I_23_24))^2*real(Z_l);
% Q_loss = (abs(I_23_24))^2*imag(Z_l);
% 
% Base = 100;
% 
% P_loss_real = P_loss *Base;
% Q_loss_real = Q_loss *Base;