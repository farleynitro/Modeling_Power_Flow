function [new_angle,new_mag,new_P,new_Q,done] = nr_method(mag,angle,P,Q)




P1  =(mag*7.1621*cosd(angle-118.6)+7.1621*cosd(-80)*(mag^2));
Q1 = (mag*7.1621*sind(angle-118.6)+7.1621*sind(-80)*(mag^2));

%Calculate derivatives
diff_P1_angle = -7.1621*mag*sind(angle-118.6); 
diff_P1_mag = 7.1621*cosd(angle-118.6)+7.1621*cosd(-80)*2*mag;

diff_Q1_angle = 7.1621*mag*cosd(angle-118.6);
diff_Q1_mag = 7.1621*sind(angle-118.6)-7.1621*sind(-80)*2*mag;

%equation to be solved
delta_flow = [(P-P1);(Q-Q1)];

%Put derivatives in Jacobian
jacob = [diff_P1_angle diff_P1_mag; diff_Q1_angle diff_P1_mag];

%Solve linear equation
x = linsolve(jacob,delta_flow);

%specify new angle and magnitude
new_angle = angle - x(1,:);
new_mag = angle - x(2,:);
new_P = P-P1
new_Q = Q-Q1


%determine if result is below threshold

if((P-P1) < 0.1)
    done = 1;
    
else 
    done = 0;
    
end
