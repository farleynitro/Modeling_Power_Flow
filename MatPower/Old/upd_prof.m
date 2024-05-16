%can be improved by varargin
function [Pd,Q_d_pos,Q_d_neg, n_bus] = upd_prof_wind(a,b, offset,wpp.dispatch,mpc.bus)

for(i=0: a) % a being the specified time period for the wind profile
     for ( j= 0:b) % 13 being the amount of strings
         n_bus = offset + j; %calculates actual bus assigned to the string
         v_w = wpp.dispatch(i,j); %gets wind speed out of dispatch received
         
         %calculate active and reactive power
         
         Pd = wpp_p_d(v_w, P_r_A,P_r_B,P_r_C,P_r_D)
         [Q_d_pos,Q_d_neg]= Q_wpp(Pd, string_n) 
         %varargin to replace bus for mpc.generator or mpc.branch?
         mpc.bus(n_bus, Pd) = Pd;
         mpc.bus(n_bus, Qmin) = Q_d_neg;
         mpc.bus(n_bus, Qmax) = Q_d_pos;
         
         j++; % j+1 until every string is updated;
         
         i++; %i+1 to go to next time period
     end
end

 %insert same for solar irradiance and temperature