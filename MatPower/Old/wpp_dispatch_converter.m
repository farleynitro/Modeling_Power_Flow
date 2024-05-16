
%wind speed conversion
function [Pd] = wpp_p_d(v_w string)

for( string(i,:)= 0:(n_strings_t-1))
    j=0;
    if(v_w <= v_c_in) %% apply boundary conditions for which equation is valid
    Pd =0;

    elseif (v_w > v_c_in && v_w < v_c_off) 
    
    Pd = string(i,j).*(P_r_A * ((v_w)^3 -(v_r)^3)/((v_r)^3-(v_c)^3))+ string(i,j+1).*(P_r_B * ((v_w)^3 -(v_r)^3)/((v_r)^3-(v_c)^3)) + string(i,j+2).*(P_r_C * ((v_w)^3 -(v_r)^3)/((v_r)^3-(v_c)^3))+ string(i,j+3).*(P_r_D * ((v_w)^3 -(v_r)^3)/((v_r)^3-(v_c)^3));

    elseif (v_w >= v_c_off)
    Pd = string(i,j).*P_r_A + string(i,j+1).*P_r_B + string(i,j+2).*P_r_C + string(i,j+3).*P_r_D;
    
    i++;
    
    


%Calculate reactive power generated at each string
function [Q_d_pos,Q_d_neg]= Q_wpp(Pd,string_n)

rc_pos = 0;
rc_neg = 0;
%string 1
if ( Pd(i) <= 3.2 && Pd(i) >= 0 && string_n = 1)
   
    rc_pos = 6.625;
    rc_neg= -rc_pos;
    Q_d_= rc_pos*Pd(i);
    Q_d_neg= rc_neg*Pd(i);
    
    
elseif ( Pd(i) <= 3.2 && Pd(i) >= 29.2 && string_n = 1)
    
    Q_d_pos= 21.2;
    Q_d_neg= -21.2;
    
elseif ( Pd(i) <= 33.6 && Pd(i) >= 29.2 && string_n = 1)
 
    rc_pos = (21.2-15.2)/(33.6-29.2);
    rc_neg = -(-14--21.2)/(33.6-29.2);
    Q_d_pos= 21.2 - rc_pos*Pd(i);
    Q_d_neg= -21.2 - rc_neg*Pd(i);
end
%string 2
if ( Pd(i+1) <= 2.8 && Pd(i+1) >= 0 && string_n = 2)
    rc_pos = 6.625;
    Q_d_pos= rc_pos*Pd(i+1);
    Q_d_neg= -rc_pos*Pd(i+1);
  
    
elseif ( Pd(i+1) <= 25.55 && Pd(i+1) >= 2.8 && string_n = 2)
    Q_d_pos= 18.55;
    Q_d_neg= -18.55;
    
elseif ( Pd(i+1) <= 25.55 && Pd(i+1) >= 29.4 && string_n = 2)
    
    rc_pos = (18.55-13.3)/(29.4-25.55);
    rc_neg = (-12.25--18.55)/(29.4-25.55);
    Q_d_pos= 18.55 - rc_pos*Pd(i+1);
    Q_d_neg= -18.55 + rc_neg*Pd(i+1);
end
%string 3 
if ( Pd(i+2) <= 1.2 && Pd(i+2) >= 0 && string_n = 3)
    
    rc_pos = 16.35/1.2;
    rc_neg = -16.35/1.2;
    Q_d_neg = rc_neg*Pd(i);
    Q_d_pos = rc_pos*Pd(i);

elseif ( Pd(i+2) <= 1.6 && Pd(i+2) >= 1.2 && string_n = 3)
    
    rc_pos = (19.15-16.35)/(1.6-1.2);
    rc_neg = (-19.15--16.35)/(1.6-1.2);
    Q_d_neg = rc_neg*Pd(i+2);
    Q_d_pos = rc_pos*Pd(i+2);
    
elseif ( Pd(i+2) <= 10.95 && Pd(i+2) >= 1.6 && string_n = 3)

    Q_d_pos = 19.15;
    Q_d_neg = -19.15;
    
elseif ( Pd(i+2) <= 12.6 && Pd(i+2) >= 10.95 && string_n = 3)

    Q_d_pos = 16.9;
    Q_d_neg = -16.6;
elseif ( Pd(i+2) <= 28.4  && Pd(i+2) >= 12.6 && string_n = 3)
    
    Q_d_pos = 11.2;
    Q_d_neg = -11.2;
    
elseif ( Pd(i+2) <= 29.4 && Pd(i+2) >= 28.4 && string_n = 3)
 
    rc_pos = (9.8-11.2)/(29.4-28.4);
    rc_neg = (-8.8--11.2)/(29.4-28.4);
    Q_d_neg = rc_neg*Pd(i+2);
    Q_d_pos = rc_pos*Pd(i+2);

end     
    
%string 4
    
if ( Pd(i+3) <= 0.8 && Pd(i+3) >= 0 && string_n = 4)
    
    rc_pos = 6.2225;
    rc_neg = -rc_pos;
    Q_d_pos= rc_pos*Pd(i+3);
    Q_d_neg= rc_neg*Pd(i+3);
    
elseif ( Pd(i+3) <= 15.4 && Pd(i+3) >= 0.8 && string_n = 4)
    
    Q_d_pos= 22.4;
    Q_d_neg= -22.4;
    
elseif ( Pd(i+3)<= 16.8 && Pd(i+3) >= 15.4 && string_n = 4)
    
    rc_pos = (10.9-6.5)/(16.8-15.4);
    rc_neg = (-5.1--10.9)/(16.8-15.4);
    Q_d_pos= 10.9 - rc_pos*Pd(i+3);
    Q_d_neg= -10.9 + rc_neg*Pd(i+3);
end
    
if ( Pd(i+4) <= 3.6 && Pd(i+4) >= 0 && string_n = 5)
    
    rc_pos = 6.22;
    rc_neg = -rc_neg;
    Q_d_pos= rc_pos*Pd(i+4);
    Q_d_neg= rc_neg*Pd(i+4);
  
    
elseif ( Pd(i+4) <= 26.6 && Pd(i+4) >= 3.6 && string_n = 5)
    
    Q_d_pos= 22.4;
    Q_d_neg= -22.4;
    
elseif ( Pd(i+4) <= 32 && Pd(i+4) >= 26.6 && string_n = 5)
    
    rc_pos = (22.4-14.4)/(32-26.6);
    rc_neg = (-12--22.4)/(32-26.6);
    Q_d_pos= 22.4 - rc_pos*Pd(i+4);
    Q_d_neg= -22.4 + rc_neg*Pd(i+4);
end
%string 6  
if ( Pd(i+5) <= 1.2 && Pd(i+5) >= 0 && string_n = 6)
    
    rc_pos = 16.35/1.2;
    rc_neg = -16.35/1.2;
    Q_d_neg = rc_neg*Pd(i+5);
    Q_d_pos = rc_pos*Pd(i+5);

elseif ( Pd(i+5) <= 1.6 && Pd(i+5) >= 1.2 && string_n = 6)
    
    rc_pos = (19.15-16.35)/(1.6-1.2);
    rc_neg = (-19.15--16.35)/(1.6-1.2);
    Q_d_neg = rc_neg*Pd(i+5);
    Q_d_pos = rc_pos*Pd(i+5);
    
elseif ( Pd(i+5) <= 10.95 && Pd(i+5) >= 1.6 && string_n = 6)

    Q_d_pos = 19.15;
    Q_d_neg = -19.15;
    
elseif ( Pd(i+5) <= 12.6 && Pd(i+5) >= 10.95 && string_n = 6)

    Q_d_pos = 16.9;
    Q_d_neg = -16.6;
    
elseif ( Pd(i+5) <= 28.4 && Pd(i+5) >= 12.6 && string_n = 6)
    
    Q_d_pos = 11.2;
    Q_d_neg = -11.2;
    
elseif ( Pd(i+5) <= 29.4 && Pd(i+5) >= 28.4 && string_n = 6)
 
    rc_pos = (9.8-11.2)/(29.4-28.4);
    rc_neg = (-8.8--11.2)/(29.4-28.4);
    Q_d_neg = rc_neg*Pd(i+5);
    Q_d_pos = rc_pos*Pd(i+5);

end     
    
    
if ( Pd(i+6) <= 1.8 && Pd(i+6) >= 0 && string_n = 7)
       
    rc_pos = 6.22;
    rc_neg = -rc_pos;
    Q_d_pos= rc_pos*Pd(i);
    Q_d_neg= rc_neg*Pd(i);
  
elseif ( Pd(i+6) <= 13.3 && Pd(i+6) >= 1.8 && string_n = 7)
    
    Q_d_pos= 12.2;
    Q_d_neg= -12.2;
    
elseif ( Pd(i+6) <= 16 && Pd(i+6) >= 13.3 && string_n = 7)
    
    rc_pos = (12.2-7.2)/(16-13.3);
    rc_neg = (-6--12.2)/(16-13.3);
    Q_d_pos= 12.2 - rc_pos*Pd(i);
    Q_d_neg= -12.2 + rc_neg*Pd(i);
end
    
if ( Pd(i+7) <= P1 && Pd(i+7) >= P2 && string_n = 8)

elseif ( Pd(i+7) <= P1 && Pd(i+7) >= P2 && string_n = 8)
    
elseif ( Pd(i+7) <= P1 && Pd(i+7) >= P2 && string_n = 8)
    
    
if ( Pd(i+8) <= 1.35 && Pd(i+8) >= 0 && string_n = 9)
    
    rc_pos = (17.847-0)/(1.35);
    rc_neg = -(17.847-0)/(1.35);
    Q_d_pos= rc_pos*Pd(i);
    Q_d_neg= rc_neg*Pd(i);

elseif ( Pd(i+8) <= 2 && Pd(i+8) >= 1.35 && string_n = 9)
    
    rc_pos = (22.4-17.847)/(2-1.35);
    rc_neg = (-22.4--17.847)/(2-1.35);
    Q_d_pos= 17.847 - rc_pos*Pd(i);
    Q_d_neg= -17.847 + rc_neg*Pd(i);
    
elseif ( Pd(i+8) <= 9.975 && Pd(i+8) >= 2 && string_n = 9)
    
    Q_d_pos = 22.4;
    Q_d_neg = -22.4;

elseif ( Pd(i+8) <= 12 && Pd(i+8) >= 9.975 && string_n = 9)
    
    Q_d_pos = 19.4;
    Q_d_neg = -19.4;

elseif ( Pd(i+8) <= 31.25 && Pd(i+8) >= 12 && string_n = 9)
    
    Q_d_pos = 14;
    Q_d_neg = -14.3; 
    
elseif ( Pd(i+8) <= 31.25 && Pd(i+8) >= 33 && string_n = 9)
    
    rc_pos = (14-12.25)/(33-31.25);
    rc_neg = (-14.3--11.3)/(33-31.25);
    Q_d_pos= 14 - rc_pos*Pd(i);
    Q_d_neg= -14.3 + rc_neg*Pd(i);
    
end 
    
if ( Pd(i+9) <= P1 && Pd(i+9) >= P2 && string_n = 10)

elseif ( Pd(i+9) <= P1 && Pd(i+9) >= P2 && string_n = 10)
    
elseif ( Pd(i+9) <= P1 && Pd(i+9) >= P2 && string_n = 10)

    
elseif ( Pd(i+9) <= P1 && Pd(i+9) >= P2 && string_n = 10)
    
elseif ( Pd(i+9) <= P1 && Pd(i+9) >= P2 && string_n = 10)
    

elseif ( Pd(i+9) <= P1 && Pd(i+9) >= P2 && string_n = 10)
    
    
    
if ( Pd(i+10) <= P1 && Pd(i+10) >= P2 && string_n = 11)

elseif ( Pd(i+10) <= P1 && Pd(i+10) >= P2 && string_n = 11)
    
elseif ( Pd(i+10) <= P1 && Pd(i+10) >= P2 && string_n = 11)
    
    
    
if ( Pd(i+11) <= 1.35 && Pd(i+11) >= 0 && string_n = 12)
    
    rc_pos = (17.847-0)/(1.35);
    rc_neg = -(17.847-0)/(1.35);
    Q_d_pos= rc_pos*Pd(i);
    Q_d_neg= rc_neg*Pd(i);

elseif ( Pd(i+11) <= 1.6 && Pd(i+11) >= 1.35 && string_n = 12)
    
    rc_pos = (19.6-17.847)/(1.6-1.35);
    rc_neg = (-19.6--17.847)/(1.6-1.35);
    Q_d_pos= 17.847 + rc_pos*Pd(i);
    Q_d_neg= -17.847 + rc_neg*Pd(i);
    
elseif ( Pd(i+11) <= 9.975 && Pd(i+11) >= 1.6 && string_n = 12)
    
    Q_d_pos = 19.6;
    Q_d_neg = -19.6;   
    
elseif ( Pd(i+11) <= 12 && Pd(i+11) >= 9.975 && string_n = 12)
    
    Q_d_pos = 16.4;
    Q_d_neg = -15.7;  
    
elseif ( Pd(i+11) <= 27.4 && Pd(i+11) >= 12 && string_n = 12)
    
    Q_d_pos = 11;
    Q_d_neg = -10.6;  
    
elseif ( Pd(i+11) <= 28.8 && Pd(i+11) >= 27.4 && string_n = 12)
    
    rc_pos = (11-9.6)/(28.8-27.4);
    rc_neg = (-8.2--10.6)/(28.8-27.4);
    Q_d_pos= 11 - rc_pos*Pd(i);
    Q_d_neg= -10.6 + rc_neg*Pd(i);
end 
    
if ( Pd(i+12) <= P1 && Pd(i+12) >= P2 && string_n = 13)
    
    rc_pos = (17.847-0)/(1.35);
    rc_neg = -(17.847-0)/(1.35);
    Q_d_pos= rc_pos*Pd(i);
    Q_d_neg= rc_neg*Pd(i);
    
elseif ( Pd(i+12) <= P1 && Pd(i+12) >= P2 && string_n = 13)
    
    rc_pos = (19.6-17.847)/(1.6-1.35);
    rc_neg = (-19.6--17.847)/(1.6-1.35);
    Q_d_pos= 17.847 + rc_pos*Pd(i);
    Q_d_neg= -17.847 + rc_neg*Pd(i);
    
elseif ( Pd(i+12) <= P1 && Pd(i+12) >= P2 && string_n = 13)
    
    Q_d_pos = 19.6;
    Q_d_neg = -19.6; 
    
elseif ( Pd(i+12) <= P1 && Pd(i+12) >= P2 && string_n = 13)
    
        
    Q_d_pos = 16.4;
    Q_d_neg = -15.7;  

elseif ( Pd(i+12) <= P1 && Pd(i+12) >= P2 && string_n = 13)
    
    Q_d_pos = 11;
    Q_d_neg = -10.6;  
     
elseif ( Pd(i+12) <= P1 && Pd(i+12) >= P2 && string_n = 13)
    
    rc_pos = (11-9.6)/(28.8-27.4);
    rc_neg = (-8.2--10.6)/(28.8-27.4);
    Q_d_pos= 11 - rc_pos*Pd(i+12);
    Q_d_neg= -10.6 + rc_neg*Pd(i+12);

end