

        %Calculate average voltage for specific Q_demand
        v_avg = mean(Vmag);
        %Calculate average minimum voltage for specific Q_demand
        v_avg_min = mean(constraints_v_low);
        
        %Calculate average losses for same Q_demand
        p_loss_tot = sum(results_loss_real);
        q_loss_tot = sum(results_loss_im);
        q_inj_tot = sum(results_inject_t);
        
        %Calculate maximum losses for same Q_demand
        p_loss_max = max(results_loss_real);
        q_loss_max = max(results_loss_im);
        q_inj_max = max(results_inject_t);
        
        %Calculate total Q_demand
        
        Q_total = sum(samples_Q,2);
        
        %Calculate total Q_demand for dispatches causing violation
        
        Q_total_constraints = 0;
                
        figure(3)
       
        plot3(p_loss_tot',Q_total, v_avg','o','Color', 'r', 'MarkerSize',5);
        
        ylabel('Q demand (MVAr)');
        xlabel('P loss total (MW)');
        zlabel('V average magnitude (p.u.)');
        
        grid on;
        hold on;
        
        figure(4)
        plot3(q_loss_tot',Q_total, v_avg','o','Color', 'r', 'MarkerSize',5);
        
        ylabel('Q demand (MVAr)');
        xlabel('Q loss total (MVAr)');
        zlabel('V average magnitude (p.u.)');
        
        grid on;
        hold on;
        figure(5)
        plot3(q_inj_tot',Q_total, v_avg','o','Color', 'r', 'MarkerSize',5);
        
        ylabel('Q demand (MVAr)');
        xlabel('Q INJ total (MVAr)');
        zlabel('V average magnitude (p.u.)');
        
        grid on;
        hold on;
  
   
        
        