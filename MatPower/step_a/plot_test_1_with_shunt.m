

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
        samples_Q_constraints = samples_Q(4,:)
        Q_total_constraints = sum(samples_Q_constraints);
         
        %Plot 3D plot P loss for minimum and average voltage
        figure(3)
       
        plot3(p_loss_tot',Q_total, v_avg','o','Color', 'r', 'MarkerSize',5);
        hold on;
        plot3(p_loss_tot',Q_total_constraints, v_avg_min,'o','Color', 'b', 'MarkerSize',5);
        grid on;
        ylabel('Q demand (MVAr)');
        xlabel('P loss total (MW)');
        zlabel('V average magnitude (p.u.)');
        
        %Plot 3D plot Q loss for minimum and average voltage
        figure(4)
        
        plot3(q_loss_tot',Q_total, v_avg','o','Color', 'r', 'MarkerSize',5);
        hold on;
        plot3(q_loss_tot',Q_total_constraints, v_avg_min','o','Color', 'b', 'MarkerSize',5);
        grid on;
        ylabel('Q demand (MVAr)');
        xlabel('Q loss total (MVAr)');
        zlabel('V average magnitude (p.u.)');
       
        %Plot 3D plot for minimum and average
        figure(5)
        
        plot3(q_inj_tot',Q_total, v_avg','o','Color', 'r', 'MarkerSize',5);
        hold on;
        plot3(q_inj_tot',Q_total_constraints, v_avg_min','o','Color', 'b', 'MarkerSize',5);
        grid on;
        ylabel('Q demand (MVAr)');
        xlabel('Q INJ total (MVAr)');
        zlabel('V average magnitude (p.u.)');
        
       
    