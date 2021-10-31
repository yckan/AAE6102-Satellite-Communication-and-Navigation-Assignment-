% Input x0 - initialization of user position and clock bias [x0;y0;z0,t0]
% Output x_u - estimated user position and clock bias

function [x_u] = solve_x_user_LS(rho,x0,x_s,y_s,z_s,Delta_t_SV)
    i = 1;
    c = 299792458.0; % speed of light (m/s)
    while 1
        % Formulate H matrix, delta_rho, and delta_x for Least Square
        % Solution
        x_array(i,:) = x0; % Store solution by iteration
        rho_hat = sqrt((x_s-x0(1)).^2 + (y_s-x0(2)).^2 + (z_s-x0(3)).^(2)) + c* x0(4) - Delta_t_SV*c; % Calculate rho_hat (pseudorange at estimated position)
        delta_rho = rho_hat - rho;
        a_x = (x_s-x0(1))./sqrt((x_s-x0(1)).^2 + (y_s-x0(2)).^2 + (z_s-x0(3)).^(2)); % Construct H Matrix components
        a_y = (y_s-x0(2))./sqrt((x_s-x0(1)).^2 + (y_s-x0(2)).^2 + (z_s-x0(3)).^(2)); % Construct H Matrix components
        a_z = (z_s-x0(3))./sqrt((x_s-x0(1)).^2 + (y_s-x0(2)).^2 + (z_s-x0(3)).^(2)); % Construct H Matrix components
        H = [a_x a_y a_z ones(size(rho,1),1)]; 
        delta_x = inv(H'*H)*H'*delta_rho; % Solution for LS method
        delta_x(4) = -delta_x(4)/c; % Convert to clock bias by dividing speed of light
        dx_array(i,:) = delta_x; % store dx in array for figure plots
        norm_dx_array(i) = norm(delta_x(1:3)); % Store residual
        x0 = x0 + delta_x; % Update solution for next iteration
        i = i+1;
        if norm(delta_x(1:3)) <= 1e-4 || i > 10  % check delta position solution smaller than threshold
            x_u = x0;
            break;
        end
    end
    
    % Plot Figure for Assignments
    %plot3(x_array(1:size(x_array,1),1),x_array(1:size(x_array,1),2),x_array(1:size(x_array,1),3),'o-')
    %hold on
    %plot3(x_array(1,1),x_array(1,2),x_array(1,3),'rp','MarkerSize',20)
    %plot3(x_array(size(x_array,1),1),x_array(size(x_array,1),2),x_array(size(x_array,1),3),'rx','MarkerSize',20)
    %grid on
    %legend('intermediate estimation','initial position','final estimation')
    %xlabel('x_u (m)'); ylabel('y_u (m)'); zlabel('z_u (m)');
    %title('User Position Estimation by Iteration (m)')
    %figure(2)
    %plot(norm_dx_array,'bo-');
    %title('norm(\bf{\Deltax_u}) Variation by Iteration (m)')
    %xlabel('No. of Iteration')
    %ylabel('norm(\Deltax_u) (m)');
    %figure(5)
    %plot(1:size(dx_array,1),dx_array(:,1),'bo-')
    %title('\Deltax_u Variation by Iteration (m)')
    %xlabel('No. of Iteration')
    %ylabel('\Deltax_u (m)');
    %figure(6)
    %plot(1:size(dx_array,1),dx_array(:,2),'bo-')
    %title('\Deltay_u Variation by Iteration (m)')
    %xlabel('No. of Iteration')
    %ylabel('\Deltay_u (m)');
    %figure(7)
    %plot(1:size(dx_array,1),dx_array(:,3),'bo-')
    %title('\Deltaz_u Variation by Iteration (m)')
    %xlabel('No. of Iteration')
    %ylabel('\Deltaz_u (m)');
    %figure(8)
    %plot(1:size(dx_array,1),dx_array(:,4),'bo-')
    %title('\Deltat_u Variation by Iteration (sec)')
    %xlabel('No. of Iteration')
    %ylabel('\Deltat_u (sec)');
end
