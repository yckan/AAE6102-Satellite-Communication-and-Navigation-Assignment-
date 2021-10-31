function Delta_t_SV = solve_SV_clock_err(a,e,E_k,af0,af1,af2,t_s,toc)

    F = -4.442807633e-10; %Relativistic correction term constant

    % Relativistic Effect
    Delta_t_r = F.*e.*sqrt(a).*sin(E_k);

    % Satellite Clock Error
    Delta_t_SV = af0 + af1.*(t_s-toc) + af2.*(t_s-toc).^2 + Delta_t_r;
    
end