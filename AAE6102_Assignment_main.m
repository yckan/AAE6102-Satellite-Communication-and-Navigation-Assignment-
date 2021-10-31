clear all; clc; close all

[x_s,y_s,z_s,a,e,af0,af1,af2,E_k,t_s,toc,pr_rcvr] = ...
    solve_SV_pos('/media/pnl/Shared_Data/AAE6102/Assignment/Assignment/Data/rcvr.dat', ...
    '/media/pnl/Shared_Data/AAE6102/Assignment/Assignment/Data/eph.dat'); % Required to change path for reading eph.dat and rcvr.dat

Delta_t_SV = solve_SV_clock_err(a,e,E_k,af0,af1,af2,t_s,toc); % Calculate Satellite clock error

x0 = [-2694685.473;-4293642.366;3857878.924;0.0]; % Initialize user position for Least Square Method
[x_u] = solve_x_user_LS(pr_rcvr,x0,x_s,y_s,z_s,Delta_t_SV); % Solve user position and clock bias by LS method
