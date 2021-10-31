% Input rcvr_path - path to rcvr.dat file
% Input eph_path - path to eph.dat file

% Output x_s - Satellite x-position in ECEF
% Output y_s - Satellite y-position in ECEF
% Output z_s - Satellite z-position in ECEF
% Output af0 - Satellite clock bias (sec)
% Output af1 - Satellite clock drift (sec/sec)
% Output af2 - Satellite frequency drift (i.e. aging) (sec/sec^2)
% Output E_k - Eccentric Anomaly
% Output t_s - Satellite time of transmission
% Output toc - reference time of clock paramters (sec)
% Output pr_rcvr - pseudorange from receiver (m)

function [x_s,y_s,z_s,a,e,af0,af1,af2,E_k,t_s,toc,pr_rcvr] = solve_SV_pos(rcvr_path,eph_path)

    % Import data from rcvr.dat and eph.dat
    rcvr=importdata(rcvr_path);
    eph=importdata(eph_path);
    
    rcvr = sortrows(rcvr,2); % sort the rcvr matrix by SVID in ascending order
    eph = sortrows(eph,2); % sort the eph matrix by SVID in ascending order
    
    % Read rcvr.dat data into associated entries
    rcvr_tow_rcvr = rcvr(:,1); % receiver time of week (s)
    svid_rcvr = rcvr(:,2); % satellite PRN number (1-32)
    pr_rcvr = rcvr(:,3); % pseudorange (m)
    cycles_rcvr = rcvr(:,4); % number of accumulated cycles
    phase_rcvr = rcvr(:,5); % to convert to (0-359.99) mult. by 360/2048
    slp_dtct_rcvr = rcvr(:,6); % 0 = no cycle slip detected; non 0 = cycle slip
    snr_dbhz_rcvr = rcvr(:,7); % signal to noise ratio (db-Hz)
    
    % Define Constant
    omega_e_dot = 7.2921151467e-5; % WGS 84 value of earth's universal gravitation constant (m^3/s^2)
    mu = 3.986005e14; % WGS 84 value of earth's rotatino rate (r/s)
    c = 299792458.0; % speed of light (m/s)
    F = -4.442807633e-10; %Relativistic correction term constant
    % Calculate satellite time of transmission by geometric range
    t_s = rcvr_tow_rcvr-pr_rcvr./c; % satellite time of transmission
    
    % Read eph.dat data into associated entries
    rcvr_tow = rcvr(:,1); % receiver time of week (s)
    svid = eph(:,2); % saatellite PRN number (1-32)
    toc = eph(:,3); % reference time of clock paramters (s)
    toe = eph(:,4); % reference time of ephemeris parameters (s)
    af0 = eph(:,5); % clock correction coefficient-group delay (s)
    af1 = eph(:,6); % clock correction coefficient (s/s)
    af2 = eph(:,7); % clock correction coefficient (s/s/s)
    ura = eph(:,8); % user range accuracy (m)
    e = eph(:,9); % eccentricity (-)
    sqrta = eph(:,10); % square root of semi-major axis a (m**1/2)
    dn = eph(:,11); % mean motion correction (r/s)
    m0 = eph(:,12); % mean anomaly at reference time (r)
    w = eph(:,13); % argument of perigee (r)
    omg0 = eph(:,14); % right ascension (r)
    i0 = eph(:,15); % inclination angle at reference time (r)
    odot = eph(:,16); % rate of right ascension (r/s)
    idot = eph(:,17); % rate of inclination angle (r/s)
    cus = eph(:,18); % argument oflatitude correction, sine (r)
    cuc = eph(:,19); % argument of latitude correction, cosine (r)
    cis = eph(:,20); % inclination correction, sine (r)
    cic = eph(:,21); % inclination correction, cosine (r)
    crs = eph(:,22); % radius correction, sine (m)
    crc = eph(:,23); % radius correction, cosine (m)
    iod = eph(:,24); % issue of data number
    
    a= sqrta.^2; % semimajor axis
    n = sqrt(mu./a.^3)+dn; % corrected mean motion
    t_k = t_s-toe; % time from ephemeris epoch
    M_k = m0 + n.*t_k; % mean anomaly
    E_k0 = M_k; % initialize eccentric anomaly
    
    % Successive Method to solve eccentric anomaly
    for i = 1:size(eph,1)
        eps = 1;
        n=1;
        while eps>=1e-5
            E_k(i,1) = e(i)*sin(E_k0(i,1))+M_k(i,1);
            eps=abs(E_k(i,1)-E_k0(i,1));
            E_k0(i,1) = E_k(i,1);
            n = n+1;
        end
    end
    
    % True Anomaly
    nu = 2*atan(sqrt((1+e)./(1-e)).*tan(E_k/2)); % true anomaly
    
    phi_k = nu+w; % Argument of latitude
    del_phi_k = cus.*sin(2*phi_k)+cuc.*cos(2*phi_k); % Argument of latitude correction
    del_r_k = crs.*sin(2*phi_k)+crc.*cos(2*phi_k); % Radius correction
    del_i_k = cis.*sin(2*phi_k)+cic.*cos(2*phi_k); % Inclinationcorrection
    u_k = phi_k + del_phi_k; % Corrected argument of latitude
    r_k = a.*(1-e.*cos(E_k))+del_r_k; % Corrected radius
    i_k = i0 + idot.*t_k+del_i_k; % Corrected inclination
    omega_k = omg0 + (odot-omega_e_dot).*t_k-omega_e_dot*toe; %Corrected longitude of the ascending node
    x_p = r_k.*cos(u_k); % In-plane x position
    y_p = r_k.*sin(u_k); % In-plane y position
    x_s = x_p.*cos(omega_k) - y_p.*cos(i_k).*sin(omega_k); % ECEF x-coordinate
    y_s = x_p.*sin(omega_k) + y_p.*cos(i_k).*cos(omega_k); % ECEF y-coordinate
    z_s = y_p.*sin(i_k); % ECEF z-coordinate
end