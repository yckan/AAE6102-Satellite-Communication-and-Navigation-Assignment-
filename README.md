# AAE6102-Satellite-Communication-and-Navigation-Assignment-
A Matlab Script to estimate user position and clock bias using gps ephemeris and receiver data.

Steps to run the script.
1. Change the path of eph.dat and rcvr.dat in AAE6102_Assignment_main.m
2. Run AAE6102_Assignment_main.m
3. x_u outputs the estiamted user positions and clock bias based on the inputted eph.dat and rcvr.dat files.

This Matlab script consists of 3 MATLAB function script and 1 main script.
1. AAE6102_Assignment_main.m <-- Run to calculate the estimated user position and clock bias
2. solve_SV_pos.m <-- A script to calculate Satellite position in ECEF based on the provided ephemeris
3. solve_SV_clock error.m <-- compute satellite clock error based on provided ephemeris af0, af1, af2, toc, and the relativistic effects.
4. solve_x_user_LS.m <-- estimated user position and clock bias by least square method

All the input and output of the individual scripts are written in the header of the script.

Two input files are required in particular format for solving user position and clock bias.
1. eph.dat <-- data files to supply ephemeris parameter
2. rcvr.dat <-- data files to supply receiver measurements

% Format of eph.dat

eph.dat is a 8 x 24 matrix containing the ephemeris data from a GPS receiver. This data is used to
estimate the orbital position of each satellite at any given time. Each row contains ephemeris data for a
single satellite. The columns of this matrix include the following data:
Column1:rcvr_tow; --receiver time of week(s)
Column 2: svid; -- satellite PRN number (1 – 32)
Column 3: toc; -- reference time of clock parameters (s)
Column 4: toe; -- reference time of ephemeris parameters (s)
Column 5: af0; -- clock correction coefficient – group delay (s)
Column 6: af1; -- clock correction coefficient (s/s)
Column 7: af2; -- clock correction coefficient (s/s/s)
Column 8: ura; -- user range accuracy (m)
Column 9: e; -- eccentricity (-)
Column 10: sqrta; -- square root of semi-major axis a (m**1/2)
Column 11: dn; -- mean motion correction (r/s)
Column 12: m0; -- mean anomaly at reference time (r)
Column 13: w; -- argument of perigee (r)
Column 14: omg0; -- right ascension (r)
Column 15: i0; -- inclination angle at reference time (r)
Column 16: odot; -- rate of right ascension (r/s)
Column 17: idot; -- rate of inclination angle (r/s)
Column 18: cus; -- argument of latitude correction, sine (r)
Column 19: cuc; -- argument of latitude correction, cosine (r)
Column 20: cis; -- inclination correction, sine (r)
Column 21: cic; -- inclination correction, cosine (r)
Column 22: crs; -- radius correction, sine (m)
Column 23: crc; -- radius correction, cosine (m)
Column 24: iod; -- issue of data number

rcvr.dat is an 8x7 matrix containing raw ranging information. Each of the 8 rows contains independent
measurements for each of 8 satellites in view at the current epoch (an epoch is simply a term refers to a
single discrete time; since our receivers provide data at approximately 1 sec. intervals, each epoch
occurs approximately 1 sec. after the prior epoch.
The columns of this matrix in clued the following data:
Column 1: rcvr_tow; -- receiver time of week (s)
Column 2: svid; -- satellite PRN number (1 – 32)
Column 3: pr; -- pseudorange (m)
Column 4: cycles; -- number of accumulated cycles
Column 5: phase; -- to convert to (0 – 359.99) mult. by 360/2048
Column 6: slp_dtct; -- 0 = no cycle slip detected; non 0 = cycle slip
Column 7: snr_dbhz; -- signal to noise ratio (dB-Hz)
