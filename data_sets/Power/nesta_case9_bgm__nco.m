%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                  %%%%%
%%%%      NICTA Energy System Test Case Archive (NESTA) - v0.7.0      %%%%%
%%%%           Optimal Power Flow - Nonconvex Optimization            %%%%%
%%%%                         05 - June - 2017                         %%%%%
%%%%                                                                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nine bus case exhibiting local minima.
%
%   Bukhsh, W. A. & Grothey, A. & McKinnon,  K. I. M. & Trodden, P. A. 
%   "Local Solutions of the Optimal Power Flow Problem"
%   IEEE Transactions on Power Systems, August, 2013
%
function mpc = nesta_case9_bgm__nco
mpc.version = '2';
mpc.baseMVA = 100.0;

%% area data
%	area	refbus
mpc.areas = [
	1	 5;
];

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 0.0	 0.0	 0.0	 0.0	 1	    0.90950	   -0.00000	 345.0	 1	    1.10000	    0.90000;
	2	 2	 0.0	 0.0	 0.0	 0.0	 1	    0.92181	   12.36699	 345.0	 1	    1.10000	    0.90000;
	3	 2	 0.0	 0.0	 0.0	 0.0	 1	    0.93880	    7.00645	 345.0	 1	    1.10000	    0.90000;
	4	 1	 0.0	 0.0	 0.0	 0.0	 1	    0.91268	   -0.39758	 345.0	 1	    1.10000	    0.90000;
	5	 1	 54.0	 18.0	 0.0	 0.0	 1	    0.91605	   -0.73239	 345.0	 1	    1.10000	    0.90000;
	6	 1	 0.0	 0.0	 0.0	 0.0	 1	    0.94259	    4.84215	 345.0	 1	    1.10000	    0.90000;
	7	 1	 60.0	 21.0	 0.0	 0.0	 1	    0.92843	    4.51554	 345.0	 1	    1.10000	    0.90000;
	8	 1	 0.0	 0.0	 0.0	 0.0	 1	    0.92910	    7.11774	 345.0	 1	    1.10000	    0.90000;
	9	 1	 75.0	 30.0	 0.0	 0.0	 1	    0.90000	   -0.63065	 345.0	 1	    1.10000	    0.90000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	 10.0	 -5.0	 300.0	 -5.0	 0.9095	 100.0	 1	 250.0	 10.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	2	 125.369	 -5.0	 300.0	 -5.0	 0.92181	 100.0	 1	 300.0	 10.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	3	 57.028	 -5.0	 300.0	 -5.0	 0.9388	 100.0	 1	 270.0	 10.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 1500.0	 0.0	 3	   0.110000	   5.000000	 150.000000;
	2	 2000.0	 0.0	 3	   0.085000	   1.200000	 600.000000;
	2	 3000.0	 0.0	 3	   0.122500	   1.000000	 335.000000;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 4	 0.0	 0.0576	 0.0	 250.0	 250.0	 250.0	 0.0	 0.0	 1	 -30.0	 30.0;
	4	 5	 0.017	 0.092	 0.158	 250.0	 250.0	 250.0	 0.0	 0.0	 1	 -30.0	 30.0;
	5	 6	 0.039	 0.17	 0.358	 150.0	 150.0	 150.0	 0.0	 0.0	 1	 -30.0	 30.0;
	3	 6	 0.0	 0.0586	 0.0	 300.0	 300.0	 300.0	 0.0	 0.0	 1	 -30.0	 30.0;
	6	 7	 0.0119	 0.1008	 0.209	 150.0	 150.0	 150.0	 0.0	 0.0	 1	 -30.0	 30.0;
	7	 8	 0.0085	 0.072	 0.149	 250.0	 250.0	 250.0	 0.0	 0.0	 1	 -30.0	 30.0;
	8	 2	 0.0	 0.0625	 0.0	 250.0	 250.0	 250.0	 0.0	 0.0	 1	 -30.0	 30.0;
	8	 9	 0.032	 0.161	 0.306	 250.0	 250.0	 250.0	 0.0	 0.0	 1	 -30.0	 30.0;
	9	 4	 0.01	 0.085	 0.176	 250.0	 250.0	 250.0	 0.0	 0.0	 1	 -30.0	 30.0;
];

% INFO    : === Translation Options ===
% INFO    : Phase Angle Bound:           30.0 (deg.)
% INFO    : AC OPF Solution File:        nesta_case9_bgm__nco.m.opf.sol
% INFO    : 
% INFO    : === Voltage Setpoint Replacement Notes ===
% INFO    : Bus 1	: V=1.0, theta=0.0 -> V=0.9095, theta=-0.0
% INFO    : Bus 2	: V=1.0, theta=0.0 -> V=0.92181, theta=12.36699
% INFO    : Bus 3	: V=1.0, theta=0.0 -> V=0.9388, theta=7.00645
% INFO    : Bus 4	: V=1.0, theta=0.0 -> V=0.91268, theta=-0.39758
% INFO    : Bus 5	: V=1.0, theta=0.0 -> V=0.91605, theta=-0.73239
% INFO    : Bus 6	: V=1.0, theta=0.0 -> V=0.94259, theta=4.84215
% INFO    : Bus 7	: V=1.0, theta=0.0 -> V=0.92843, theta=4.51554
% INFO    : Bus 8	: V=1.0, theta=0.0 -> V=0.9291, theta=7.11774
% INFO    : Bus 9	: V=1.0, theta=0.0 -> V=0.9, theta=-0.63065
% INFO    : 
% INFO    : === Generator Setpoint Replacement Notes ===
% INFO    : Gen at bus 1	: Pg=0.0, Qg=0.0 -> Pg=10.0, Qg=-5.0
% INFO    : Gen at bus 1	: Vg=1.0 -> Vg=0.9095
% INFO    : Gen at bus 2	: Pg=163.0, Qg=0.0 -> Pg=125.369, Qg=-5.0
% INFO    : Gen at bus 2	: Vg=1.0 -> Vg=0.92181
% INFO    : Gen at bus 3	: Pg=85.0, Qg=0.0 -> Pg=57.028, Qg=-5.0
% INFO    : Gen at bus 3	: Vg=1.0 -> Vg=0.9388
% INFO    : 
% INFO    : === Writing Matpower Case File Notes ===
