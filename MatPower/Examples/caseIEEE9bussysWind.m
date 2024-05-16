function mpc = caseIEEE9bussysWind
%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from the following book:
%   Joe H. Chow (Editor), "Time-Scale Modeling of Dynamic Networks with Applications to Power Systems", 
%                          Springer, Berlin, Heidelberg, 1982, Chapter 4, page 70.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	    0	    0	0	1	1	0	16.5 1	1.05 0.95;
	2	2	0	    0	    0	0	1	1	0	18   1	1.05 0.95;
	3	2	0	    0	    0	0	1	1	0	13.8 1	1.05 0.95;
	4	1	0	    0	    0	0	1	1	0	230	 1	1.05 0.95;
	5	1	1.0*90	30	    0	0	1	1	0	230	 1	1.05 0.95;
	6	1	0	    0	    0	0	1	1	0	230	 1	1.05 0.95;
	7	1	1.0*100	35	    0	0	1	1	0	230	 1	1.05 0.95;
	8	1	0	    0	    0	0	1	1	0	230	 1	1.05 0.95;
	9	1	1.0*125	50	    0	0	1	1	0	230	 1	1.05 0.95;
   10	1	0	    0	    0	0	1	1	0	110	 1	1.05 0.95;
   11	1	0	    0	    0	0	1	1	0	110	 1	1.05 0.95;
   12	2	0	    0       0	0	1	1	0	33	 1	1.05 0.95;
%    12	1	-180	-59.16  0	0	1	1	0	33	 1	1.05 0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax  Qmin	 Vg    Sbase  Status  Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	200	  -200	1.040	100	     1	  250	 10	  0	 0	0	0	0	0	0	0	0	0	0;
	2	200	0	240   -240	1.025	100	     1	  300	 10	  0	 0	0	0	0	0	0	0	0	0	0;
	3	150	0	150	  -150	1.025	100	     1	  270	 10	  0	 0	0	0	0	0	0	0	0	0	0;
    12	180	0	0	     0	1.000	100	     1	  180	180	  0	 0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	     rateA	rateB rateC	  ratio	angle status angmin	angmax
mpc.branch = [
	1	4	0	    0.0576	0	      250	250	  250	   0	 0	   1	 -360	360;
	4	5	0.017	0.092	0.158	  250	250	  250	   0	 0	   1	 -360	360;
	5	6	0.039	0.17	0.358	  150	150	  150	   0	 0	   1	 -360	360;
	3	6	0	    0.0586	0	      300	300	  300	   0	 0	   1	 -360	360;
	6	7	0.0119	0.1008	0.209	  150	150	  150	   0	 0	   1	 -360	360;
	7	8	0.0085	0.072	0.149	  250	250	  250	   0	 0	   1	 -360	360;
	8	2	0	    0.0625	0	      250	250	  250	   0	 0	   1	 -360	360;
	8	9	0.032	0.161	0.306	  250	250	  250	   0	 0	   1	 -360	360;
	9	4	0.01	0.085	0.176	  250	250	  250	   0	 0	   1	 -360	360;
    7	10  0.0016  0.0640  0         200   200   200      1.01	 0     1     -360	360;
    10	11  0.0080  0.0347  0.2179    200   200   200      0     0     1     -360	360;
 	11	12  0.0016  0.0827  0         200   200   200      1.01	 0     1     -360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500	0	3	0.11	5	150;
	2	2000	0	3	0.085	1.2	600;
	2	3000	0	3	0.1225	1	335;
    2	3000	0	3	0       0	  0;
];
