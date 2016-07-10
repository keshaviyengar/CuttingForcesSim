% Project 2: Simulating End Mill Forces %
% By: Keshav Iyengar %

% m file to simulate helical end milling forces %
clear;

% Angle conversions
d2r = pi/180; % degree to radian
r2d = 180/pi; % radian to degree

calib_factor=200; % force calibration factor [N/V]

% options %
phi_st = 0;
phi_ex = pi/3; % quarter immersion entry and exit angles [rad]
%phi_ex = pi/2; % half immersion entry and exit angles [rad]
%phi_ex = 2*pi/3; % three quarter immersion entry and exit angles [rad]

% Load Data file %
title_text = 'Quarter Immersion - Down-milling';
load DM_4mm_1_4C.dat; tmp=DM_4mm_1_4C; clear DM_4mm_1_4C
% load New_Quarter01.txt; tmp=New_Quarter01.txt; clear New_Quarter01.txt
% OTHER DATA FILES
% title_text = 'Half Immersion - Down-milling';
% load DM_4mm_2_4C.dat;

% title_text = '3/4 Immersion - Down-milling';
% load DM_4mm_3_4C.dat; 

% title_text = 'Full Immersion ';
% load DM_4mm_4_4C.dat; 



% measured forces %
t_meas = tmp(:,1);
Fx_meas = calib_factor*(tmp(:,3));
Fy_meas = calib_factor*(tmp(:,2));
F_meas = sqrt(Fx_meas.*Fx_meas + Fy_meas.*Fy_meas);
N_meas = length(Fx_meas);

% Workpeice data %
% AL 6061 %
tau_s = 209.24; % shear strength [N/mm2]
phi_n = 27.67; % shear angle [rad]
beta_n = 32.27; % Friction angle [rad]
Kte = 23.66; % tangential edge force coefficent [N/mm]
Kfe = 5.57; % normal edge force coefficent [N/mm]

% tool %
% 4 flute HSS cylindrical endmill %
D = 25.4; % diame   ter [mm]
beta = d2r*30.0; % helix angle [rad]
alpha_n = d2r*5.0; % normal rake angle [rad]
Nf = 2; % number of flutes []

% Machining conditions %
n = 1200; % spindle speed [rpm]
a = 4.0; % axis depth of cut [mm]
f = 1200; % feedrate [mm/min]
i = beta; % oblique angle = helix angle [rad]
eta = i; % chip flow angle = oblique angle [rad] [Stabler's chip flow rule]

% Cutting coefficents
den = sqrt(cos(phi_n+beta_n-alpha_n)^2 + tan(eta)^2*sin(beta_n)); % Denominator

% Discritization %
Ts = 0.0001; % Sampling period [sec]
Trev = (n/60); % period of one revolution [sec]

rot_number = 2; % number of revolutions to consider []
t = 0:Ts:rot_number*Trev; % time array
%t_meas = (0:ts:(N_meas-1)*Ts)' % time array for measured forces [sec]
dphi = 2*pi*Ts/Trev; % angular increment
phi = (0:dphi:rot_number*dphi)'; % rotation array for bottom of flute 1 [rad]
Nt = length(phi); % number of time (or rotation) samples to consider

dz = 0.1; % vertical integration increment [mm]
z = (0:dz:a)'; % vertical integration array
Nz = a/dz; % number of vertical segments

phi_p = 2*pi/Nf; % pitch angle [rad]
c = f/(Nf*n); % feed per tooth [mm]

% allocate array for %
Fx = zeros(Nt,1); % x axis force history [N]
Fy = zeros(Nt,1); % y axis force history [N]
Ft = zeros(Nt,1); % tangential force history [N]
F = zeros(Nt,1); % resultant force history (in x-y plane) [N]
T = zeros(Nt,1); % torque history [N*m]
