% Project 2: Simulating End Mill Forces %
% Half Immersion %
% By: Keshav Iyengar %

% m file to simulate helical end milling forces %
clear;

% Angle conversions
d2r = pi/180; % degree to radian
r2d = 180/pi; % radian to degree

calib_factor = 200; % force calibration factor [N/V]

% options %
phi_st = 0;
%phi_ex = pi/3; % quarter immersion entry and exit angles [rad]
phi_ex = pi/2; % half immersion entry and exit angles [rad]
%phi_ex = 2*pi/3; % three quarter immersion entry and exit angles [rad]
%phi_ex = pi; % full immersion entry and exit angles [rad]

% Load Data file %
% title_text = 'Quarter Immersion - Down-milling';
% load DM_4mm_1_4C.dat; tmp=DM_4mm_1_4C; clear DM_4mm_1_4C
% load New_Quarter01.txt; tmp=New_Quarter01.txt; clear New_Quarter01.txt
% OTHER DATA FILES
title_text = 'Half Immersion - Down-milling';
load DM_4mm_2_4C.dat; tmp=DM_4mm_2_4C; clear DM_4mm_2_4C

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
kte = 23.66; % tangential edge force coefficent [N/mm]
kfe = 5.57; % normal edge force coefficent [N/mm]

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
ktc = tau_s*(cos(beta_n-alpha_n) + tan(i)*tan(eta)*sin(beta_n))/(sin(phi_n)*den);
kfc = tau_s*sin(beta_n-alpha_n)/(sin(phi_n)*cos(i)*den);
krc = tau_s*(cos(beta_n-alpha_n)*tan(i)-tan(eta)*sin(beta_n))/(sin(phi_n)*den);

% Discritization %
Ts = 0.0001; % Sampling period [sec]
Trev = (n/60); % period of one revolution [sec]

rot_number = 2; % number of revolutions to consider []
t = (0:Ts:rot_number*Trev)'; % time array
%t_meas = (0:ts:(N_meas-1)*Ts)' % time array for measured forces [sec]
dphi = 2*pi*Ts/Trev; % angular increment
phi = (0:dphi:rot_number*2*pi)'; % rotation array for bottom of flute 1 [rad]
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

for kt = 1:Nt, % rotation counter
    
    for kz = 1:Nz, % axial depth counter
        
        for kf = 1:Nf, % flute counter
            
            % Immersion angle of current flute at current height along
            % cutter
            phi_cur = phi(kt) - phi_p*(kf-1) - 2*tan(beta)*z(kz)/D;
            
            % express current angle between 0 and 2*pi. Simple mapping.
            if phi_cur>=0,
                phi_cur2 = phi_cur - 2*pi*floor(phi_cur/(2*pi));
            else
                phi_cur2 = phi_cur - 2*pi*floor(-phi_cur/(2*pi));
            end
            
            % consider contribution to cutting process if current position
            % of tool is cutting.
            if (phi_cur2 >= phi_st) && (phi_cur2 <= phi_ex),
                h = c*sin(phi_cur2);
                dFt = dz * (ktc * h + kte); % tangential force contribution [N]
                dFr = dz * (krc * h + kte); % radial force contribution [N]
                dFx = -dFt*cos(phi_cur2) - dFr*sin(phi_cur2); % x axis force contribution [N]
                dFy = dFt*sin(phi_cur2) - dFr*cos(phi_cur2); % y axis force contribution [N]
                % (dF dFR dFx dFy)
                
                Ft(kt) = Ft(kt) + dFt; % Intregrate tangential force contribution
                Fx(kt) = Fx(kt) + dFx; % Intregrate x axis force contribution
                Fy(kt) = Fy(kt) + dFy; % Intregrate y axis force contribution
            end
         % tool counter loop terminates
        end
        % axis depth loop terminates 
    end
    
    F(kt) = sqrt(Fx(kt)^2 + Fy(kt)^2); % resultant Cutting forces history [N]
    T(kt) = 1e-3*(D/2)*Ft(kt); % cutting torque history [N*m]
end

%plot results%
figure(1);clf;zoom on;
subplot(2,1,1); plot(t,Fx,'r',t_meas,Fx_meas,'r--'); title('Measured and Predicted Forces'); ylabel('Fx [N]'); xlabel('Time [s]'); tmp = axis; axis([0 max(t) tmp(3) tmp(4)])
subplot(2,1,2); plot(t,Fy,'b',t_meas,Fy_meas,'b--'); title('Measured and Predicted Forces'); ylabel('Fy [N]'); xlabel('Time [s]'); tmp = axis; axis([0 max(t) tmp(3) tmp(4)])

figure(2);clf;zoom on;
subplot(2,1,1); plot(t,F,'k',t_meas,F_meas,'k--'); title('Resultant Force and Torque'); ylabel('F [N]'); xlabel('Time [s]');
tmp = axis; axis([0 max(t) tmp(3) tmp(4)]);
subplot(2,1,2); plot(t,T,'k'); title('Resultant Force and Torque'); ylabel('T [Nm]'); xlabel('Time [s]');
tmp = axis; axis([0 max(t) tmp(3) tmp(4)]);

