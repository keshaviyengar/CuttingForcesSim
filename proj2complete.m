% Project 2: Simulating End Mill Forces %
% Quarter Immersion %
% By: Keshav Iyengar: 20466956 %

% m file to simulate helical end milling forces %
clear;

% Angle conversions
d2r = pi/180; % degree to radian
r2d = 180/pi; % radian to degree

% Workpeice data %
% AL 6061 %
tau_s = 209.24; % shear strength [N/mm2]
phi_c = 27.67*d2r; % shear angle [rad]
B_a = 32.27*d2r; % Friction angle [rad]
kte = 23.66; % tangential edge force coefficent [N/mm]
kfe = 5.57; % normal edge force coefficent [N/mm]

% 4 flute HSS cylindrical endmill %
D = 25.4; % diameter [mm]
r = D/2; % radius [mm]
B = d2r*30.0; % helix angle [rad]
a_n = d2r*5.0; % normal rake angle [rad]
Nf = 2; % number of flutes []

% Machining conditions %
n = 1200; % spindle speed [rpm]
a = 4.0; % axis depth of cut [mm]
f = 1200; % feedrate [mm/min]
calib_factor = 200; % force calibration factor [N/V]

% Assumptions %
i = B; % oblique angle = helix angle [rad]
eta = i; % chip flow angle = oblique angle [rad] [Stabler's chip flow rule]
B_n=atan(tan(B_a)*cos(eta)); %[rad] Friction Coefficient Ba and Shear Stress, Ts, are the same in both orthogonal & Oblique Cutting
phi_n = phi_c; %[rad] Orthogonal Shear Angle = Normal Shear Angle in Oblique
a_r = a_n; %[rad] Orthogonal Normal Rake Angle = Rake Angle in Oblique 
kre = kfe; %[N/mm] Normal Edge Force Coefficient in Orthogonal Cutting is Equivalent to Radial Edge Force in Milling
% Feed Per Tooth
c = f/(Nf*n); % feed per tooth [mm]

% Cutting coefficents
den = sqrt(cos(phi_n+B_n-a_n)^2 + tan(eta)^2*sin(B_n)^2); % Denominator
ktc = (tau_s/sin(phi_c))*(cos(B_n-a_n)+tan(i)*tan(eta)*sin(B_n))/(den);
krc = (tau_s/sin(phi_c))*(cos(B_n-a_n)*tan(i)-tan(eta)*sin(B_n))/(den);
kfc = (tau_s/(sin(phi_c)*cos(i)))*((sin(B_n-a_n))/(den));

% Discritization %
Ts = 0.0001; % Sampling period [sec]
Trev = (n/60); % period of one revolution [sec]
rot_number = 4; % number of revolutions to consider []
dphi = Trev*2*pi*Ts; % angular increment
phi = (0:dphi:rot_number*2*pi); % rotation array for bottom of flute 1 [rad]
Nt = length(phi); % number of time (or rotation) samples to consider
t = (0:Ts:rot_number/Trev); % time array
%t_meas = (0:ts:(N_meas-1)*Ts)' % time array for measured forces [sec]

phi_p = (2*pi)/Nf; % Angle between teeth [rad]

% Vertical Integration %
dz = 0.1; % vertical integration increment [mm]
Nz = a/dz; % number of vertical segments
z = (dz:dz:a); % [mm] vertical integration array

% Immersion %
im = [0.25,0.5,0.75,1]*D; %[mm] Downmilling 1/4,1/2,3/4,1 

% allocate array; columns are each immersion%
Fx = zeros(length(im),Nt); % x axis force history [N]
Fy = zeros(length(im),Nt); % y axis force history [N]
Ft = zeros(length(im),Nt); % t axis force history [N]
F = zeros(length(im),Nt); % resultant force history (in x-y plane) [N] 
T = zeros(length(im),Nt); % torque history [N*m]
P = zeros(length(im),Nt); % Power history

% Main Loop for Simulation %
for count=1:length(im)
    phi_st = pi - (acos((r-im(count))/r)); % [rad] start angle
    phi_ex = pi; %[rad] exit angle
    for kt=1:Nt % Rotation counter
        for kz=1:Nz
            for kf=1:Nf
                phi_bottom = phi(kt) + phi_p*(kz-1); % Immersion angle for tooth
                si = ((2*tan(B)*z(kf))/D); % Correction factor for helx
                phi_cur = phi_bottom-si; % Current phi

                % Ensure current angle between 0 and 2pi
                if phi_cur>=0
                    phi_cur2=phi_cur-2*pi*floor(phi_cur/(2*pi));
                else
                    phi_cur2=phi_cur+2*pi*floor(-phi_cur/(2*pi)+1);
                end

                if (phi_cur2>=phi_st) && (phi_cur2 <=phi_ex)
                    h = c*sin(phi_cur2);

                    dFt = dz * (ktc * h + kte); % tangential force contribution [N]
                    dFr = dz * (krc * h + kre); % radial force contribution [N]
                    dFx = -dFt*cos(phi_cur2) - dFr*sin(phi_cur2); % x axis force contribution [N]
                    dFy = dFt*sin(phi_cur2) - dFr*cos(phi_cur2); % y axis force contribution [N]
                    % (dF dFR dFx dFy)

                    Ft(count,kt) = Ft(count,kt) + dFt; % Intregrate tangential force contribution
                    Fx(count,kt) = Fx(count,kt) + dFx; % Intregrate x axis force contribution
                    Fy(count,kt) = Fy(count,kt) + dFy; % Intregrate y axis force contribution
                end
            end
        end
        F(count,kt) = sqrt(Fx(count,kt)^2 + Fy(count,kt)^2); % resultant Cutting forces history [N]
        T(count,kt) = 1e-3*(D/2)*Ft(count,kt); % cutting torque history [N*m]
        P(count,kt) = T(count,kt)*n*2*pi/60; % Power [W]
    end
end

% Load Data file %
load DM_4mm_1_4C.dat; tmp1=DM_4mm_1_4C; clear DM_4mm_1_4C;
t_meas1 = tmp1(:,1) + 0.10; % Phase shift
Fx_meas1 = calib_factor*tmp1(:,3);
Fy_meas1 = calib_factor*tmp1(:,2);
F_meas1 = sqrt(Fx_meas1.*Fx_meas1 + Fy_meas1.*Fy_meas1); 
N_meas1 = length(Fx_meas1);

load DM_4mm_2_4C.dat; tmp2=DM_4mm_2_4C; clear DM_4mm_2_4C;
t_meas2 = tmp2(:,1) + 0.11015; % Phase shift
Fx_meas2 = calib_factor*tmp2(:,3);
Fy_meas2 = calib_factor*tmp2(:,2);
F_meas2 = sqrt(Fx_meas2.*Fx_meas2 + Fy_meas2.*Fy_meas2); 
N_meas2 = length(Fx_meas2);

load DM_4mm_3_4C.dat; tmp3=DM_4mm_3_4C; clear DM_4mm_3_4C;
t_meas3 = tmp3(:,1) + 0.1145; % Phase shift
Fx_meas3 = calib_factor*tmp3(:,3);
Fy_meas3 = calib_factor*tmp3(:,2);
F_meas3 = sqrt(Fx_meas3.*Fx_meas3 + Fy_meas3.*Fy_meas3); 
N_meas3 = length(Fx_meas3);

load DM_4mm_4_4C.dat; tmp4=DM_4mm_4_4C; clear DM_4mm_4_4C;
t_meas4 = tmp4(:,1) + 0.1041; % Phase shift
Fx_meas4 = calib_factor*tmp4(:,3);
Fy_meas4 = calib_factor*tmp4(:,2);
F_meas4 = sqrt(Fx_meas4.*Fx_meas4 + Fy_meas4.*Fy_meas4); 
N_meas4 = length(Fx_meas4);

% Plots Measured and Simulated %
% 1/4 Immersion Forces%
figure (1)  
subplot(2,1,1) 
plot(t,Fx(1,:),'r',t_meas1,Fx_meas1,'b--')
title('1/4 Immersion Measured and Predicted Forces'); 
ylabel('Fx [N]'); 
axis tight
subplot(2,1,2); 
plot(t,Fy(1,:),'r',t_meas1,Fy_meas1,'b--')
ylabel('Fy [N]'); 
xlabel('Time [sec]'); 
axis tight
legend ('Simulated','Measured')

%1/2 Immersion Forces %
figure(2)
subplot(2,1,1) 
plot(t,Fx(2,:),'r',t_meas2,Fx_meas2,'b--')
title('1/2 Immersion Measured and Predicted Forces'); 
ylabel('Fx [N]'); 
axis tight
subplot(2,1,2); 
plot(t,Fy(1,:),'r',t_meas1,Fy_meas1,'b--')
ylabel('Fy [N]'); 
xlabel('Time [sec]'); 
axis tight
legend ('Simulated','Measured')

%3/4 Immersion Forces %
figure(3)
subplot(2,1,1) 
plot(t,Fx(3,:),'r',t_meas3,Fx_meas3,'b--')
title('3/4 Immersion Measured and Predicted Forces'); 
ylabel('Fx [N]'); 
axis tight
subplot(2,1,2); 
plot(t,Fy(3,:),'r',t_meas3,Fy_meas3,'b--')
ylabel('Fy [N]'); 
xlabel('Time [sec]'); 
axis tight
legend ('Simulated','Measured')

%Full Immersion Forces %
figure(4)
subplot(2,1,1) 
plot(t,Fx(4,:),'r',t_meas4,Fx_meas4,'b--')
title('Full Immersion Measured and Predicted Forces'); 
ylabel('Fx [N]'); 
axis tight
subplot(2,1,2); 
plot(t,Fy(4,:),'r',t_meas4,Fy_meas4,'b--')
ylabel('Fy [N]'); 
xlabel('Time [sec]'); 
axis tight
legend ('Simulated','Measured')

% Plot Total Force %
figure (5)
subplot(4,1,1)
plot(t,F(1,:),'r', t_meas1,F_meas1,'r--')
title('Total Force Vs. Time')
legend ('1/4 Imm. Simulated', '1/4 Imm. Measured')
ylabel('Force [N]')
subplot(4,1,2)
plot(t,F(2,:),'b',(t_meas2+0.0025),F_meas2,'b--')
ylabel('Force [N]')
legend ('1/2 Imm. Simulated', '1/2 Imm. Measured')
subplot(4,1,3)
plot(t,F(3,:),'g',t_meas3,F_meas3,'g--')
ylabel('Force [N]')
legend ('3/4 Imm. Simulated', '3/4 Imm. Measured')
subplot(4,1,4)
plot(t,F(4,:),'k',t_meas4,F_meas4,'k--')
ylabel('Force [N]')
xlabel('Time [s]')
legend ('Full Imm. Simulated', 'Full Imm. Measured')


% Plot Torque %
figure (6)
subplot(4,1,1)
plot(t,T(1,:),'r')
title('Simulated Torque Vs. Time')
legend ('1/4 Immersion')
ylabel('Torque [Nm]')
subplot(4,1,2)
plot(t,T(2,:),'b')
ylabel('Torque [Nm]')
legend ('1/2 Immersion')
subplot(4,1,3)
plot(t,T(3,:),'g')
ylabel('Torque [Nm]')
legend ('3/4 Immersion')
subplot(4,1,4)
plot(t,T(4,:),'k')
ylabel('Torque [Nm]')
xlabel('Time [s]')
legend ('Full Immersion')
 
% Plot Power %
figure (7)
subplot(4,1,1)
plot(t,P(1,:),'r')
title('Simulated Power Vs. Time')
legend ('1/4 Immersion')
ylabel('Power [W]')
axis tight
subplot(4,1,2)
plot(t,P(2,:),'b')
ylabel('Power [W]')
legend ('1/2 Immersion')
axis tight
subplot(4,1,3)
plot(t,P(3,:),'g')
ylabel('Power [W]')
legend ('3/4 Immersion')
axis tight
subplot(4,1,4)
plot(t,P(4,:),'k')
ylabel('Power [W]')
xlabel('Time [s]')
legend ('Full Immersion')
axis tight
