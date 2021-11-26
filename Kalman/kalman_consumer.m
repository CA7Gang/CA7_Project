clc; clear; close all;

N = 10000
w = 2
S_meas = [0 0 0 0 0;
	0 0 0 -w 0;
	0 0 0 0 -w/2;
	0 w 0 0 0;
	0 0 w/2 0 0];

B_meas = [1;0;0;0;0];
C_meas = [1, 1, 1, 0, 0];

ts = 1/1000;

S_meas_d = (eye(5) + S_meas.*ts)
abs(eig(S_meas_d))


%%
k0 = 1.2, k1 = 1, k2 = 1, k3 = 1, k4 = 1
var_w = 0.001 % measurement uncertainty
var_z = 0.02 % Model/prediction uncertainty
var_w_guess = var_w*1000 % guess of measurement uncertainty
var_z_guess = var_z/10
% initialisation to allow for wanted phaseshift

x = [k0;
	k1*cos(w*0-pi);
	k2*cos(w/2*0 - pi/2);
	k3*sin(w*0-pi);
	k4*sin(w/2*0-pi/2)]

y_meas = [2]

for i = 2:N
	x(:,i) = S_meas_d*x(:,i-1);
	y_meas(i) = C_meas*x(:,i) + randn(1,1)*sqrt(var_w);

end

% figure()
% plot(y_meas)
% 
% t = 0:ts:N*ts
% y2 = k0 + k1*cos(w*t-pi) + k2*cos(w/2*t-pi/2)
% figure()
% plot(t,y2)

%% kalman
% N = 1000
% w = 1
S = [0 0 0;
	0 0 -w;
	0 w 0]

B = [0; 0; 0]
C = [1, 1, 0]

% ts = 1/100

% sys_ct = ss(S,B,C,0)
% sys_dt = c2d(sys_ct,ts,'zoh')

S_d = (eye(3) + S.*ts)
% B_d = B*ts

Qz = [0 0 0; 0 0 0; 0 0 var_z_guess];
Qw = var_w_guess;

% =================== initialising ===================
x_est(:,1) = [0, 0, 0];
Pn = eye(3);
x_est_m(:,1) = x_est(:,1);


for n=2:N
    % =================== prediction ===================
    x_est_m(:,n) = S_d*x_est(:,n-1);
    Pn_m = S_d*Pn*S_d'+Qz;
	y_est_m(:,n) = C*x_est_m(:,n);
    
    % =================== estimation ===================
    Kn = Pn_m*C' * inv(C*Pn_m*C'+Qw);   
    x_est(:,n) = x_est_m(:,n) + Kn*(y_meas(n)-y_est_m(:,n));
    Pn = (eye(3) - Kn*C)*Pn_m;
%  	Pn_compl = (eye(3) - Kn*C)*Pn_m*(eye(3)-Kn*C)' + Kn*Qv*Kn'
end


figure()
plot(y_meas)
hold on
plot(0:N-1, x_est(1,:)+x_est(2,:),'LineWidth',2)
legend('y_{est} - estimated consumer flow', 'y_{measured} = d_p + d_\tau')
