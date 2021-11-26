% Actual system/consumer behaviour simulation
% ------------------------------------------------------------------------
clc; clear; close all;
% Figure settings
fig_pos_reg		= [1600 500 900 600];
fig_color = 'white';

% ------------------------------------------------------------------------
N = 10000;		% # of samples
ts = 1/1000;	% Sampling time
w = 2;			% omega


% Simulation of actual system/behaviour
% -----------------------------------------
% A-Matrix
S_meas =	[0 0 0 0 0
			0 0 0 -w 0
			0 0 0 0 -w/2
			0 w 0 0 0
			0 0 w/2 0 0];

B_meas = [1;0;0;0;0];
C_meas = [1, 1, 1, 0, 0];

% Discretized actual system
S_meas_d = (eye(5) + S_meas.*ts)
abs(eig(S_meas_d)) % Discrete poles


var_w = 0.02 % measurement uncertainty
var_z = 0.02 % Model/prediction uncertainty
var_w_guess = var_w*1	% Guess of measurement uncertainty
var_z_guess = var_z/10	% Guess of Model/prediction uncertainty

k0 = 1.2; k1 = 1; k2 = 1; % Constants

% state initialisation to allow for wanted phaseshift
x = [k0
	k1*cos(w*0-pi)
	k2*cos(w/2*0 - pi/2)
	k1*sin(w*0-pi)
	k2*sin(w/2*0-pi/2)]

y_meas = [2];

% Running simulation of acutal system w. a flow measurement as output
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
% ------------------------------------------------------------------------

% N = 10000;		% # of samples
% ts = 1/1000;	% Sampling time
% w = 2;			% omega

% Model of system (not the right model)
% -----------------------------------------
S = [0 0 0
	0 0 -w
	0 w 0];

B = [0; 0; 0];
C = [1, 1, 0];

% Discretized model of system
S_d = (eye(3) + S.*ts)
% B_d = B*ts

% Testing that manual discretization is correct
% sys_ct = ss(S,B,C,0)
% sys_dt = c2d(sys_ct,ts,'zoh')


% Model uncertainty covariance matrix 
Qz = [0 0 0; 0 0 0; 0 0 var_z_guess];
% Measurement uncertainty7
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


f = figure()
f.Position = fig_pos_reg;
f.Color = fig_color;
plot(y_meas)
hold on
plot(0:N-1, x_est(1,:)+x_est(2,:),'LineWidth',2)
legend('y_{measured} = d_p + d_\tau', 'y_{est} - estimated consumer flow')
