clc
clear all
close all

% Define nominal system and system with modelling error, then discretize

K = 0.0217;
a = 6.87;
Ac = [0 1; 0 -a];
Bc = [0;1];
Cc = [K 0];
Dc = 0;

A_error = [0 1; 0 -2*a];
B_error = [0;0.5];
C_error = [0.4*K 0];


nomsys = ss(Ac,Bc,Cc,Dc);
contsys = ss(A_error,B_error,C_error,Dc);
% contsys = nomsys;

fs = 0.1; Ts = 1/fs;

dNomSys = c2d(nomsys,Ts);
dSys = c2d(contsys,Ts);

n = size(dSys.A,1);
m = size(dSys.B,2);
y = size(dSys.C,1);

x0 = zeros(n+y,1);
ref = 5;


% Construct velocity-form system matrices:

% Av = [dSys.A zeros(n,y); dSys.C eye(y,y)];
% Bv = [dSys.B ; zeros(y,m)];
% Cv = [dSys.C eye(y,y)];

Av = [dSys.A zeros(n,y); dSys.C*dSys.A eye(y,y)];
Bv = [dSys.B ; zeros(y,m)];
Cv = [zeros(y,n) eye(y,y)];

VSys = ss(Av,Bv,Cv,[],Ts);

% Make an LQR gain matrix

Q = Cv'*Cv; % Reference deviation cost
R = 0.01.*eye(m,m); % Actuation cost

[K,P,e] = lqr(VSys,Q,R);

%%


x = zeros(3,1);
dU(1) = 0; % Control input delta
x_real = zeros(n,1) + [1; -3];
refval(1) = 5;
uLQR(1) = 50;
iters = 1000;
yLQR(1) = 0;

for ii = 1:iters
    
   dU(ii+1) = uLQR(ii)+dU(ii);
    
   x = Av*x+Bv*uLQR(ii); 
   yLQR(ii+1) = Cv*x;
 
   x_real(:,ii+1) = dNomSys.A*x_real(:,ii)+dNomSys.B*dU(ii+1);
   y_real(ii+1,1) = dNomSys.C*x_real(:,ii+1);
   x(3) = y_real(ii+1,1)-refval(ii);
   
   uLQR(ii+1) = -K*x;
   
   if (ii < ceil(2*iters/3)) && (ii > ceil(iters/3))
       refval(ii+1) = 2.5;
   elseif ii > ceil(2*iters/3)
       refval(ii+1) = 2.5;
   else
       refval(ii+1) = refval(ii);
   end
   
end

figure(1)
subplot(2,2,1)
ylabel('Output')
xlabel('Samples')
plot(yLQR)
title('Velocity-form system')
legend('Tracking error','interpreter','latex')
subplot(2,2,2)
plot(uLQR)
title('Differential control input')

subplot(2,2,3)
ylabel('Output')
xlabel('Samples')
plot(y_real)
hold on
plot(refval,'--r')
hold off
title('Real System')
legend('Process value','Reference')
subplot(2,2,4)
plot(dU)
title('Full control input')