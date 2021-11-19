clc; clear; close all;


% =========================================================================
% ====================== our own linearisation model ======================
% =========================================================================
A=[-0.3236   -0.0406   -0.1577   10.0671   -0.7623    0.0000
   -0.1429   -0.3189    0.2176   -1.3489   -6.0984         0
   -0.0275    0.0687   -0.3968    7.1758    1.5246         0
    0.1089    0.0687    0.0196  -20.2792    1.5246   -0.0000
   -0.0551   -0.0443   -0.0272    1.4425  -15.3466   -0.0999
    0.1486    0.0817   -0.2877   11.1436    8.2419   -0.4174]

B=[0.0982    0.0000
    0.0078         0
    0.1147         0
   -0.0408         0
   -0.0087   -0.0193
   -0.0651   -0.0715]

C = [0 0 1 0 0 0;
	0 0 -1 -1 -1 -1]

D = zeros(2,2)

states = {'q_{c1}' 'q_{c2}' 'd_{p1}' 'd_{c1}' 'd_{c2}' 'd_{t}'};
inputs = {'\omega_1' '\omega_2'};
outputs = {'d_{1}' 'd_{13}'};
sys = ss(A,B,C,D,'statename', states, 'inputname', inputs, 'outputname', outputs)

% RGA = A.*inv(A);

x      = 4;   % Screen position
y      = 3;   % Screen position
width  = 30; % Width of figure
height = 25; % Height of figure (by default in pixels)

% figure(1)
figures = []
figures = figure( 'Color', 'white', 'Units','centimeters','Position', [x y width height])
h = bodeplot(sys)
grid on
setoptions(h,'PhaseVisible','off');
ylim = {[-100 1]; [-100 1]; [-100 1]; [-100 1]}
setoptions(h,'Ylim',ylim);
hold off


% p = getoptions(h)
%%
% collect tf 11, and tf 22
[Num_pump1, Denom_pump1] = ss2tf(A,B,C,D,1);
[Num_pump2, Denom_pump2] = ss2tf(A,B,C,D,2);

G_pump11 = zpk(tf(Num_pump1(1,:),Denom_pump1))          %tf from in1 to out1
G_pump22 = zpk(tf(Num_pump2(2,:),Denom_pump2))          %tf from in2 to out2
[num_Delay,den_Delay] = pade(4,8)
g_delay = tf(num_Delay, den_Delay)


figure(21)
step(G_pump11)
hold on
step(G_pump11*g_delay)

% figure(9)
% ax1 = subplot(211)
% pzmap(G_pump11)
% ax2 = subplot(212)
% linkaxes([ax1 ax2], 'x')
% pzmap(G_pump22)
pole(G_pump11)'
pole(G_pump22)'
zero(G_pump11)'
zero(G_pump22)'
s = tf('s');
G_pump1 = zpk(0.1147 * (s+0.1944)/((s+0.393)*(s+0.1597)))

figure(22)
step(G_pump11,G_pump1)
legend(["True tf", "Approximation"])

figure(23)
bode(G_pump11,G_pump1)
legend(["True tf", "Approximation"])






%%
controlSystemDesigner('rlocus',G_pump1);

%% rlocus figure for worksheets
figures = [figures figure( 'Color', 'white', 'Units','centimeters','Position', [x y width height])]
% PI = 
rlocus(((s+0.05)/s)*G_pump1*pade(4,1))
%%
[num_Delay,den_Delay] = pade(4,10)
g_delay = tf(num_Delay, den_Delay)
% 1/timedelay > BW 
% BW = 0.1. 
% 1/Ti 10 times less than BW. 1/Ti = 0.01
% gain obtained with root locus: 4.5129 (s+0.01)
% PI_slow = 1.8033*(s+0.0497)/s		%bw = 0.05 stable
PI_slow = 1.8033*(s+0.0497)/s		%bw = 0.05 stable
PI_med = 3.0756*(s+0.0497)/s		%bw = 0.1 stable
PI_fast = 4.0756*(s+0.0497)/s		%bw = 0.2 stable but OS


figure(31)
step(feedback(G_pump11*PI_slow,1),feedback(G_pump11*PI_slow*g_delay,1))
hold on
step(feedback(G_pump11*PI_med,1),feedback(G_pump11*PI_med*g_delay,1))
step(feedback(G_pump11*PI_fast,1),feedback(G_pump11*PI_fast*g_delay,1))

figure(32)
margin(G_pump11*PI_slow)
hold on
margin(G_pump11*PI_med)
margin(G_pump11*PI_fast)
legend("Slow","Med","Fast")

%% save figures
savepath = 'C:\Users\kaspe\Documents\Git\Repos\CA7_Writings\CA7_Writings_Worksheets\Pictures'
filename = ["PumpMagPlot.png"]

% filename = []
% for i1=1:length(folders)
%     for i2=1:length(TCs) 
%         temp_fil = append("fig_accep_", TCs(i2), ".png")
%         filename = [filename, temp_fil]
%     end
% end

for i=1:length(filename)
    f = fullfile(savepath,filename(i))

    exportgraphics(figures(i), f, 'Resolution', 400)
end