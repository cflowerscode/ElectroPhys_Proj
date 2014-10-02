




%% Project Part 1 %%%%%
%  Author: Caitlin Flowers
%  email: cflowers7@gatech.edu

%% Constant Variables %%%%
%h = 1/100;      %step value
%tTot = 100;  %total simulation time in milliseconds 
%t = 0:h:100; 
% Conductance values for K, Na, and L respectively mS/cm^2
gK = 36.0;
gNa  = 120.0;
gL = .3;

%Ion potentials (K,Na,L) in mV
eK = -12.0;
eNa = 115.0;
eL = 10.6;
vM = -70.0;      %Resting membrane potential

%% Equations %%%%

%Gating variables
aM0 = .1*((25-vM)/(exp((25-vM)/10)-1));
bM0 = 4*exp(-vM/18);
aN0 = .01*((10-vM)/(exp((10-vM)/10)-1));
bN0 = .125*exp(-vM/80);
aH0 = .07*exp(-vM/20);
bH0 = 1/(exp((30-vM)/10)+1);


%%  Part 1 %%%%
% Simulate a steady state neuron and plot the membrane potential and the
% sodium and potassium conductance values

% Nothing should change

% Initialize arrays 

%V_Mem = zeros(size(t));
%CondK = zeros(size(t));
%CondNa = zeros(size(t));

%V_Mem(1) = vM;

%for n = 1:length(t)
%    V_Mem(n+1) = V_Mem(n)+h*(iIon);
%    CondK(n) = gK;
%    CondNa(n) = gNa;
%end

%Plot Membrane Voltage
%figure 
%plot(t, V_Mem, 'B')
%title('Steady State Membrane Voltage')
%xlabel('Time ms')
%ylabel('Membrane Voltage mV')

%Plot conductance
%figure
%plot(t,CondK,'R',t,CondNa, 'B')
%title('Steady State Conductance Values')
%xlabel('Time ms')
%ylabel('Conductancs mS/cm^2')
%axis([0, tTot, 0, 150])


%% Part 2 %%%%
% Introduce an injection current to help the axon live!

% Initial conditions (?)
m0 = aM0/(aM0+bM0);
n0 = aN0/(aN0+bN0);
h0 = aH0/(aH0+bH0);



iInj = 50*10^-3; %Injection current
iNa = 0;
iK = 0;
iL = 0;
iIon = 0; 

st = 1/10;
t = 0:st:100;

M = zeros(size(t));
M(1) = m0;
N = zeros(size(t));
N(1) = n0;
H = zeros(size(t));
H(1) = h0;

I = ones(size(t))*iInj;
%I(:) = iInj;


V_Mem = zeros(size(t));
V_Mem(1) = vM;
%Use Euler's now

for l = 1:length(t)-1
    %integrate gating variables
    v = V_Mem(l);
    m = M(l);
    n = N(l);
    h = H(l);
    aM = (.1*((25-v)/(exp((25-v)/10)-1)));
    bM = (4*exp(-v/18));
    aN = (.01*((10-v)/(exp((10-v)/10)-1)));
    bN = (.125*exp(-v/80));
    aH = (.07*exp(-v/20));
    bH = (1/(exp((30-v)/10)+1));
    M(l+1) = m + st*((aM*(1-m)-bM*m));
    N(l+1) = n + st*((aN*(1-n)-bN*n));
    H(l+1) = h + st*((aH*(1-h)-bH*h));
    
    %calculate currents
    iNa = (m^3)*gNa*h*(v-eNa);
    iK = (n^4)*gK*(v-eK);
    iL = gL*(v - eL);
    iIon = I(l) - iK - iNa - iL;
    
    V_Mem(l+1) = v + st*(iIon);
end

%Plot Membrane Voltage
figure 
plot(t, V_Mem)
title('Steady State Membrane Voltage')
xlabel('Time ms')
ylabel('Membrane Voltage mV')


    





