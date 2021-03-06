clear all; clc;
q       = 1.6e-10;          %1.6e-19;            %electron charge in A/ns
vg      = 8.36e4;           %8.36e9;             %group velocity
d       = 24;               %24e-9;              %thickness of active layer in cm
l       = 1000;             %1e-6;               %length of cavity in VCSEL cm 
ao      = 3.5e-2;           %3.5e-16;            %differential gain in cm2
beta    = 1;                %1;                  %spontaneous emission rate
ntr     = 1.33e-3;             %1.33e18;            %transparency carrier density
nsp     = 2;                %2;                  %inversion factor
gamma_1 = 0.63;             %0.75;               %confinement factor for mode 1
gamma_2 = 0.35 ;            %0.35;               %confinement factor for mode 2
ks      = 8.6e-7;           %8.6e-7;             %gain compression coefficient
r1      = 4000;             %4e-6;               %radius of mode 1
r2      = 7500;             %7.5e-6;             %radius of mode 2
alpha   = 1e-6;             %1000;               %material loss of the active layer
rf      = 0.9991;           %0.9991;             %reflectivity of the front fece
rb      = 0.9998;           %0.9998;             %reflectivity of the back face
taup1   = 2;                %2e-12;              %photon lifetime of mode 1
taup2   = 1.88;             %1.88e-12;           %photon lifetime of mode 2
taue    = 3000;             %3e-9;               %carrier lifetime
H       = 6.626e-22;        %                    %planks constant in nm g and ns

v1=pi*(r1^2)*d;
v2=pi*(r2^2)*d;


% g1l=(vg*(d/l)*ao*(c1-ntr));              % modal gain
% g2l=(vg*(d/l)*ao*(c2-ntr));              % modal gain
% 
% g1=((vg*(d/l)*ao*(c1-ntr))*(1-ks*(gamma_1*p1+gamma_2*p2)));             %linear modal gain
% g2=((vg*(d/l)*ao*(c2-ntr))*(1-ks*((p1*(1-gamma_1))+(p2*(1-gamma_2))))); %linear modal gain



% fc1= @carrier1 %(t, bias,c1,p1,p2) (bias/(v1*q)) - (c1/taue) - ((gamma_1*p1 + gamma_2*p2)*(((vg*ao*(c1-ntr)*(d/l))*(1-ks*(gamma_1*p1+gamma_2*p2)))/v1));
% fc2= @carrier2 %(t, bias,c2,p1,p2) (bias/(v2*q)) - (c2/taue) - (((1-gamma_1)*p1 + (1-gamma_2)*p2)*(((vg*(d/l)*ao*(c2-ntr))*(1-ks*((p1*(1-gamma_1))+(p2*(1-gamma_2)))))/v2));
% fp1= @photon1  %(t, p1,p2,c1,c2) (gamma_1*((vg*(d/l)*ao*(c1-ntr))*(1-ks*(gamma_1*p1+gamma_2*p2))) + (1-gamma_1)*((vg*(d/l)*ao*(c2-ntr))*(1-ks*((p1*(1-gamma_1))+(p2*(1-gamma_2))))) - (1/taup1))*p1 + ((beta*nsp)/taup1);
% fp2= @photon2  %(t, p1,p2,c1,c2) (gamma_2*((vg*(d/l)*ao*(c1-ntr))*(1-ks*(gamma_1*p1+gamma_2*p2))) + (1-gamma_2)*((vg*(d/l)*ao*(c2-ntr))*(1-ks*((p1*(1-gamma_1))+(p2*(1-gamma_2))))) - (1/taup2))*p2 + ((beta*nsp)/taup2);

%step size

h=5.1333e-15;%1e-14;

% modal wavelengths

nu = 30/1510e-7;  % photon frequency mode 1
nu2  = 30/1550e-7;  % photon frequency mode 2

%initial values
bias=[0:0.1:120]*1e-2;
% bias=[0:0.01:3];
% bias =(1+ sin(2*pi*nu*((0.1:0.001:.35)*1e-3))+sin(2*pi*nu2*((0.1:0.001:.35)*1e-3)));
% bias=[0.00E+00 5.00E-02 1.00E-01 1.50E-01 2.00E-01 2.50E-01 3.00E-01 3.50E-01 4.00E-01 4.50E-01 5.00E-01 5.50E-01 6.00E-01 6.50E-01 7.00E-01 7.50E-01 8.00E-01 8.50E-01  9.00E-01 9.50E-01 1.00E+00 1.05E+00 1.10E+00 1.15E+00 1.20E+00 1.25E+00 1.30E+00 1.35E+00 1.40E+00 1.45E+00 1.50E+00 1.55E+00 1.60E+00 1.65E+00 1.70E+00 1.75E+00 1.80E+00 1.85E+00 1.90E+00 1.95E+00 2.00E+00 2.05E+00 2.10E+00 2.15E+00 2.20E+00 2.25E+00 2.30E+00 2.35E+00 2.40E+00];
% bias=0:0.1:2.8;
c1(1)=0;
c2(1)=0;
p1(1)=0;
p2(1)=0;
t(1)=0;
o1(1)=0;
o2(1)=0;

for i=1:length(bias)
    t(i+1)=t(i)+h;
    
   % expression for evaluating k1s
   k1c1= carrier11(t(i), bias(i),c1(i),p1(i),p2(i));
   k1c2= carrier22(t(i), bias(i),c2(i),p1(i),p2(i));
   k1p1= photon11(t(i), p1(i),p2(i),c1(i),c2(i));
   k1p2= photon22(t(i), p1(i),p2(i),c1(i),c2(i));
   
   % expression for evaluating k2s
   k2c1= carrier11(t(i)+h/2, bias(i)        ,c1(i)+h/2+k1c1 ,p1(i)+h/2+k1p1 ,p2(i)+h/2+k1p2);
   k2c2= carrier22(t(i)+h/2, bias(i)        ,c2(i)+h/2+k1c2 ,p1(i)+h/2+k1p1 ,p2(i)+h/2+k1p2);
   k2p1= photon11(t(i)+h/2 , p1(i)+h/2+k1p1 ,p2(i)+h/2+k1p2 ,c1(i)+h/2+k1c1 ,c2(i)+h/2+k1c2);
   k2p2= photon22(t(i)+h/2 , p1(i)+h/2+k1p1 ,p2(i)+h/2+k1p2 ,c1(i)+h/2+k1c1 ,c2(i)+h/2+k1c2);
   
   % expression for k3s
   k3c1= carrier11(t(i)+h/2, bias(i), c1(i)+h/2+k2c1   ,p1(i)+h/2+k2p1   ,p2(i)+h/2+k2p2);
   k3c2= carrier22(t(i)+h/2, bias(i), c2(i)+h/2+k2c2   ,p1(i)+h/2+k2p1   ,p2(i)+h/2+k2p2);
   k3p1= photon11(t(i)+h/2          , p1(i)+h/2+k2p1   ,p2(i)+h/2+k2p2   ,c1(i)+h/2+k2c1  ,c2(i)+h/2+k2c2);
   k3p2= photon22(t(i)+h/2          , p1(i)+h/2+k2p1   ,p2(i)+h/2+k2p2   ,c1(i)+h/2+k2c1  ,c2(i)+h/2+k2c2);
  
   % expression for k4s
   k4c1= carrier11(t(i)+h, bias(i), c1(i)+h+k3c1   ,p1(i)+h+k2p1   ,p2(i)+h+k2p2);
   k4c2= carrier22(t(i)+h, bias(i), c2(i)+h+k3c2   ,p1(i)+h+k2p1   ,p2(i)+h+k2p2);
   k4p1= photon11(t(i)+h          , p1(i)+h+k3p1   ,p2(i)+h+k2p2   ,c1(i)+h+k2c1  ,c2(i)+h+k2c2);
   k4p2= photon22(t(i)+h          , p1(i)+h+k3p1   ,p2(i)+h+k2p2   ,c1(i)+h+k2c1  ,c2(i)+h+k2c2);
   
   c1(i+1)= c1(i) + (h*(k1c1+2*k2c1+2*k3c1*k4c1))/6;
   c2(i+1)= c2(i) + h/6 *(k1c2+2*k2c2+2*k3c2*k4c2);
   p1(i+1)= p1(i) + h/6 *(k1p1+2*k2p1+2*k3p1*k4p1);
   p2(i+1)= p2(i) + h/6 *(k1p2+2*k2p2+2*k3p2*k4p2);
   
%    ps1(i)=sqrt(0.1*nu)*(abs(fft(p1(i))));
%    ps2(i)=sqrt(0.1*nu2)*(abs(fft(p2(i))));
   
   
   o1(i+1)=p1(i)*H*nu*(rb-rf)/taup1;
   o2(i+1)=p2(i)*H*nu2*(rb-rf)/taup2;
   
   
   newo1(i+1)=(p1(i)*H*nu*(1/2*l)*log(1/rf)) /(alpha + (1/2*l)*log(1/(rf*rb))*taue);
end
for i=1:length(p1)
fs1(i)=((-1+2*rand))*sqrt((2*beta*nsp*p1(i)*0.1)/(0.1*taup1));
np1(i)=p1(i)+fs1(i);
fs2(i)=((-1+2*rand))*sqrt((2*beta*nsp*p2(i)*0.1)/(0.1*taup2));
np2(i)=p2(i)+fs2(i);
end

for i=1:length(p1)
    no1(i)=np1(i)*H*nu*(rb-rf)/taup1;
    no2(i)=np2(i)*H*nu2*(rb-rf)/taup1;
end
figure;
plot(t,p1,'displayname','Photon LP01');legend('-dynamiclegend');xlabel('time (in sec)');
ylabel('Photon Concentration');title('Photon Concentration in mode LP01');
figure;
plot(t,p2,'displayname','Photon LP11');legend('-dynamiclegend');xlabel('time (in sec)');
ylabel('Photon Concentration');title('Photon Concentration in mode LP11');
figure;
plot(t,o1,'displayname','Output Power LP01');legend('-dynamiclegend');xlabel('time (in sec)');
ylabel('Power (in W)');title('Output Power in mode LP01');
hold on;
% plot(t,no1,'displayname','Distorted LP01');legend('-dynamiclegend');
figure;
plot(t,o2,'displayname','Output Power LP11');hold on;legend('-dynamiclegend');xlabel('time (in sec)');
ylabel('Power (in W)');title('Output Power in mode LP11');
% plot(t,no2,'displayname','Distorted LP11');legend('-dynamiclegend');
