function p1=photon11(t, p1,p2,c1,c2)
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
eps11   =2e-7;              %self-gain saturation coefficient mode 1
eps22   =2e-7;              %self-gain saturation coefficient mode 2
eps12   =0.5e-7;            %cross- gain saturation coefficient mode 1
eps21   =0.5e-7;            %cross- gain saturation coefficient mode 2


v1=pi*(r1^2)*d;
v2=pi*(r2^2)*d;

g1l=(vg*(d/l)*ao*(c1-ntr));              % modal gain
g2l=(vg*(d/l)*ao*(c2-ntr));              % modal gain

g1=((vg*(d/l)*ao*(c1-ntr))*(1-(p1*eps11)-(eps12*p2))*(1-ks*(gamma_1*p1+gamma_2*p2)));             %linear modal gain
g2=((vg*(d/l)*ao*(c2-ntr))*(1-(p2*eps22)-(eps21*p1))*(1-ks*((p1*(1-gamma_1))+(p2*(1-gamma_2))))); %linear modal gain
if p1<0
    p11=-p1;
else
    p11=p1;
end
p11;
fs1=((-1+2*rand))*sqrt((2*beta*nsp*p11*0.1)/(0.1*taup1));


p1=(gamma_1*g1 + (1-gamma_1)*g2 - (1/taup1))*p11 + ((beta*nsp)/taup1)+fs1;
end