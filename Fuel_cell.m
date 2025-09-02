clear all
clc
W = zeros;
V = zeros;
F = zeros;
d = zeros;
g = zeros;
a = 1:100000;
tmem = 0.017;
za = 1.5;
zc = 3;
aa = 0.5;
ac = 1;
R = 8.314;
f = 96485;
P = 3;
N = 2;
XA = 0;
XC = 3.76;
imax = 2;
HHV = 286000;
i = input('Enter the step size required for I ');
j = input('Enter the step size required for T ');
k = input('Enter the step size required for A ');
l = input('Enter the step size required for No ');
I = 0.01:i:2;
T = 325:j:355;
A = 400:k:500;
No = 80:l:120;
s = length(I);
ss = length(T);
sss = length(A);
ssss = length(No);
for n=1:ss
    for m=1:s
        for o=1:sss
            for p=1:ssss
                io = 1.08*10^(-21)*exp(0.086*T(n));
                Psat = 10^(-2.1794+0.02953*T(n)-9.1837*10^(-5)*T(n)^2+1.4454*10^(-7)*T(n)^3)/100000;
                XHOA = Psat/P;
                XHOC = Psat/P;
                XH = (1-XHOA)/(1+(XA/2)*(1+(za/(za-1))));
                PH = P*XH;
                XO = (1-XHOC)/(1+(XC/2)*(1+(zc/(zc-1))));
                PO = P*XO;
                F(n,m) = -I(m)*((ac-aa)/(ac*aa))*((R*T(n))/(N*f))*log(I(m)/io);
                if ((PO/0.1173)+Psat)<2
                    B1 = (7.16*10^(-4)*T(n)-0.622)*((PO/0.1173)+Psat)+(-1.45*10^(-3)*T(n)+1.68);
                else
                    B1 = (8.66*10^(-5)*T(n)-0.068)*((PO/0.1173)+Psat);
                end
                B2 = 2.0;
                d(m) = -I(m)*I(m)*(B1*(I(m)/imax))^B2;
                a1 = XHOA*(P/Psat);
                lmem = 0.043+17.81*a1-39.85*a1^2+39.85*a1^3;
                g(n,m) = -I(m)*I(m)*tmem/((0.005139*lmem-0.00326)*exp(1268*((1/303)-(1/T(n)))));
                c(n,m) = (1.229-8.5*10^(-4)*(T(n)-298.15)+4.3085*10^(-5)*T(n)*(log(PH)+0.5*log(PO)));
                FH(m,o) = (No(p)*A(o)*I(m))/(2*f);
                W(n,m,o,p) = A(o)*No(p)*(I(m)*c(n,m)+d(m)+F(n,m)+g(n,m));
                V(n,m) = c(n,m)+d(m)/I(m)+F(n,m)/I(m)+g(n,m)/I(m);
                E(n,m,o,p) = W(n,m,o,p)/(FH(m,o)*HHV);
                Q(n,m,o,p) = (1-(E(n,m,o,p)/100))*W(n,m,o,p);
                if W(n,m,o,p)>5000 && W(n,m,o,p)<5005
                        Power(a) = W(n,m,o,p);
                        Temp(a) = T(n);
                        Area(a) = A(o);
                        Current(a) = I(m);
                        Flow_Rate(a) = FH(m);
                        Flow(a) = Flow_Rate(a)*22.4*3600;
                        Voltage(a) = V(n,m);
                        Number(a) = No(p);
                        Efficiency(a) = E(n,m,o,p);
                        Heat_Loss(a) = Q(n,m,o,p);
                end
                a=a+1;
            end
        end
    end
end
plot(I,V,'*');
xlabel('Current (A/cm^2)');
ylabel('Voltage (V)');