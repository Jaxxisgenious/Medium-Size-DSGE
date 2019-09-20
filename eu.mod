%%% This code was greated by Jaxx Yu Kang based on the medium sized DSGE
%%% model published by Prof. Wenli Xu
var lamb C i pi R Z I v pssi w w_star h1 h2 N_d K_hat K 
    mc pi_star x1 x2 Y G A v_p 
    tauc taun tauk
    B M
    i_m q BC i_bc u S betta delta1 cm cbc cbond pint pi_m pi_bc C_Y I_Y B_Y BC_Y g;

varexo e_i e_z e_v e_pssi
       e_tauc e_taun e_tauk e_g e_m e_a e_cm;
       
parameters cb alfa delta0 delta2 pi_target ki fiw epsilonw ka
           etaw epsilonp etap fip omegga fipi fiy rhoa rhoz rhog rhov rhopsi
           psis rhotauc taucs rhotaun tauns rhotauk tauks  fin fik 
           cgamma_m chi_m c2 rhom rho1 rho2 rho3 rhocm cm_ss spread;
cb=0.86;
alfa=0.25;
pi_target=0.00025;
delta0 = 0.045;
%delta0 = 0.08;
%delta0 = 0.12;
delta2=0.05;
ki=2.5;
epsilonw=1.8;
%epsilonp=200;
epsilonp=1000;
ka=2;
etaw=0.2;
etap=0.2;
fiw=0.5;
fip=0.5;
fipi=2;
fiy=0.5;
rhoa=0.92;
rhoz=0.94;
rhog=0.72;
rhov=0.9;
rhopsi=0.7;
psis=1;
%psis = 0.1;
rhotauc=0.9;
taucs=0.2;
rhotaun=0.9;
tauns=0.44;
%tauns=0;
rhotauk=0.9;
tauks=0.3;
fin=0.5;
fik=0.5;
omegga = .23;
cgamma_m = 0.01;
chi_m = 2;
c2 = 0.74;
rhom = 0.8;
%rho1 = 0.7;
%rho2 = 0.05;
%rho3 = 0.005;
rho1 = 0.2;
rho2 = 0.1;
rho3 = 0.0018;
rhocm = 0.9;
cm_ss = 0;
spread = 0.0125;
model;

[name='FOC for Consumption']
(1+tauc)*lamb=v/(C-cb*C(-1))-betta*cb*v(1)/(C(1)-cb*C);

[name='FOC for Bond Holding']
lamb=betta*lamb(1)*(1+i)*(1+pi(1))^(-1);

[name='FOC for CBDC']
i_m = 1/betta - 1/lamb*(v*cgamma_m*(M^chi_m))-1;

[name='FOC for Wage']
N_d=((epsilonw-1)/epsilonw*w_star/(psis*v)*(1-taun)/(1+tauc)/(1-cb)*(1-betta*cb)*((1-omegga)*(alfa*mc/R)^(alfa/(1-alfa))-delta0*(alfa*mc/R)^(1/(1-alfa)))^(-1)*(w/w_star)^(-epsilonw*ka)*((1-fiw*betta*(1+pi)^(epsilonw*(1+ka)*(1-etaw)))/(1-fiw*betta*(1+pi)^((epsilonw-1)*(1-etaw)))))^(1/(1+ka));

[name='Def h1']
h1=v*pssi*(w/w_star)^(epsilonw*(1+ka))*N_d^(1+ka)+fiw*betta*(1+pi)^(-etaw*epsilonw*(1+ka))*(1+pi(1))^(epsilonw*(1+ka))*(w_star(1)/w_star)^(epsilonw*(1+ka))*h1(1);

[name='Def h2']
h2=(1-taun)*lamb*(w/w_star)^epsilonw*N_d+fiw*betta*(1+pi)^(etaw*(1-epsilonw))*(1+pi(1))^(epsilonw-1)*(w_star(1)/w_star)^epsilonw*h2(1);

[name='Def betta']
betta = (1+pi_target)/(1+(rho3/(1-rho2))*(1-cm));

[name='FOC for K & N']
w/R=(1-alfa)/alfa*(K_hat(-1)/N_d);

[name='FOC for Price']
1+pi_star=epsilonp/(epsilonp-1)*(1+pi)*x1/x2;

[name='Real Wage']
w^(1-epsilonw)=(1-fiw)*w_star^(1-epsilonw)+(1+pi(-1))^(etaw*(1-epsilonw))*fiw*(1+pi)^(epsilonw-1)*w(-1)^(1-epsilonw);

[name='Def Price Diffusion Index']
v_p=(1-fip)*((1+pi_star)/(1+pi))^(-epsilonp)+(1+pi(-1))^(-etap*epsilonp)*(1+pi)^epsilonp*fip*v_p(-1);

[name='Def x1']
x1=lamb*mc*Y+fip*betta*(1+pi)^(-etap*epsilonp)*(1+pi(1))^epsilonp*x1(1);

[name='Def x2']
x2=lamb*Y+fip*betta*(1+pi)^(etap*(1-epsilonp))*(1+pi(1))^(epsilonp-1)*x2(1);

[name='Def mc']
mc=w/((1-alfa)*A*(K_hat(-1)/N_d)^alfa);

[name = 'FOC for capital utility']
q * (delta1+delta2*(u-1)) = (1-tauk) * R;

[name = 'FOC for capital accumulation']
R = q/(1-tauk )*((1+i-pi)/betta+delta0-1);

[name = 'Def BC']
i_bc = i + spread;

[name='Def K_hat']
K_hat=K*u;

[name='Def Capital Accumulation']
%K=Z*(1-ki/2*(I/I(-1)-1)^2)*I+(1-delta0)*K(-1);
K=Z*(1-ki/2*(I/I(-1)-1)^2)*I+(1-delta0-delta1*(u-1)-delta2/2*(u-1)^2)*K(-1);
%K(1)=Z(1)*(1-ki/2*(I(1)/I-1)^2)*I(1)+(1-delta0-delta1*(u(1)-1)-delta2/2*(u(1)-1)^2)*K;

[name='Def capital utility']
q = 1/((1+pi(1))*(1+i-pi));

[name='Def delta1']
delta1 = ((1+pi)/betta-pi)/betta + delta0 -1;

[name='Def inflation rate']
(1+pi)^(1-epsilonp)=(1-fip)*(1+pi_star)^(1-epsilonp)+(1+pi(-1))^(etap*(1-epsilonp))*fip;

[name='Def Taylor rule of Monetary Policy']
i=(rho2-rho1)*steady_state(i)+rho1*i(-1)+rho3*(1-M(-1)/Y(-1))+(rho2-rho1)*(fipi*(pi-pi_target)+fiy*(log(Y)-log(Y(-1))))+e_i;

[name='CBDC holding constraint']
M =(cm* Y(-1))*S;

[name='Def cm']
cm = (1-rhocm)*cm_ss+rhocm*cm(-1)+ e_cm;

[name='Government expenditure']
log(G)=(1-rhog)*log(omegga*Y)+rhog*log(G(-1))+e_g;

[name = 'Bond Holding constraint']
(B+M)/Y = c2;

[name='Def productivity shock']
log(A)=rhoa*log(A(-1))+e_a;

[name='Def marginal efficiency of investment shock']
log(Z)=rhoz*log(Z(-1))+e_z;

[name='Def intratemporal preference shock']
log(v)=rhov*log(v(-1))+e_v;

[name='Def intratemporal preference (labor supply) shock']
log(pssi)=(1-rhopsi)*log(psis)+rhopsi*log(pssi(-1))+e_pssi;

[name='Def tauc']
tauc=(1-rhotauc)*taucs+rhotauc*tauc(-1)+e_tauc;

[name='Def taun']
taun=(1-rhotaun)*(tauns+fin*log(w*N_d/(steady_state(w)*steady_state(N_d))))+rhotaun*taun(-1)+e_taun;

[name='Def tauk']
tauk=(1-rhotauk)*(tauks+fik*log(R*K_hat(-1)/(steady_state(R)*steady_state(K_hat))))+rhotauk*tauk(-1)+e_tauk;

[name='Goods Market Clear']
Y=C+I+G;

[name='Goods Market Clear II']
Y*v_p=A*K_hat(-1)^alfa*N_d^(1-alfa);

[name='Credit Market Clear']
BC=I;

[name='Def CBDC supply shock']
log(S)=rhom*log(S(-1))+e_m;

[name='Def cbc']
cbc = BC/Y;

[name='Def cbond']
cbond = B/Y;

[name='Def interest rate in percentage term']
pint = i*100;

[name='Def CBDC interest rate in percentage term']
pi_m = i_m*100;

[name='Def BC interest rate in percentage term']
pi_bc = i_bc*100;

[name='Def C_Y']
C_Y = C/Y;

[name='Def I_Y']
I_Y = I/Y;

[name='Def B_Y']
B_Y = B/Y;

[name='Def BC_Y']
BC_Y = BC/Y;
   
[name='Def g']
g = (Y-Y(-1))/Y(-1)*100;
end;

initval;
A=1;
Z=1;
v=1;
u=1;
S=1;
pssi=psis;
tauc=taucs;
taun=tauns;
tauk=tauks;
pi=pi_target;
cm = cm_ss;
betta = (1+pi_target)/(1+(rho3/(1-rho2))*(1-cm));
delta1 = ((1+pi_target)/betta-pi_target)/betta + delta0 -1;
i = rho3/(1-rho2)*(1-cm);
q = 1/((1+pi)*(1+i-pi));
i_bc = i + 0.0125;
R = q*delta1/(1-tauk);
pi_star=(((1+pi)^(1-epsilonp)-fip*(1+pi)^(etap*(1-epsilonp)))/(1-fip))^(1/(1-epsilonp))-1;
v_p=(1-fip)*((1+pi_star)/(1+pi))^(-epsilonp)/(1-fip*(1+pi)^(epsilonp*(1-etap)));
mc=(epsilonp-1)/epsilonp*(1+pi_star)/(1+pi)*(1-fip*betta*(1+pi)^(epsilonp*(1-etap)))/(1-fip*betta*(1+pi)^((1-epsilonp)*(etap-1)));
w=(1-alfa)*mc*(alfa*mc/R)^(alfa/(1-alfa));
w_star=w*((1-fiw*(1+pi)^((epsilonw-1)*(1-etaw)))/(1-fiw))^(1/(1-epsilonw));
N_d=((epsilonw-1)/epsilonw*w_star/(psis*v)*(1-taun)/(1+tauc)/(1-cb)*(1-betta*cb)*((1-omegga)*(alfa*mc/R)^(alfa/(1-alfa))-delta0*(alfa*mc/R)^(1/(1-alfa)))^(-1)*(w/w_star)^(-epsilonw*ka)*((1-fiw*betta*(1+pi)^(epsilonw*(1+ka)*(1-etaw)))/(1-fiw*betta*(1+pi)^((epsilonw-1)*(1-etaw)))))^(1/(1+ka));
K_hat=(alfa*mc/R)^(1/(1-alfa))*N_d;
K=K_hat;
I=delta0*K;
Y=A*K_hat^alfa*N_d^(1-alfa)/v_p;
M = cm * Y;
B = c2*Y - M;
G= omegga*Y;
BC=I;
C = Y - I - G;
lamb = 1/C*(1-betta*cb)/(1-cb)/(1+tauc);
%i_m = 1/lamb - (v*cgamma_m*(M^chi_m)/lamb)-1;
i_m = 1/betta - (v*cgamma_m*(M^chi_m)/lamb)-1;
x1=lamb*mc*Y/(1-fip*betta*(1+pi)^((1-etap)*epsilonp));
x2=lamb*Y/(1-fip*betta*(1+pi)^((etap-1)*(1-epsilonp)));
h1=psis*(w/w_star)^(epsilonw*(1+ka))*N_d^(1+ka)/(1-fiw*betta*(1+pi)^(epsilonw*(1+ka)*(1-etaw)));
h2=(1-taun)*lamb*(w/w_star)^epsilonw*N_d/(1-fiw*betta*(1+pi)^((epsilonw-1)*(1-etaw)));
pint = i*100;
pi_m = i_m*100;
pi_bc = i_bc*100;
C_Y=C/Y;
I_Y=I/Y;
B_Y=B/Y;
BC_Y=BC/Y;
end;


steady;
check;



model_diagnostics;
shocks;
%var e_a;
%periods 1:40 40:50 50:500;
%values 0 0.0001 0;
var e_cm;
periods 1:50 50:500;
values 0 0.015;
end;

simul(periods = 500);
%stoch_simul(order=2);




