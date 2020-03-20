/*
Costa Junior, C. J., Garcia-Cintado, A. C., & Usabiaga, C. 
Fiscal Adjustments and the Shadow Economy in an Emerging Market. Macroeconomic Dynamics, 1-35.

Celso Costa
26/02/2018
*/

var Y Ym Yu C I G
K Km Ku N Nm Nu 
Wm Wu PIWm PIWu R RB P PI MC
S_P S_L Am Au pr S_G S_T S_tau_c S_tau_corp S_tau_s S_tau_l S_m
tau_c tau_l tau_s tau_corp B TE TAX T varphim MCm MCu
D_C D_G D_I D_Lm SELIC IPCA By D_ICM D_IRPJ D_IRPF D_Y D_Ym D_Yu D_L D_Lu Gr_Yu Gr_Nu D_TAX;

varexo e_P e_L e_Am e_Au e_pr e_G e_T e_tau_c e_tau_corp e_tau_s e_tau_l e_m e_varphim e_preco;

parameters alpha beta delta gammaB gammaMCm gammaMCu prss psim psiu s sigma tau_css tau_corpss tau_lss tau_sss theta thetaW varphi
rhoP rhoL rhoAm rhoAu rhopr rhoG rhoT rhotau_c rhotau_corp rhotau_s rhotau_l rhom rhovarphim
gammaR gammaPI gammaY gammaG gammaT gammatau_c gammatau_corp gammatau_l gammatau_s
phiG phiT phitau_c phitau_corp phitau_l phitau_s varphimss Omegam
Pss gammaCss gammaIss gammaGss gammaKss gammaWss;

//Calibrated Parameters ***************************************************
beta = 0.985;         //Cavalcanti and Vereda (2010)
delta = 0.025;          
alpha = 0.39;         //Kanczuc (2002)
prss = 0.03;
s = 1.3;

tau_css = 0.1594;     //Araújo and Ferreira (1999)
tau_lss = 0.1730;     //Araújo and Ferreira (1999)
tau_sss = 0.105;      //Cavalcanti and Silva (2010)
tau_corpss = 0.17;

varphimss = 7.758;    //Pablo Burriel, Jesús Fernández-Villaverde e Juan F. Rubio-Ramírez: MEDEA: a DSGE model for the Spanish economy 
varphi = 8.8;         //Castro et al (2010)
thetaW = 0.75;        //Castro et al (2010)
theta = 0.74;         //Castro et al (2010)
    
gammaR = 0.79;        //Castro et al (2010)
gammaPI = 2.43;       //Castro et al (2010)
gammaY = 0.16;        //Castro et al (2010)
gammaMCm = 0.4760;    //PNAD
gammaMCu = 1 - gammaMCm;

//Parameters estimated by the Bayesian method******************************
Omegam = 0.5;
psim = 1.5;
psiu = 1.7;
sigma = 2;

rhoP = 0.5;
rhoL = 0.5;
rhoAm = 0.5;
rhoAu = 0.5;
rhopr = 0.5;
rhoG = 0.5;
rhoT = 0.5;
rhotau_c = 0.5;
rhotau_corp = 0.5;
rhotau_s = 0.5;
rhotau_l = 0.5;
rhom = 0.5;
rhovarphim = 0.5;

gammaG = 0.1;
gammaT = 0; 
gammatau_c  = 0.1;
gammatau_corp  = 0;
gammatau_l  = 0;
gammatau_s = 0;
phiG  = -0.1;
phiT  = 0;
phitau_c  = 0.1;
phitau_corp  = 0;
phitau_l  = 0;
phitau_s = 0;

Pss = 1;

gammaCss = 0.61;
gammaIss = 0.20;
gammaGss = 0.19;
gammaKss = 2.6;
gammaB = 0.35;
gammaWss = 0.7;

model(linear);
//Steady state
#RBss = 1/beta;
#Rss = Pss*(1+tau_css)*((1/beta)-(1-delta));
#Yss = 1;
#Ymss = gammaMCm*Yss;
#Yuss = (1-gammaMCm)*Yss;
#Css = gammaCss*Yss;
#Iss = gammaIss*Yss;
#Gss = gammaGss*Yss;
#Kss = gammaKss*Yss;
#Kmss = gammaMCm*Kss;
#Kuss = (1-gammaMCm)*Kss;
#Bss = gammaB*Yss;
#Nmss = (Ymss*(Kmss^(-alpha)))^(1/(1-alpha));
#Nuss = (Yuss*(Kuss^(-alpha)))^(1/(1-alpha));
#Nss = Nmss + Nuss;
#Wmss = (varphimss/(varphimss-1))*(1/(1-beta*theta))*((Omegam*(Css^sigma)*(Nmss^psim)*Pss)/(1-tau_lss));
#Wuss = gammaWss*Wmss;
#TAXss = Pss*Gss - prss*s*tau_corpss*Yuss*Pss + Bss*(1-(1/RBss));
#Tss = TAXss - (tau_css*Pss*(Css+Iss) + tau_corpss*Ymss*Pss + (tau_lss+tau_sss)*Wmss*Nmss);
#MCmss = (Kmss/Ymss)*(Rss/alpha) + tau_corpss*Pss;
#MCuss = (Kuss/Yuss)*(Rss/alpha) + tau_corpss*Pss*prss*s;
#MCss = MCmss + MCuss;
#TEss = (tau_lss+tau_sss)*Wuss*Nuss + (1-prss)*tau_corpss*Yuss*Pss;
#A1 = 1/(MCuss-prss*s*tau_corpss*Pss);
#A2 = (MCuss - prss*s*tau_corpss*Pss)/(MCmss - tau_corpss*Pss);
//*************************************************************************
//Structural Model
//1-Capital movement law
K = (1-delta)*K(-1) + delta*I;
//2-Productivity shock
S_P = rhoP*S_P(-1) + e_P;
//3-Labor supply shock
S_L = rhoL*S_L(-1) + e_L;
//4-Supply of labor in the informal sector
S_L + sigma*C + psiu*Nu = Wu - P - (tau_css/(1+tau_css))*tau_c;
//5-Euler equation in relation to capital assets
(sigma/beta)*(C(+1)-C) = S_P(+1) - S_P + (Rss/(Pss*(1+tau_css)))*(R(+1) - P - (tau_css/(1+tau_css))*tau_c);
//6-Euler equation in relation to the government bond
sigma*(C(+1)-C) + PI(+1) + (tau_css/(1+tau_css))*(tau_c(+1)-tau_c) = S_P(+1) - S_P + RB;
//7-Phillips equation for wages (labor supply in the regular sector)
PIWm = beta*PIWm(+1) + ((1-thetaW)*(1-beta*thetaW)/thetaW)*((1/(1-varphimss))*varphim + S_L + sigma*C 
+ psim*Nm - (Wm-P) + (tau_lss/(1+tau_lss))*tau_l - (tau_css/(1+tau_css))*tau_c);
//8-Production function of the regular sector
Ym = Am + alpha*Km(-1) + (1-alpha)*Nm;
//9-Productivity shock in the regular sector
Am = rhoAm*Am(-1) + e_Am;
//10-Production function of the informal sector
Yu = Au + alpha*Ku(-1) + (1-alpha)*Nu;
//11-Productivity shock in the informal sector
Au = rhoAu*Au(-1) + e_Au;
//12-Shock in tax inspection
pr = rhopr*pr(-1) + e_pr;
//13-Tradeoff between labor market sector and the informal sector
Nm - Nu = A1*(A2*(MCmss*MCm - tau_corpss*Pss*(tau_corp + P)) - MCuss*MCu
+prss*s*tau_corpss*(pr + tau_corp + P)) + Ym - Yu + Wu - Wm - (tau_sss/(1+tau_sss))*tau_s;
//14-Tradeoff between capital of the market sector and the informal sector
Km(-1) - Ku(-1) = A1*(A2*(MCmss*MCm - tau_corpss*Pss*(tau_corp + P)) - MCuss*MCu
+prss*s*tau_corpss*(pr + tau_corp + P)) + Ym - Yu;
//15-Households budget constraint
Pss*(1+tau_css)*(Css*(C+P+(tau_css/(1+tau_css))*tau_c)+Iss*(I+P+(tau_css/(1+tau_css))*tau_c))
+ (Bss/RBss)*(B-RB(-1)) = (1-tau_lss)*Wmss*Nmss*(Wm + Nm - (tau_lss/(1-tau_lss))*tau_l) + Wuss*Nuss*(Wu + Nu)
+ Rss*Kss*(R + K(-1)) + Bss*B(-1) - Tss*T;
//16-Total marginal cost
MCss*MC = MCmss*MCm + MCuss*MCu;
//17-Phillips equation for prices
PI = beta*PI(+1) + ((1-theta)*(1-beta*theta)/theta)*(MC-P) + e_preco;
//18-Government budget constraint
(Bss/RBss)*(B-RB) - Bss*B(-1) = Pss*Gss*(P + G) - TAXss*TAX - prss*s*tau_corpss*Yuss*Pss*(pr + tau_corp + Yu + P);
//19-Total government revenue
TAXss*TAX = tau_css*Pss*(Css*(tau_c+P+C)+Iss*(tau_c+P+I))
+tau_corpss*Ymss*Pss*(tau_corp+Ym+P)+ Wmss*Nmss*(tau_lss*(tau_l+Wm+Nm)+tau_sss*(tau_s+Wm+Nm))+ Tss*T;
//20-Government spending movement law
P + G = gammaG*(G(-1) + P(-1)) + (1-gammaG)*phiG*(B(-1) - Ym(-1) - P(-1)) + S_G;
//21-Lumpsum tax movement law
T = gammaT*T(-1) + (1-gammaT)*phiT*(B(-1) - Ym(-1) - P(-1)) + S_T;
//22-Tax on consumption movement law
tau_c = gammatau_c*tau_c(-1) + (1-gammatau_c)*phitau_c*(B(-1) - Ym(-1) - P(-1)) + S_tau_c;
//23-Tax on companies movement law
tau_corp = gammatau_corp*tau_corp(-1) + (1-gammatau_corp)*phitau_corp*(B(-1) - Ym(-1) - P(-1)) + S_tau_corp;
//24-Labor movement tax law
tau_l = gammatau_l*tau_l(-1) + (1-gammatau_l)*phitau_l*(B(-1) - Ym(-1) - P(-1)) + S_tau_l;
//25-Labor contribution movement law
tau_s = gammatau_s*tau_s(-1) + (1-gammatau_s)*phitau_s*(B(-1) - Ym(-1) - P(-1)) + S_tau_s;
//26-Government spending shock
S_G = rhoG*S_G(-1) + e_G;
//27-Lumpsum tax shock
S_T = rhoT*S_T(-1) - e_T;
//28-Tax on consumption shock
S_tau_c = rhotau_c*S_tau_c(-1) - e_tau_c;
//29-Tax on companies shock
S_tau_corp = rhotau_corp*S_tau_corp(-1) - e_tau_corp;
//30-Labor movement shock
S_tau_l = rhotau_l*S_tau_l(-1) - e_tau_l;
//31-Labor contribution shock
S_tau_s = rhotau_s*S_tau_s(-1) - e_tau_s;
//32-Total tax evasion
TEss*TE = Wuss*Nuss*(tau_lss*(tau_l + Wu + Nu)+tau_sss*(tau_s + Wu + Nu))
+ tau_corpss*Yuss*Pss*(tau_corp + Yu + P -prss*(pr+tau_corp + Yu + P));
//33-Taylor rule
RB = gammaR*RB(-1) + (1-gammaR)*(gammaY*Ym + gammaPI*PI) + S_m;
//34-Monetary policy shock
S_m = rhom*S_m(-1) + e_m;
//35-Inflation definition
PI(+1) = P(+1) - P;
//36-Definition of wage inflation
PIWm = Wm - Wm(-1);
//37-Equilibrium condition in the goods market
Yss*Y = Css*C + Iss*I + Gss*G;
//38-Product aggregation (informal and regular sectors)
Yss*Y = Ymss*Ym + Yuss*Yu;
//39-Capital aggregation (informal and regular sectors)
Kss*K = Kmss*Km + Kuss*Ku;
//40-Labor aggregation (informal and regular sectors)
Nss*N = Nmss*Nm + Nuss*Nu;
//41-Defining the informal sector inflation rate
PIWu = Wu - Wu(-1);
//42-Shock in the elasticity of work in the regular sector
varphim = rhovarphim*varphim(-1) + e_varphim;
//43-Marginal cost of the regular sector
MCmss*MCm = ((((1+tau_sss)*Wmss)/(1-alpha))^(1-alpha))*((Rss/alpha)^alpha)
*((1-alpha)*(Wm+(tau_sss/(1+tau_sss))*tau_s) + alpha*R - Am) + tau_corpss*Pss*(tau_corp + P);
//44-Marginal cost of the informal sector
MCuss*MCu = (((Wuss/(1-alpha))^(1-alpha))*((Rss/alpha)^alpha))*((1-alpha)*Wu + alpha*R - Au) 
+ prss*s*tau_corpss*Pss*(pr + tau_corp + P);
//Endogenized Variables ************************************************* 
D_C = C + P - C(-1) - P(-1) + 0.024764992;
D_G = G + P - G(-1) - P(-1) + 0.026338249;
D_I = I + P - I(-1) - P(-1) + 0.021386032;
D_Lm = Nm - Nm(-1) - 0.00271715;
SELIC = RB + 1.032523462;
IPCA = PI - PI(-1) + 1.016183723;
By = B - Y + 0.275906154;
D_ICM = tau_css*Pss*(Css*(C+P+tau_c)+Iss*(I+P+tau_c)) 
- (tau_css*Pss*(Css*(C(-1)+P(-1)+tau_c(-1))+Iss*(I(-1)+P(-1)+tau_c(-1)))) + 0.02357596;
D_IRPF = (Wm+Nm+tau_l) - (Wm(-1)+Nm(-1)+tau_l(-1)) + 0.034165954;
D_IRPJ = tau_corpss*Ymss*Pss*((tau_corp+Ym+P)-(tau_corp(-1)+Ym(-1)+P(-1))) + 0.031640654;
D_Y = Y - Y(-1) + 0.025;
D_Ym = Ym - Ym(-1) - 0.00271715;
D_Yu = Yu - Yu(-1) - 0.016;
D_L = N - N(-1) + 0.025;
D_Lu = Nu - Nu(-1) - 0.016;
Gr_Yu = Yu - Y;
Gr_Nu = Nu - N;
D_TAX = TAX - TAX(-1);
end;

varobs D_C D_G D_I D_Lm SELIC IPCA By D_ICM D_IRPJ D_IRPF;

shocks;
var e_P; stderr 1;
var e_L ; stderr 1;
var e_Am; stderr 1;
var e_Au; stderr 1;
var e_pr; stderr 1;
var e_G; stderr 1;
var e_T; stderr 1;
var e_tau_c; stderr 1;
var e_tau_corp; stderr 1;
var e_tau_s; stderr 1;
var e_tau_l; stderr 1;
var e_m; stderr 1;
var e_varphim; stderr 1;
var e_preco; stderr 1;
end;

estimated_params;
//prss, uniform_pdf, , , 0.01, 0.03; //Renzo Orsi,Davi de Raggia, Francesco Turino: Size,trend,and policy implications of the undergroundeconomy
//s, uniform_pdf, , , 1, 1.3;        //Renzo Orsi,Davi de Raggia, Francesco Turino: Size,trend,and policy implications of the undergroundeconomy
//gammaMCm, beta_pdf, 0.4760, 0.1;  
//gammaB, uniform_pdf, , , 1, 2;  
//Omegam, uniform_pdf, , , 0.4, 0.5;
//Pss, uniform_pdf, , , 0.1, 2;
//sigma, gamma_pdf, 2, 0.1;
//varphi, gamma_pdf, 9, 0.1;
//thetaW, beta_pdf, 0.65, 0.1;
//theta, beta_pdf, 0.65, 0.05;
//gammaR, beta_pdf, 0.6, 0.15;       //Castro et al (2010)
//gammaPI, normal_pdf, 3, 0.1;       //Castro et al (2010)
//gammaY, beta_pdf, 0.25, 0.05;      //Castro et al (2010)
varphimss, gamma_pdf, 8, 0.5; 
psim, uniform_pdf, , , 1.4, 1.5;
psiu, uniform_pdf, , , 1.51, 1.8;
tau_css, beta_pdf, 0.1594, 0.01;     //Araújo and Ferreira (1999)
tau_lss, beta_pdf, 0.1730, 0.01;     //Araújo and Ferreira (1999)
tau_sss, beta_pdf, 0.105, 0.01;      //Cavalcanti and Silva (2010)
tau_corpss, uniform_pdf, , , 0.25, 0.35;
gammaG, uniform_pdf, , , 0.01, 0.8;
gammaT, uniform_pdf, , , 0.01, 0.99; 
gammatau_c, uniform_pdf, , , 0.5, 0.8;
gammatau_corp, uniform_pdf, , , 0.01, 0.99;
gammatau_l, uniform_pdf, , , 0.01, 0.99;
gammatau_s, uniform_pdf, , , 0.01, 0.99;
phiG, uniform_pdf, , , -0.8, -0.3;
phiT, uniform_pdf, , , 0, 1;
phitau_c, uniform_pdf, , , 0.3, 0.5;
phitau_corp, uniform_pdf, , , 0, 0.3;
phitau_l, uniform_pdf, , , 0, 1;
phitau_s, uniform_pdf, , , 0, 1;
//Componentes autoregressivos
rhoP, beta_pdf, 0.5, 0.25;
rhoL, beta_pdf, 0.5, 0.25;
rhoAm, beta_pdf, 0.5, 0.25;
rhoAu, beta_pdf, 0.5, 0.25;
rhopr, beta_pdf, 0.5, 0.25;
rhoG, beta_pdf, 0.5, 0.25;
rhoT, beta_pdf, 0.5, 0.25;
rhotau_c, beta_pdf, 0.5, 0.25;
rhotau_corp, beta_pdf, 0.5, 0.25;
rhotau_s, beta_pdf, 0.5, 0.25;
rhotau_l, beta_pdf, 0.5, 0.25;
rhom, beta_pdf, 0.5, 0.25;
rhovarphim, beta_pdf, 0.5, 0.25;
//Desvios padrão
stderr e_P, inv_gamma_pdf, 1, inf;
stderr e_L , inv_gamma_pdf, 1, inf;
stderr e_Am, inv_gamma_pdf, 1, inf;
stderr e_Au, inv_gamma_pdf, 1, inf;
stderr e_pr, inv_gamma_pdf, 1, inf;
stderr e_G, inv_gamma_pdf, 1, inf;
stderr e_T, inv_gamma_pdf, 1, inf;
stderr e_tau_c, inv_gamma_pdf, 1, inf;
stderr e_tau_corp, inv_gamma_pdf, 1, inf;
stderr e_tau_s, inv_gamma_pdf, 1, inf;
stderr e_tau_l, inv_gamma_pdf, 1, inf;
stderr e_m, inv_gamma_pdf, 1, inf;
stderr e_varphim, inv_gamma_pdf, 1, inf;
stderr e_preco, inv_gamma_pdf, 1, inf;
end;

identification;
//dynare_sensitivity;

// /*
estimation(datafile=BaseBrutaTrimestral_Final, mode_check, smoother, plot_priors = 0, mh_nblocks = 2, mh_drop=0.5, 
mh_jscale=0.3,mh_replic=1000000, mode_compute=6, bayesian_irf,irf=12);

shock_decomposition Gr_Yu Gr_Nu;

Resultados260218 = oo_;
save Resultados260218;
// */
