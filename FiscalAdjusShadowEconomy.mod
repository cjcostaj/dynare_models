/*
Costa Junior, C. J., Garcia-Cintado, A. C., & Usabiaga, C. 
Fiscal Adjustments and the Shadow Economy in an Emerging Market. Macroeconomic Dynamics, 1-35.

Celso Costa
10/09/2018
*/

var Y C I G IGm IGu YF
K Km Ku N Nm Nu KGm KGu INPD INPm INPu INPDD INPDF INPFD
Wm Wu PIWm PIWu R RB P PI MC MCm MCu RN RF BF Poup LAMBDA S PD PF Rf
S_P S_L Am Au pr S_G S_T S_tau_c S_tau_corp S_tau_s S_tau_l S_m S_IGm S_IGu
tau_c tau_l tau_s tau_corp B TE TAX T varphim
DC DG DI DX DIMP DRB DIPCA DRF DIRPF DIRPJ DTAX DPREV DYF DPF DNm DS DICM DB DBF DTE Gr_Nu;

varexo e_P e_L e_Am e_Au e_pr e_G e_IGm e_IGu e_T e_tau_c e_tau_corp e_tau_s e_tau_l e_m e_varphim e_PF e_RF e_YF e_preco e_gov;

parameters alpha1 alpha2 alpha3 beta delta gammaB gammaMCm gammaMCu prss psim psiu s sigma 
tau_css tau_corpss tau_lss tau_sss theta thetaW varphi
rhoP rhoL rhoAm rhoAu rhopr rhoG rhoT rhotau_c rhotau_corp rhotau_s rhotau_l rhom rhovarphim rhoIGm rhoIGu
gammaR gammaPI gammaY gammaG gammaT gammatau_c gammatau_corp gammatau_l gammatau_s
phiG phiT phitau_c phitau_corp phitau_l phitau_s varphimss Omegam
Pss gammaCss gammaIss gammaGss gammaWss gammaImss gammaIuss gammaBFss
gammaIGm gammaIGu phiIGm phiIGu
omegaD omegaF psiD psiF phic chiBF rhoPF rhoRF rhoYF
PDss PFss omegaRf;

//Parâmetros Calibrados ***************************************************
sigma = 2;                                                                 //Cavalcanti and Silva (2010)
beta = 0.989;                                                              //Castro et al (2010)
delta = 0.015;                                                             //Castro et al (2010)
prss = 0.03;                                                               //Renzo Orsi,Davi de Raggia, Francesco Turino: Size,trend,and policy implications of the undergroundeconomy
s = 1.3;                                                                   //Renzo Orsi,Davi de Raggia, Francesco Turino: Size,trend,and policy implications of the undergroundeconomy
alpha2 = 0.6;                                                              //Mussolini (2011)
alpha1 = 0.3;                                                              //Mussolini (2011)                                                              
alpha3 = 0.4-alpha1;

varphimss = 7.758;                                                         //Pablo Burriel, Jesús Fernández-Villaverde e Juan F. Rubio-Ramírez: MEDEA: a DSGE model for the Spanish economy 
varphi = 11;                                                               //Castro et al (2010)
thetaW = 0.75;                                                             //Castro et al (2010)
theta = 0.74;                                                              //Castro et al (2010)
    
gammaR = 0.79;                                                             //Castro et al (2010)
gammaPI = 2.43;                                                            //Castro et al (2010)
gammaY = 0.16;                                                             //Castro et al (2010)
gammaMCm = 0.4760;                                                         //PNAD
gammaMCu = 1 - gammaMCm;

Pss = 1;                                                                   //NUMERARIO

gammaCss = 0.61;                                                           //CONTAS NACIONAIS
gammaIss = 0.17;                                                           //CONTAS NACIONAIS
gammaGss = 0.20;                                                           //CONTAS NACIONAIS
gammaB = 0.35;
gammaWss = 0.7;
gammaIuss = 0.015;                                                         //Mussolini e Teles
gammaImss = 0.015;                                                         //Mussolini e Teles
gammaBFss = -0.1;

omegaD = 0.15;                                                             //Participação dos insumos importados na produção do bem doméstico
omegaF = 0.02;                                                             //Participação dos insumos produzidos na economia doméstica na produção do bem resto do mundo
psiD = 2;
psiF = 3;

phic = 0.74;                                                               //Castro et al (2010)

chiBF = -0.003;
PDss = 4;
PFss = 6;

//Parâmetros estimados por método bayesiano
Omegam = 0.5;
psim = 1.5;
psiu = 1.7;

tau_css = 0.1594;     //Araújo e Ferreira (1999)
tau_lss = 0.1730;     //Araújo e Ferreira (1999)
tau_sss = 0.105;      //Cavalcanti and Silva (2010)
tau_corpss = 0.17;

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
rhoPF = 0;
rhoRF = 0;
rhoYF = 0;

rhoIGm = 0.5;
rhoIGu = 0.5;

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

gammaIGm = 0; 
gammaIGu = 0; 
phiIGm = 0; 
phiIGu = 0;

omegaRf = 0.5;

model(linear);
//Estado estacionário
#RNss = 1/beta;
#Rfss = (1-omegaRf) + omegaRf*RNss;
#RBss = 1/beta;
#RFss = 1/beta;
#Rss = Pss*(1+tau_css)*((1/beta)-(1-delta));
#Yss = 1.67;                                                               //C+I+G no segundo trimestre de 2018 (1.667.756,1964 R$ milhões)
#Css = gammaCss*Yss;
#Iss = gammaIss*Yss;
#IGmss = gammaImss*Yss;
#IGuss = gammaIuss*Yss;
#Gss = gammaGss*Yss;
#Kss = Iss/delta;
#Kmss = gammaMCm*Kss;
#Kuss = (1-gammaMCm)*Kss;
#KGmss = 0.2*Kmss;                                                         //Mussolini e Teles
#KGuss = 0.2*Kuss;                                                         //Mussolini e Teles
#Bss = gammaB*Yss;
#BFss = gammaBFss*Yss;
#INPDss = (1-omegaD)*Yss;
#INPDDss = 0.85*INPDss;
#INPDFss = 0.15*INPDss;
#YFss = (1/omegaF)*Yss;
#INPFDss = (1-omegaF)*Yss;
#INPmss = gammaMCm*INPDss;
#INPuss = gammaMCu*INPDss;
#Sss = (PDss/PFss)*((omegaD/(1-omegaD))*(INPDDss/INPFDss))^(1/psiD);
#Nmss = (INPmss*(Kmss^(-alpha1)))^(1/alpha2);
#Nuss = (INPuss*(Kuss^(-alpha1)))^(1/alpha2);
#Nss = Nmss + Nuss;
#Wmss = (varphimss/(varphimss-1))*(1/(1-beta*theta))*((Omegam*(Css^sigma)*(Nmss^psim)*Pss)/(1-tau_lss));
#Wuss = gammaWss*Wmss;
#TAXss = Pss*Gss - prss*s*tau_corpss*INPuss*Pss + Bss*(1-(1/RBss));
#Tss = TAXss - (tau_css*Pss*(Css+Iss) + tau_corpss*INPmss*Pss + (tau_lss+tau_sss)*Wmss*Nmss);
#MCmss = (Kmss/INPmss)*(Rss/alpha1) + tau_corpss*Pss;
#MCuss = (Kuss/INPuss)*(Rss/alpha1) + tau_corpss*Pss*prss*s;
#MCss = MCmss + MCuss;
#TEss = (tau_lss+tau_sss)*Wuss*Nuss + (1-prss)*tau_corpss*INPuss*Pss;
#A1 = 1/(MCuss-prss*s*tau_corpss*Pss);
#A2 = (MCuss - prss*s*tau_corpss*Pss)/(MCmss - tau_corpss*Pss);
#Poupss = 0.4*Yss;
//*************************************************************************
//Modelo Estrutural
//FAMÍLIAS XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//1-Restrição orçamentária da família
Pss*(1+tau_css)*(Css*(C+P+(tau_css/(1+tau_css))*tau_c)+Iss*(I+P+(tau_css/(1+tau_css))*tau_c))
+ (Bss/RBss)*(B-RB(-1)) + RFss*BFss*Sss*(RF(-1) + BF(-1) + S) + Poupss*Poup 
= (1-tau_lss)*Wmss*Nmss*(Wm + Nm - (tau_lss/(1-tau_lss))*tau_l) + Wuss*Nuss*(Wu + Nu)
+ Rss*Kss*(R + K(-1)) + Bss*B(-1) + BFss*BF + Poupss*RNss*(RN(-1) + Poup(-1)) - Tss*T;
//2-Lei de movimento do capital
K = (1-delta)*K(-1) + delta*I;
//3-Choque de preferência
S_P = rhoP*S_P(-1) + e_P;
//4-Choque de oferta de trabalho
S_L = rhoL*S_L(-1) + e_L;
//5-Escolha do consumo
(1-phic*beta)*(LAMBDA + P + (tau_css/(1+tau_css))*tau_c) = S_P - (sigma/(1-phic))*(C - phic*C(-1)) 
- (beta*phic/(1-phic))*(S_P(+1) - sigma*(C(+1) - phic*C));
//6-Escolha do trabalho no setor informal
LAMBDA = S_P + S_L + psiu*Nu - Wu;
//7-Equação de Euler do bem de capital
(Pss*(1+tau_css)/beta)*(LAMBDA-LAMBDA(+1) + P + (tau_css/(1+tau_css))*tau_c) 
= (1-delta)*Pss*(1+tau_css)*(P(+1) + (tau_css/(1+tau_css))*tau_c) + Rss*R(+1); 
//8-Equação de Euler (Título Público)
LAMBDA - RB = LAMBDA(+1);
//9-Equação de Euler (Título Extrangeiro)
RF + LAMBDA(+1) + S(+1) = LAMBDA + S - chiBF*BFss*BF;
//10-Equação de Euler para o título emitidos pelas firmas
LAMBDA - RN = LAMBDA(+1);
//11-Equação de Phillips para os salários (oferta de trabalho no setor regular)
PIWm = beta*PIWm(+1) + ((1-thetaW)*(1-beta*thetaW)/thetaW)*((1/(1-varphimss))*varphim 
+ S_P + S_L + psim*Nm - LAMBDA + (tau_lss/(1+tau_lss))*tau_l - Wm);                                     
//12-Choque na elasticidade do trabalho no setor regular
varphim = rhovarphim*varphim(-1) + e_varphim;
//FIRMAS XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//13-Função de produção do setor regular
INPm = Am + alpha1*Km(-1) + alpha2*Nm + alpha3*KGm(-1);
//14-Choque de produtividade no setor regular
Am = rhoAm*Am(-1) + e_Am;
//15-Função de produção do setor informal
INPu = Au + alpha1*Ku(-1) + alpha2*Nu + alpha3*KGu(-1);
//16-Choque de produtividade no setor informal
Au = rhoAu*Au(-1) + e_Au;
//17-Choque na fiscalização do fisco
pr = rhopr*pr(-1) + e_pr;
//18-Tradeoff entre trabalho do setor mercado e do setor informal
Nm - Nu = A1*(A2*(MCmss*MCm - tau_corpss*Pss*(tau_corp + P)) - MCuss*MCu
+prss*s*tau_corpss*(pr + tau_corp + P)) + INPm - INPu + Wu - (Wm + Rf) - (tau_sss/(1+tau_sss))*tau_s;
//19-Tradeoff entre capital do setor mercado e do setor informal
Km(-1) - Ku(-1) = A1*(A2*(MCmss*MCm - tau_corpss*Pss*(tau_corp + P)) - MCuss*MCu
+prss*s*tau_corpss*(pr + tau_corp + P)) + INPm - INPu;
//20-Custo do empréstimo da firma
Rfss*Rf = omegaRf*RNss*RN;
//21-Custo marginal total (Preço do insumo doméstico)
PDss*PD = MCmss*MCm + MCuss*MCu;
//22-Custo marginal do setor mercado
MCmss*MCm = ((((1+tau_sss)*Wmss*RNss)/alpha2)^alpha2)*((Rss/alpha1)^alpha1)*(1/KGmss^alpha3)
*(alpha2*(Wm + Rf +(tau_sss/(1+tau_sss))*tau_s) + alpha1*R - Am - alpha3*KGm(-1)) + tau_corpss*Pss*(tau_corp + P);
//23-Custo marginal do setor informal
MCuss*MCu = (((Wuss/alpha2)^alpha2)*((Rss/alpha1)^alpha1))*(1/KGuss^alpha3)*(alpha2*Wu + alpha1*R - Au - alpha3*KGu(-1)) 
+ prss*s*tau_corpss*Pss*(pr + tau_corp + P);
//24-Insumos produzidos domesticamente entre os setores
INPDss*INPD = INPmss*INPm + INPuss*INPu;
//25-Insumos produzidos domesticamente
INPDss*INPD = INPDDss*INPDD + INPDFss*INPDF;
//26-Função de Produção do Produto Intermediário
(Yss^((psiD-1)/psiD))*Y = ((1-omegaD)^(1/psiD))*(INPDDss^((psiD-1)/psiD))*INPDD 
+ (omegaD^(1/psiD))*(INPFDss^((psiD-1)/psiD))*INPFD;
//27-Demanda por insumos domésticos
INPDD = psiD*(MC - PD) + Y;
//28-Demanda por insumos importados
INPFD = psiD*(MC -S - PF) + Y;
//29-Custo Marginal da Firma Produtora de Bens Intermediários
(MCss^(1-psiD))*MC = omegaD*(PDss^(1-psiD))*PD + (1-omegaD)*((Sss*PFss)^(1-psiD))*(S+PF);
//30-Equação de Phillips para os preços
PI = beta*PI(+1) + ((1-theta)*(1-beta*theta)/theta)*(MC-P) + e_preco;
//GOVERNO XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//31-Restrição orçamentária do governo
(Bss/RBss)*(B-RB) - Bss*B(-1) = Pss*Gss*(P + G) + Pss*IGmss*(P + IGm) + Pss*IGuss*(P + IGu) 
- TAXss*TAX - prss*s*tau_corpss*INPuss*Pss*(pr + tau_corp + INPu + P)
- Wmss*Nmss*tau_sss*(tau_s+Wm+Nm) + e_gov;
//32-Receita total do governo
TAXss*TAX = tau_css*Pss*(Css*(tau_c+P+C)+Iss*(tau_c+P+I))
+tau_corpss*INPmss*Pss*(tau_corp+INPm+P)+ Wmss*Nmss*tau_lss*(tau_l+Wm+Nm)+ Tss*T;
//33-Lei de movimento do gasto do governo
P + G = gammaG*(G(-1) + P(-1)) + (1-gammaG)*phiG*(B(-1) - Y(-1) - P(-1)) + S_G;
//34-Lei de movimento do investimento do governo no setor formal
P + IGm = gammaIGm*(IGm(-1) + P(-1)) + (1-gammaIGm)*phiIGm*(B(-1) - Y(-1) - P(-1)) + S_IGm;
//35-Lei de movimento do investimento do governo no setor informal
P + IGu = gammaIGu*(IGu(-1) + P(-1)) + (1-gammaIGu)*phiIGu*(B(-1) - Y(-1) - P(-1)) + S_IGu;
//36-Lei de movimento do tributo lumpsum
T = gammaT*T(-1) + (1-gammaT)*phiT*(B(-1) - Y(-1) - P(-1)) + S_T;
//37-Lei de movimento do tributo sobre consumo
tau_c = gammatau_c*tau_c(-1) + (1-gammatau_c)*phitau_c*(B(-1) - Y(-1) - P(-1)) + S_tau_c;
//38-Lei de movimento do tributo sobre as empresas
tau_corp = gammatau_corp*tau_corp(-1) + (1-gammatau_corp)*phitau_corp*(B(-1) - Y(-1) - P(-1)) + S_tau_corp;
//39-Lei de movimento do tributo sobre trabalho
tau_l = gammatau_l*tau_l(-1) + (1-gammatau_l)*phitau_l*(B(-1) - Y(-1) - P(-1)) + S_tau_l;
//40-Lei de movimento da contribuição trabalhista
tau_s = gammatau_s*tau_s(-1) + (1-gammatau_s)*phitau_s*(B(-1) - Y(-1) - P(-1)) + S_tau_s;
//41-Choque no gasto do governo
S_G = rhoG*S_G(-1) - e_G;
//42-Choque no investimento do governo no setor formal
S_IGm = rhoIGm*S_IGm(-1) - e_IGm;
//43-Choque no investimento do governo no setor formal
S_IGu = rhoIGu*S_IGu(-1) - e_IGu;
//44-Choque no tributo lumpsum
S_T = rhoT*S_T(-1) + e_T;
//45-Choque no tributo sobre o consumo
S_tau_c = rhotau_c*S_tau_c(-1) + e_tau_c;
//46-Choque no tributo sobre as empresas
S_tau_corp = rhotau_corp*S_tau_corp(-1) + e_tau_corp;
//47-Choque no tributo sobre o trabalho
S_tau_l = rhotau_l*S_tau_l(-1) + e_tau_l;
//48-Choque no contribuição trabalhista
S_tau_s = rhotau_s*S_tau_s(-1) + e_tau_s;
//49-Lei de movimento do capital público do setor formal
KGm = (1-delta)*KGm(-1) + delta*IGm;
//50-Lei de movimento do capital público do setor formal
KGu = (1-delta)*KGu(-1) + delta*IGu;
//51-Evasão fiscal total
TEss*TE = Wuss*Nuss*(tau_lss*(tau_l + Wu + Nu)+tau_sss*(tau_s + Wu + Nu))
+ tau_corpss*INPuss*Pss*(tau_corp + INPu + P - prss*(pr+tau_corp + INPu + P));
//52-Regra de Taylor
RB = gammaR*RB(-1) + (1-gammaR)*(gammaY*Y + gammaPI*PI) + S_m;
//53-Choque de política monetária
S_m = rhom*S_m(-1) + e_m;
//54-Definição de inflação
PI(+1) = P(+1) - P;
//55-Definição de inflação de salários
PIWm = Wm - Wm(-1);
//56-Definição da taxa de inflação do setor informal
PIWu = Wu - Wu(-1);
//SETOR EXTERNO XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//57-Demanda pelo insumo doméstico no resto do mundo
INPDF = psiF*(S + PF - PD) + YF;
//58-Balanço de Pagamentos
Sss*BFss*(S+BF)-Sss*RFss*BFss*(S+RF(-1)+BF(-1)) = PFss*Sss*INPFDss*(PF + S + INPFD) - PDss*INPDFss*(PD + INPDF);
//59-Choque no produto do resto do mundo
YF = rhoYF*YF(-1) + e_YF;
//60-Choque nas Taxas de Juros Internacional
RF = rhoRF*RF(-1) + e_RF;
//61-Choque no Nível de Preços Internacional
PF = rhoPF*PF(-1) + e_PF;
//CONDIÇÃO DE EQUILÍBRIO XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//62-Condição de equilíbrio no mercado de bens
Yss*Y = Css*C + Iss*I + Gss*G  + IGmss*IGm + IGuss*IGu + e_YF;
//63-Agregação do capital (setores informal e regular)
Kss*K = Kmss*Km + Kuss*Ku;
//64-Agregação do trabalho (setores informal e regular)
Nss*N = Nmss*Nm + Nuss*Nu;
//Variáveis Endogenizadas (65-85)****************************************** 
DC = C - C(-1) + 0.030762799;
DG = G - G(-1) + 0.019556704;
DI = I - I(-1) + 0.026821303;
DX = INPDF - INPDF(-1) + 0.042695399;
DIMP = INPFD - INPFD(-1) + 0.067699417;
DRB = RB - RB(-1) - 0.000688522;
DIPCA = PI - PI(-1) - 0.000122768;
DRF = RF - RF(-1) - 0.000035304;
DIRPF = (Wm+Nm+tau_l) - (Wm(-1)+Nm(-1)+tau_l(-1)) + 0.030261618;
DIRPJ = tau_corpss*INPmss*Pss*((tau_corp+INPm+P)-(tau_corp(-1)+INPm(-1)+P(-1))) + 0.020697347;
DTAX = TAX - TAX(-1) + 0.028281255;
DPREV = Wmss*Nmss*tau_sss*((tau_s+Wm+Nm) - ((tau_s(-1)+Wm(-1)+Nm(-1)))) + 0.027611405;
DYF = YF - YF(-1) + 0.004850033;
DPF = PF - PF(-1) + 0.005220996;
DNm = Nm - Nm(-1) - 0.001073551;
DS = S - S(-1) + 0.001562958;
DICM = tau_css*Pss*(Css*(C+P+tau_c)+Iss*(I+P+tau_c)) 
- (tau_css*Pss*(Css*(C(-1)+P(-1)+tau_c(-1))+Iss*(I(-1)+P(-1)+tau_c(-1)))) + 0.023084938;
DB = B - Y - 0.488974011;
DBF = BF - Y - 0.076629944;
DTE = TE - TE(-1);
Gr_Nu = Nu - N - 0.017870;
end;

varobs DC DG DI DRB DIPCA DRF DIRPF DIRPJ DTAX DYF DPF DNm DS DICM DB DBF;

//DX DIMP        
//steady;
//check;

estimated_params;
//omegaD, uniform_pdf, , , 0.07, 0.2;
//omegaF, uniform_pdf, , , 0.005, 0.025;
//thetaY, uniform_pdf, , , 0.01, 10;
//Omegam, uniform_pdf, , , 0.4, 0.5;
//PDss, uniform_pdf, , , 3, 4.5;  //xxxxxxxxxxxx
//PFss, uniform_pdf, , , 3, 8;    //xxxxxxxxxxxx
//chiBF, uniform_pdf, , , -0.005, -0.0001;  //xxxxxxxxxxxx
//prss, uniform_pdf, , , 0.01, 0.03; //Renzo Orsi,Davi de Raggia, Francesco Turino: Size,trend,and policy implications of the undergroundeconomy
//s, uniform_pdf, , , 1, 1.3; //Renzo Orsi,Davi de Raggia, Francesco Turino: Size,trend,and policy implications of the undergroundeconomy
//gammaMCm, beta_pdf, 0.4760, 0.1;  
//gammaB, uniform_pdf, , , 1, 2;  
//Pss, uniform_pdf, , , 0.1, 2;
//sigma, gamma_pdf, 2, 0.1;
//varphimss, gamma_pdf, 8, 0.5; //xxxxxxxxxxxx
//varphi, gamma_pdf, 9, 0.1;
//thetaW, beta_pdf, 0.65, 0.1;
//theta, beta_pdf, 0.65, 0.05;
//gammaR, beta_pdf, 0.6, 0.15;          //Castro et al (2010)
//gammaPI, normal_pdf, 3, 0.1;       //Castro et al (2010)
//gammaY, beta_pdf, 0.25, 0.05;        //Castro et al (2010)
psiD, gamma_pdf, 1, 0.5;                                                   //Castro et al
psiF, gamma_pdf, 1, 0.5;
omegaRf, beta_pdf, 0.5, 0.25;
psim, uniform_pdf, , , 1.4, 1.5;  //xxxxxxxxxxxx
psiu, uniform_pdf, , , 1.51, 1.7;  //xxxxxxxxxxxx
tau_css, beta_pdf, 0.1594, 0.01;     //Araújo e Ferreira (1999)
tau_lss, beta_pdf, 0.1730, 0.01;     //Araújo e Ferreira (1999)
tau_sss, beta_pdf, 0.105, 0.01;      //Cavalcanti and Silva (2010)
tau_corpss, uniform_pdf, , , 0.25, 0.35;
gammaG, uniform_pdf, , , 0.01, 0.99;
gammaIGm, uniform_pdf, , , 0.01, 0.99;
gammaIGu, uniform_pdf, , , 0.01, 0.99;
gammaT, uniform_pdf, , , 0.01, 0.99; 
gammatau_c, uniform_pdf, , , 0.01, 0.99;
gammatau_corp, uniform_pdf, , , 0.5, 0.99;
gammatau_l, uniform_pdf, , , 0.01, 0.99;
gammatau_s, uniform_pdf, , , 0.01, 0.99;
phiG, uniform_pdf, , , -1, 0;
phiIGm, uniform_pdf, , , -1, 0;
phiIGu, uniform_pdf, , , -1, 0;
phiT, uniform_pdf, , , 0, 1;
phitau_c, uniform_pdf, , , 0, 1;
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
rhoIGm, beta_pdf, 0.5, 0.25;
rhoIGu, beta_pdf, 0.5, 0.25;
rhoT, beta_pdf, 0.5, 0.25;
rhotau_c, beta_pdf, 0.5, 0.25;
rhotau_corp, beta_pdf, 0.5, 0.25;
rhotau_s, beta_pdf, 0.5, 0.25;
rhotau_l, beta_pdf, 0.5, 0.25;
rhom, beta_pdf, 0.5, 0.25;
rhovarphim, beta_pdf, 0.5, 0.25;
rhoPF, beta_pdf, 0.5, 0.25;
rhoRF, beta_pdf, 0.5, 0.25;
rhoYF, beta_pdf, 0.5, 0.25;
//Desvios padrão
stderr e_P, inv_gamma_pdf, 1, inf;
stderr e_L , inv_gamma_pdf, 1, inf;
stderr e_Am, inv_gamma_pdf, 1, inf;
stderr e_Au, inv_gamma_pdf, 1, inf;
stderr e_pr, inv_gamma_pdf, 1, inf;
stderr e_G, inv_gamma_pdf, 1, inf;
stderr e_IGm, inv_gamma_pdf, 1, inf;
stderr e_IGu, inv_gamma_pdf, 1, inf;
stderr e_T, inv_gamma_pdf, 1, inf;
stderr e_tau_c, inv_gamma_pdf, 1, inf;
stderr e_tau_corp, inv_gamma_pdf, 1, inf;
stderr e_tau_s, inv_gamma_pdf, 1, inf;
stderr e_tau_l, inv_gamma_pdf, 1, inf;
stderr e_m, inv_gamma_pdf, 1, inf;
stderr e_varphim, inv_gamma_pdf, 1, inf;
stderr e_PF, inv_gamma_pdf, 1, inf; 
stderr e_RF, inv_gamma_pdf, 1, inf; 
stderr e_YF, inv_gamma_pdf, 1, inf;
stderr e_preco, inv_gamma_pdf, 1, inf;
stderr e_gov, inv_gamma_pdf, 1, inf;
corr e_IGm, e_IGu, 0.5, , , beta_pdf, 0, 0.3, -1, 1;
end;

//identification;
//dynare_sensitivity;

// /*
estimation(datafile=BASE_TRI_14092018, mode_check, smoother, plot_priors = 0, mh_nblocks = 2, mh_drop=0.5, 
mh_jscale=0.3,mh_replic=500000, mode_compute=6, bayesian_irf,irf=20);

shock_decomposition Gr_Nu;

ResultadosBrasil2709 = oo_;
save ResultadosBrasil2709;
// */
