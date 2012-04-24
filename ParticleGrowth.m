function [dCv dCp] = ParticleGrowth(Cv, Cp, CondSink, j, KelvinYN,ROG,Dp)

global modelAtm

%j is the population number
n=modelAtm.NumBins;
KelvinYN = modelAtm.KelvinYN;
CStarBasis = modelAtm.CStarBasis;
alpha = modelAtm.SOA.alphaProd;
dimerbin = 8;

CpAging = zeros(1,n);
CvAging = zeros(1,n);

%k_age_f = 0.0464;
k_age_f = 0.093*1e10*16/(1*3600)*1e-3;  %THESE ONES ARE WORKING
k_age_r = 5e-5*1e10*16/(1*3600)*1e-3; %/s %THESE ONES ARE WORKING


K = 5000; % from a = 0.005 and b = 0.15.... K = b/a^2
k_age_r = 5e-5; %kr is fixed
k_age_f = K*k_age_r;

k_age_f = 0.00001;
k_age_f = 0.025;
k_age_f = 0.0001;
%%%k_age_f = 0.5;
k_age_f2 = k_age_f;

modelAtm.k_age_f = k_age_f;

%if sum(Cp)>1e-5 && Dp > 15*1e-9
if sum(Cp)>1e-5 
    %R_age_f = k_age_f*(Cp(9)/sum(Cp))^2;
    R_age_f = k_age_f*Cp(dimerbin)^2;
    R_age_f2 = k_age_f2*Cp(dimerbin)^2;
    %R_age_r = k_age_r*(Cp(1)/sum(Cp));
    R_age_r = k_age_r*Cp(1);
else
    R_age_f = 0;
    R_age_r = 0;
end


modelAtm.R_age_f = R_age_f;

CpAging(1) = 2*R_age_f-R_age_r;
CpAging(1) = R_age_f - R_age_r;
CpAging(dimerbin) = -1*CpAging(1);
CpAging(dimerbin) = 2*R_age_r - 2*R_age_f;
modelAtm.bla = CpAging;
bla = modelAtm.bla;
if (CpAging(dimerbin)*1 + Cp(dimerbin)) < 1e-5
    CpAging(dimerbin) = 0;
    CpAging(1) = 0;
end

if (CpAging(1)*1 + Cp(1)) < 1e-5
    CpAging(dimerbin) = 0;
    CpAging(1) = 0;
end

UptakeFact = modelAtm.UptakeFact;
%UptakeFact = 1; % NO MASS XFER LIMITATIONS
    if KelvinYN==1
        Kelvin = KelvinTerm(j,Dp);
    else
        Kelvin = 1;     %do not consider Kelvin effect       
    end
    
for i=1:n
    if j==1 %suspended particles
        P(i) = alpha(i)*ROG;  %Production rate of species i (vapor), [=] ug/m3-s
    else
        P(i) = 0;              % Only produce vapors once
    end

    if sum(Cp)<=0 
        MolFrac(i) = 0;
    else
        MolFrac(i)=Cp(i)/(sum(Cp));  %Calculate mole fraction of organic i in particle phase
    end
    
    Phi(i)=UptakeFact*CondSink*(Cv(i)-Kelvin*MolFrac(i)*CStarBasis(i)); % calculte flux assuming no other limitations [=]ug/m3-s    
%    Phi(i) = -Kelvin*MolFrac(i)*CStarBasis(i)*UptakeFact*CondSink; % evap into vaccuum
%   Phi(i) = QPhi(i);
   if Dp*1e9 < 5
       Phi(i) = 0;
   end
   if Cp(i) < 1e-13
       Phi(i) = 0; 
   end
   
  CvAging(i) = 0; % no aging
  dCv(i)=P(i)-Phi(i)+CvAging(i);
  %dCp(i)=Phi(i)+CvAging(i);
  dCp(i) = Phi(i);    
      %Change in gas-phase concentration, or production minus flux to particles, [=] ug/m3-s
      %Change in particle-phase concentration, or flux to particle, [=] ug/m3-s
end
    



