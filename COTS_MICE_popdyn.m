function [cots_num, fgc_biomass,fgc_eatenbycots, rho,fgc_effcotsmort, fgc_growth] = COTS_MICE_popdyn( years  )

% Simple COTS-coral population dynamics model v2
% Implementation of a published MICE model by Morello et al. 2014, MEPS
% Karlo Hock, University of Queensland, v1 2014; v2 2018

% Variable and parameter names used in Morello et al. given in comments


cots_num=zeros(years+1,3);%N per year per age, number of COTS
cots_num(1,1)=0.505*exp(2*2.56);%N per year of year 0 COTS
cots_num(1,2)=0.505*exp(2.56);%N per year of year 1 COTS
cots_num(1,3)=0.505;%N per year of year 2+ COTS

resid=zeros(years,1);%epsilon per year
resid(3,1)=4.307;%epsilon for year 1996

imig=zeros(years+1,1);%eta per year
imig(1,1)=4.292;%eta for year 1994

fgc_biomass=zeros(years,1);%C f per year
sgc_biomass=zeros(years,1);%C m per year
fgc_biomass(1,1)=2500;%C f initial, same as carrying capacity K f
sgc_biomass(1,1)=500;%C m initial, same as carrying capacity K m
rho=zeros(years,1);%rho, switch function
fgc_eatenbycots=zeros(years,1);%Q f per year
sgc_eatenbycots=zeros(years,1);%Q m per year
cots_selfrec=zeros(years+1,1);%R per year
cots_selfrec(1,1)=1;%R f for 1994, using the term from equation 1a
fgc_effcotsmort=zeros(years,1);%f(C f per year)
fgc_growth=zeros(years,1);%fast coral growth, second term in equation 2
sgc_growth=zeros(years,1);%slow coral growth, second term in equation 3

for yr=1:years
    rho(yr)=1+exp((-5)*(fgc_biomass(yr)/2500));
    fgc_eatenbycots(yr)=(1-rho(yr))*(0.129*(cots_num(yr,2)+cots_num(yr,3))*fgc_biomass(yr))/(1+exp((-1*(cots_num(yr,2)+cots_num(yr,3)))/10));%
    sgc_eatenbycots(yr)=rho(yr)*(0.268*(cots_num(yr,2)+cots_num(yr,3))*sgc_biomass(yr))/(1+exp((-1*(cots_num(yr,2)+cots_num(yr,3))/8)));
    %if yr>1
    fgc_growth(yr)=0.5*fgc_biomass(yr)*(1-(fgc_biomass(yr)/2500));
    sgc_growth(yr)=0.1*sgc_biomass(yr)*(1-(sgc_biomass(yr)/500));
    %end
    fgc_biomass(yr+1)=fgc_biomass(yr)+fgc_growth(yr)+fgc_eatenbycots(yr);%*23.5
    
    sgc_biomass(yr+1)=sgc_biomass(yr)+sgc_growth(yr)+sgc_eatenbycots(yr);
    %if yr>1
        cots_selfrec(yr+1,1)=(cots_num(yr,3)/cots_num(yr,3))*exp(resid(yr));
   %end
    cots_num(yr+1,1)=cots_selfrec(yr+1)+exp(imig(yr+1));
    fgc_effcotsmort(yr,1)=1-(0.258*(fgc_biomass(yr)/(1+fgc_biomass(yr))));
    cots_num(yr+1,2)=(cots_num(yr,1)*exp(-1*fgc_effcotsmort(yr)*2.56));
    cots_num(yr+1,3)=(cots_num(yr,2)+cots_num(yr,3))*exp(-1*fgc_effcotsmort(yr)*2.56);
    
end
end