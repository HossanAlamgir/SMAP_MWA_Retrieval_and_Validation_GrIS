function [TBH_incoh,TBV_incoh,TBH_coh,TBV_coh,NPR_incoh,NPR_coh,NPR_ref]=IceSheet_EM_Model_nlayer(f,T0,d,eps,TBH,TBV,theta_0)

% This code calculates the TB for a n layer media using both 
% the coherent (wave) approach and incoherent radiative transfer (RT)
% approaches. This model works for both Melt and Freeze seasons. In this
% model, there are n layers and n-1 boundaries.
%%%---Assumptions for this model:
%%% 1. Volume scattering is negligible in both wet and dry snow layers 
%%%   in the frequency of operation (Zeroth Order RT).
%%% 2. Surfaces are electrically smooth in the frequency of operation.
%%% 3. Atmospheric attenuation is negligible in the frequency of operation.
%%% 4. All layers are homogenous medium with uniform 
%%%   temperature profile.
%%% 5. Last layer from the top looks like a semi-infinite medium in the frequency
%%%   of operation.
%%%
%%%--- INPUTS:
%%%--- f: frequency of operation in GHz
%%%--- T0: Physical Temperature of each layer in vector format (K)
%%%--- d: location of boundaties in vector format in cm (d>0)
%%%--- eps: relative permittivity of each layer in vector format
%%%--- TBH: SMAP H-pol brightness temperature of freeze season (K)
%%%--- TBV: SMAP V-pol brightness temperature of freeze season (K)
%%%--- theta_0: incident (observation) angle in degrees 
%%%
%%%
%%%--- OUTPUTS:
%%%--- TBH_incoh: incoherent H-pol brightness temperature (K)
%%%--- TBV_incoh: incoherent V-pol brightness temperature (K)
%%%--- TBH_coh: coherent H-pol brightness temperature (K)
%%%--- TBV_coh: coherent V-pol brightness temperature (K)
%%%--- NPR_incoh: incoherent normalized polarization ratio (NPR)
%%%--- NPR_coh: coherent normalized polarization ratio (NPR)
%%%--- NPR_ref: refernce NPR calculated from freeze season TB

%author: Mohammad Mousavi (MM), April 2020, NASA-JPL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=3e8;
f=f.*1e9;
k0=2*pi.*(f)./c;




NPR_ref=(TBV-TBH)/(TBV+TBH);% it's just a NPR value as a reference based on SMAP TB inputs during WREF

n=length(eps);
d=d.*1e-2;

alfa=-k0.*imag(sqrt(eps));
kan=2.*alfa;%wet snow power absoprtion coeff

thetan=zeros(n,1);
sin_thetan=zeros(n,1);
cos_thetan=zeros(n,1);
thetan(1)=theta_0*(pi/180);%SMAP Incident angle
sin_thetan(1)=sin(thetan(1));
cos_thetan(1) = sqrt(1 - sin_thetan(1).^2);
for ii=2:n
    sin_thetan(ii) = sqrt(eps(ii-1))./sqrt(eps(ii)).*sin(thetan(ii-1));
    cos_thetan(ii) = sqrt(1 - sin_thetan(ii).^2);
    thetan(ii)=(asin(sqrt(eps(ii-1))./sqrt(eps(ii)).*sin(thetan(ii-1))));
end

rhoh=zeros(n-1,1);
rhov=zeros(n-1,1);
rh=zeros(n-1,1);
rv=zeros(n-1,1);

for ii=1:n-1


    
    
    rhoh(ii) = (sqrt(eps(ii)).*cos(thetan(ii))-sqrt(eps(ii+1)).*cos(thetan(ii+1))) ./...
        (sqrt(eps(ii)).*cos(thetan(ii))+sqrt(eps(ii+1)).*cos(thetan(ii+1)));
    rhov(ii) = (sqrt(eps(ii)).*cos(thetan(ii+1))-sqrt(eps(ii+1)).*cos(thetan(ii))) ./...
        (sqrt(eps(ii)).*cos(thetan(ii+1))+sqrt(eps(ii+1)).*cos(thetan(ii)));
    
    
    
    rh(ii)=(abs(rhoh(ii))).^2;
    rv(ii)=(abs(rhov(ii))).^2;


    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n==2
    TBH_incoh=(1-rh(1))*T0(2);
    TBV_incoh=(1-rv(1))*T0(2);    
elseif n>2
mat_c_h=zeros(2*(n-2),2*(n-2));
mat_c_v=zeros(2*(n-2),2*(n-2));
T_y_h=zeros(2*(n-2),1);
T_y_v=zeros(2*(n-2),1);



thetan=real(thetan);

tmp2_h=zeros(size(mat_c_h,1),1);
tmp2_v=zeros(size(mat_c_v,1),1);
tmp2_h(1,1)=-rh(1).*exp(-kan(2)*d(2)*sec(thetan(2)));
tmp2_v(1,1)=-rv(1).*exp(-kan(2)*d(2)*sec(thetan(2)));


for ii=2:size(mat_c_h,1)
    if mod(ii,2)==0
       ind_tmp=(ii/2) +1;
       tmp2_h(ii)= -rh(ind_tmp).*exp(kan(ind_tmp)*(-d(ind_tmp)+d(ind_tmp-1))*sec(thetan(ind_tmp)));
       tmp2_v(ii)= -rv(ind_tmp).*exp(kan(ind_tmp)*(-d(ind_tmp)+d(ind_tmp-1))*sec(thetan(ind_tmp)));
    elseif mod(ii,2)~=0
       ind_tmp=floor(ii/2) +1;
       tmp2_h(ii)= -rh(ind_tmp).*exp(-kan(ind_tmp+1)*(-d(ind_tmp)+d(ind_tmp+1))*sec(thetan(ind_tmp+1)));
       tmp2_v(ii)= -rv(ind_tmp).*exp(-kan(ind_tmp+1)*(-d(ind_tmp)+d(ind_tmp+1))*sec(thetan(ind_tmp+1)));

    end
end




tmp1_h=zeros(size(mat_c_h,1)-1,1);
tmp1_v=zeros(size(mat_c_h,1)-1,1);

for ii=1:size(tmp1_h,1)
    if mod(ii,2)~=0
        tmp1_h(ii)=1;
        tmp1_v(ii)=1;
    elseif mod(ii,2)==0
        ind_tmp=ii/2 +1;
        tmp1_h(ii)=-(1-rh(ind_tmp)).*exp(kan(ind_tmp).*(-d(ind_tmp)+d(ind_tmp-1)).*sec(thetan(ind_tmp)));
        tmp1_v(ii)=-(1-rv(ind_tmp)).*exp(kan(ind_tmp).*(-d(ind_tmp)+d(ind_tmp-1)).*sec(thetan(ind_tmp)));
    end   
    
end

tmp3_h=zeros(size(mat_c_h,1)-1,1);
tmp3_v=zeros(size(mat_c_h,1)-1,1);

for ii=1:size(tmp3_h,1)
    if mod(ii,2)~=0
        tmp3_h(ii)=1;
        tmp3_v(ii)=1;
    elseif mod(ii,2)==0
        ind_tmp=ii/2 +1;
        tmp3_h(ii)=-(1-rh(ind_tmp)).*exp(-kan(ind_tmp+1).*(-d(ind_tmp)+d(ind_tmp+1)).*sec(thetan(ind_tmp+1)));
        tmp3_v(ii)=-(1-rv(ind_tmp)).*exp(-kan(ind_tmp+1).*(-d(ind_tmp)+d(ind_tmp+1)).*sec(thetan(ind_tmp+1)));
    end
  
end



mat_c_h=diag(tmp1_h,-1)+diag(tmp2_h)+diag(tmp3_h,1);
mat_c_v=diag(tmp1_v,-1)+diag(tmp2_v)+diag(tmp3_v,1);


% T0(1)=0;

T_y_h(1)=-rh(1).*(1-exp(-kan(2)*d(2)*sec(thetan(2)))).*T0(2);
T_y_v(1)=-rv(1).*(1-exp(-kan(2)*d(2)*sec(thetan(2)))).*T0(2);

N=length(eps)-1;
T_y_h(end)=rh(N).*(1-exp(kan(N).*(-d(N)+d(N-1)).*sec(thetan(N)))).*T0(N) +...
            (1-rh(N)).*T0(N+1);
T_y_v(end)=rv(N).*(1-exp(kan(N).*(-d(N)+d(N-1)).*sec(thetan(N)))).*T0(N) +...
            (1-rv(N)).*T0(N+1);


for ii=2:size(T_y_h,1)-1
    if mod(ii,2)==0
        ind_tmp=ii/2 +1;
        T_y_h(ii)=rh(ind_tmp).*(1-exp(kan(ind_tmp).*(-d(ind_tmp)+d(ind_tmp-1)).*sec(thetan(ind_tmp)))).*T0(ind_tmp) +...
            (1-rh(ind_tmp)).*(1-exp(-kan(ind_tmp+1).*(-d(ind_tmp)+d(ind_tmp+1)).*sec(thetan(ind_tmp+1)))).*T0(ind_tmp+1);
        T_y_v(ii)=rv(ind_tmp).*(1-exp(kan(ind_tmp).*(-d(ind_tmp)+d(ind_tmp-1)).*sec(thetan(ind_tmp)))).*T0(ind_tmp) +...
            (1-rv(ind_tmp)).*(1-exp(-kan(ind_tmp+1).*(-d(ind_tmp)+d(ind_tmp+1)).*sec(thetan(ind_tmp+1)))).*T0(ind_tmp+1);
    elseif mod(ii,2)~=0
        ind_tmp=floor(ii/2) +1;
        T_y_h(ii)=rh(ind_tmp).*(1-exp(-kan(ind_tmp+1).*(-d(ind_tmp)+d(ind_tmp+1)).*sec(thetan(ind_tmp+1)))).*T0(ind_tmp+1) +...
            (1-rh(ind_tmp)).*(1-exp(kan(ind_tmp).*(-d(ind_tmp)+d(ind_tmp-1)).*sec(thetan(ind_tmp)))).*T0(ind_tmp);
        T_y_v(ii)=rv(ind_tmp).*(1-exp(-kan(ind_tmp+1).*(-d(ind_tmp)+d(ind_tmp+1)).*sec(thetan(ind_tmp+1)))).*T0(ind_tmp+1) +...
            (1-rv(ind_tmp)).*(1-exp(kan(ind_tmp).*(-d(ind_tmp)+d(ind_tmp-1)).*sec(thetan(ind_tmp)))).*T0(ind_tmp);
    end
  
end


% warning('off','all')
AB_vec_h=mat_c_h \ T_y_h;
AB_vec_v=mat_c_v \ T_y_v;
% AB_vec_h=inv(mat_c_h) * T_y_h;
% AB_vec_v=inv(mat_c_v) * T_y_v;
% AB_vec_h=pinv(mat_c_h,rcond(mat_c_h)) * T_y_h;
% AB_vec_v=pinv(mat_c_v,rcond(mat_c_v)) * T_y_v;

% AB_vec_h=lsqminnorm(mat_c_h, T_y_h);
% AB_vec_v=lsqminnorm(mat_c_v,T_y_v);
% AB_vec_h(isnan(AB_vec_h))=0;
% AB_vec_v(isnan(AB_vec_v))=0;

TBH_incoh=(AB_vec_h(1).*exp(-kan(2)*d(2)*sec(thetan(2))) + (1-exp(-kan(2)*d(2)*sec(thetan(2)))).*T0(2))*(1-rh(1));
TBV_incoh=(AB_vec_v(1).*exp(-kan(2)*d(2)*sec(thetan(2))) + (1-exp(-kan(2)*d(2)*sec(thetan(2)))).*T0(2))*(1-rv(1));
end
% if isnan(TBH_incoh)
%     TBH_incoh=-T_y_h(1);
% end
% 
% if isnan(TBV_incoh)
%     TBV_incoh=-T_y_v(1);
% end


NPR_incoh=(TBV_incoh-TBH_incoh)./(TBV_incoh+TBH_incoh);

%%%% Coherent part will be updated!
TBH_coh=NaN;
TBV_coh=NaN;
NPR_coh=NaN;
NPR_ref=(TBV-TBH)./(TBV+TBH);
end



