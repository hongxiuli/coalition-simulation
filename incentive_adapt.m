function [gamma_i, gamma_temp_i] = incentive_adapt(alpha,beta,phi, noncoal, coal, r, r_S)
  %S = [0,0,1,1,0,0,1,1,0,0];
  e_bar = alpha./beta ; % a n*1 column vector
  E_bar = sum(e_bar);
  %phi = v- theta_rel; % a n*1 column vector

  psi = phi./beta; % a n*1 column vector

  phi_S = sum(phi(coal)); %a scarlar
      
  psi_Sj = phi_S./beta; % a n*1 column vector
  psi_S = sum(psi_Sj(coal));  %a scarlar
  psi_O = sum(psi)-sum(psi(coal)); %a scarlar


  lamda = 1+psi_O+psi_S; %a scarlar
  lamda_i = lamda+psi_Sj+r.*phi.*sum(1./beta(coal))-((1-r).*psi); %vector

  e_O = e_bar - psi./lamda*E_bar; % a nonmember's emissions n*1 column vector
  e_S = e_bar - (psi_Sj+r.*psi)./lamda_i*E_bar; % a member's emissions n*1 column vector

  gamma_all = (e_O-e_S).*[alpha-beta./2.*(e_O+e_S)]-phi./2.*[1/lamda^2-r./lamda_i.^2].*(E_bar^2);
  gamma_i = gamma_all(noncoal);

  %gamma_new = beta/2.*(E_bar^2).*( ((psi_Sj+r.*psi).^2+ r.*psi)./(lamda_i.^2)- (psi.^2+psi)./(lamda^2))
  %gammanew_i = gamma_new(noncoal);

  %first part of the gamma is a function of r
  %gamma_temp = phi./2.*(E_bar^2)./(lamda_i.^3).* ((2.*(psi_Sj+r.*psi)+ 1).*lamda_i) - 2.*((psi_Sj+r.*psi).^2+ r.*psi).*(beta.*sum(1./beta(coal))+1);
  gamma_temp = beta./2.*(E_bar^2)./(lamda_i.^3).* ((((2.*psi.*(psi_Sj+r.*psi)+ psi)).*lamda_i) - 2.*((psi_Sj+r.*psi).^2+ r.*psi).*(phi.*sum(1./beta(coal))+psi));

  gamma_temp_i = gamma_temp(noncoal);
  %((psi_Sj+r.*psi).^2+ r.*psi)%./(lamda_i.^2)

end



