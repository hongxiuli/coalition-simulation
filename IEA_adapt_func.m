function [intcoal,extcoal,coalition,E_S,E,W] = IEA_adapt_func(alpha,beta,phi, index, n)
  %s is the coalition vector, indicating which countries are in the coaltion
  s = size(index,2);
  e_bar = alpha./beta ; % a n*1 column vector
  E_bar = sum(e_bar);% a scarlar

  psi = phi./beta; % a n*1 column vector
  coalition = zeros(size(index));
  intcoal = zeros(size(index));
  extcoal = zeros(size(index));
  E_S = zeros(size(index,1),1);
  E = zeros(size(index,1),1);
  W = zeros(size(index,1),1);

  %Do this loop for each row in the index matrix
  if (1<s) && (s<n)
    for x = 1:nchoosek(n,s)
      index_x = index(x,:)'; %index_x is a column vector of a selection of s countries
      
      phi_S = sum(phi(index_x)); %a scarlar
      
      psi_Sj = phi_S./beta; % a n*1 column vector
      psi_S = sum(psi_Sj(index_x));  %a scarlar
      psi_O = sum(psi)-sum(psi(index_x)); %a scarlar
      
      lamda1 = 1+psi_O+psi_S; %a scarlar
      
      e_S = e_bar - psi_Sj./lamda1*E_bar; % a member's emissions n*1 column vector
      e_O = e_bar - psi/lamda1*E_bar; % a nonmember's emissions n*1 column vector
      e = e_O;
      e(index_x) = e_S(index_x);
      E_S_x = sum(e_S(index_x)) ;
      E_x = E_bar/lamda1;%a scarlar
      w = alpha.*e-beta/2.*e.^2-1/2*phi.*E_x^2;
      
      %Internal stability test for each signatories in the index_x
      %generate internal stab column vector. 
      intstab = zeros(size(index_x,1),1);

      for y = 1:size(index_x,1)
        index_y = index_x;
        index_y(y)=[]; %the y member is now out of the coalition
        
        phi_S_y = sum(phi(index_y)); %a scarlar
    
        psi_Sj_y = phi_S_y./beta; % a n*1 column vector
        psi_S_y = sum(psi_Sj_y(index_y));  %a scarlar
        psi_O_y = sum(psi)-sum(psi(index_y)); %a scarlar
    
        lamda1_y = 1+psi_O_y+psi_S_y; %a scarlar
        %lamda2_y = 1+psi_O_y+psi_S_y+2*psi-psi_Sj_y-sum(1./beta).*phi;  %a n*1 column vector
        E_y = E_bar/lamda1_y; %a scarlar
        
        %the incentive for the y element in the index_x to deviate from the
        %IEA
        gamma_y = (psi*E_y-psi_Sj*E_x).*(alpha-beta/2.*(2*e_bar-psi*E_y-psi_Sj*E_x))...
                -1/2*phi.*(E_x^2-E_y^2);
        
        % if the incentive is positive for the element, then this element
        % is stable.
        if (gamma_y(index_x(y)) >= 0)
           intstab(y) = 1;
        end
      end
      
      if isequal(intstab,ones(size(intstab)))
         intcoal(x,:) = index_x;
      end
      
      % External stability test for each nonsgnatories
      ns = n-size(index_x,1);
      extstab = zeros(ns,1);
      nslist = [1:n]'; %all countries list
      nslist(index_x) = []; % all nonsignatories
      for z = 1:ns
        index_z = nslist;
        index_x_z = [index_x;index_z(z)]; %all counterfactual signatories
        index_z(z) = []; % all counterfactual nonsignatories
        
        phi_S_z = sum(phi(index_x_z)); %a scarlar
    
        psi_Sj_z = phi_S_z./beta; % a n*1 column vector
        psi_S_z = sum(psi_Sj_z(index_x_z));  %a scarlar
        psi_O_z = sum(psi)-sum(psi(index_x_z)); %a scarlar
        lamda1_z = 1+psi_O_z+psi_S_z; %a scarlar
        %lamda2_y = 1+psi_O_y+psi_S_y+2*psi-psi_Sj_y-sum(1./beta).*phi;  %a n*1 column vector
        E_z = E_bar/lamda1_z; %a scarlar
        
        %the incentive for the y element in the index_x to deviate from the IEA
        gamma_z = (psi*E_x-psi_Sj_z*E_z).*(alpha-beta/2.*(2*e_bar-psi*E_x-psi_Sj_z*E_z))...
                -1/2*phi.*(E_z^2-E_x^2);
        % if the incentive is negative for the element, then this element
        % is external stable.
        if (gamma_z(nslist(z)) <= 0)
           extstab(z) = 1;
        end   
      end
    
      if isequal(extstab,ones(size(extstab)))
         extcoal(x,:) = index_x;
      end
  
      if isequal(intstab,ones(size(intstab))) && isequal(extstab,ones(size(extstab)))
      coalition(x,:) = index_x';
      E_S(x,:) = E_S_x;
      E(x,:) = E_x;
      W(x,:) = sum(w);
      end
    end
  elseif s==1
      psi_O = sum(psi); %a scarlar
      lamda1 = 1+psi_O+0; %a scarlar
      
      e_O = e_bar - psi./lamda1.*E_bar; % a nonmember's emissions n*1 column vector
      e = e_O;
      E = E_bar/lamda1;%a scarlar
      w = alpha.*e-beta/2.*e.^2-1/2*phi.*E.^2;
      W = sum(w);  
  else
      phi_S = sum(phi); %a scarlar
       
      psi_Sj = phi_S./beta; % a n*1 column vector
      psi_S = sum(psi_Sj);  %a scarlar
      psi_O = 0; %a scarlar
      
      lamda1 = 1+psi_O+psi_S; %a scarlar
       
      e_S = e_bar - psi_Sj./lamda1*E_bar; % a member's emissions n*1 column vector
      e = e_S;
      E_S = sum(e_S);
      E = E_bar/lamda1;%a scarlar
      w = alpha.*e-beta/2.*e.^2-1/2*phi.*E.^2;
      W = sum(w);
      
      intstab = zeros(size(index,1),1);
      
      %test interal stable of the grand coalition
      for y = 1: n
        index_y = index;
        index_y(y)=[]; %the y member is now out of the coalition
        
        phi_S_y = sum(phi(index_y)); %a scarlar
    
        psi_Sj_y = phi_S_y./beta; % a n*1 column vector
        psi_S_y = sum(psi_Sj_y(index_y));  %a scarlar
        psi_O_y = 0; %a scarlar
    
        lamda1_y = 1+psi_O_y+psi_S_y; %a scarlar
        %lamda2_y = 1+psi_O_y+psi_S_y+2*psi-psi_Sj_y-sum(1./beta).*phi;  %a n*1 column vector
        E_y = E_bar/lamda1_y; %a scarlar
        
        %the incentive for the y element in the index_x to deviate from the
        %IEA
        gamma_y = (psi*E_y-psi_Sj*E).*(alpha-beta/2.*(2*e_bar-psi*E_y-psi_Sj*E))...
                -1/2*phi.*(E^2-E_y^2);
        % if the incentive is positive for the element, then this element
        % is stable.
        if (gamma_y(index(y)) >= 0)
           intstab(y) = 1;
        end
      end
      
      if isequal(intstab,ones(size(intstab)))
         intcoal = index';
         coalition = index';
      end   
  end
end
%We have e_S and e_O, given to the S coaltion. 
%for the first s countries, if gamma>=0, then we have internal stability
%for the last n-s countries, if gamma<=0, then we have external stability
%stability




