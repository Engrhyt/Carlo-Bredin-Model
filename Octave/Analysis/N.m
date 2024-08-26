 function pdf = W1(b,x,ub,lb)
  
  %constrain
 %for i=1:4 
 %if(b(i)>ub(i))
 %    b(i)=ub(i);
 %endif
 %if(b(i)<lb(i))
 %    b(i)=lb(i);
 %endif
 %endfor
 
  b(4)=abs(b(4));
  %modify gaussian function
  mu=b(1);
  sigma=b(2);
  tau=b(3);
  if license('test', 'Symbolic_Toolbox')
      useVPA = 1; % compute with additional precision
  else
      useVPA = 0;
  end
  if (tau > 0 && sigma > 0) %check for correctness
      % preliminary computations
      if useVPA
        sigmaPerTau = vpa(sigma./tau);
      else
        sigmaPerTau = sigma./tau;
      end
      muMinusX = mu - x;
      normalPart = 1 - normcdf(muMinusX./sigma + sigmaPerTau);
      expPart = muMinusX./tau + (sigmaPerTau.^2)./2;
      if (useVPA)
          pdf = double((1/tau).*normalPart.*exp(expPart));
      else
          %if (expPart > 25)
              %warning('The combination of the input values may result in numerical instability. Please change your units to make (mu-x)/tau < 25 and sigma/tau < 5. Then please re-run the script')
          %end
          pdf = (1/tau).*normalPart.*exp(expPart);
      end
  else
      pdf = zeros(length(x), 1);
  end
  pdf(pdf==Inf) = 0;
  pdf=pdf./max(pdf).*b(4);
  
end