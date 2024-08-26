function pdf=M(b,x,ub,lb)
  
  pdf=abs(b(3)).*exp(-abs(b(2)).*(x-b(1)).^2);
  pdf(pdf==Inf) = 0;
  
end