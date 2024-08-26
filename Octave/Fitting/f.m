function y=f(b,x)
  dlmwrite('pa.csv',b,',');
  dlmwrite('in.csv',[x x],',');
  system("/home/engrhyt/Desktop/carlofitfast/fit/conver2");
  out=dlmread('out.csv',',');
  y=out(:,2);
  
figure('visible','off')
h=plot(out(:,1),out(:,3).*1000,'g');
hold on
h=plot(out(:,1),out(:,4).*1000,'r');
hold on
h=plot(out(:,1),out(:,5).*1000,'b');
ylabel('Probe Current [mA]','fontsize',14)
xlabel('Probe Voltage [V]','fontsize',14)
xlim([-40 40])
%ylim([-5 30])
grid on
set(gca,'fontsize',14)
saveas(h,['fit.png'])
close
end