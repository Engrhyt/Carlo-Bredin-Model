%--------------------------------------------------------------------------------------------------------------------------------
   figure
   LP_D=zeros(12,1);
   LP_D_err=zeros(12,1);
   
   LP_H=zeros(12,1);
   LP_H_err=zeros(12,1);
   
   ParcD=0;
   ParcD_err=0;
   
   ParcH=0;
   ParcH_err=0;
   
   for position_test=1:12
     
    %Deuterium
    LP_D(position_test)=mean(collect_point(1:count_collect(position_test,1),position_test,2,1));
    %LP_D_err(position_test)=std(collect_point(1:count_collect(position_test,1),position_test,2,1),1)^2;
    %LP_D_err(position_test)=LP_D_err(position_test)+mean(collect_point(1:count_collect(position_test,1),position_test,3,1))^2;
    %LP_D_err(position_test)=LP_D_err(position_test)^0.5;
    ParcD=ParcD+mean(collect_point(1:count_collect(position_test,1),position_test,1,1));
    %ParcD_err=std(collect_point(1:count_collect(position_test,1),position_test,1,1),1)^2;
    %ParcD_err=ParcD_err+mean(collect_point(1:count_collect(position_test,1),position_test,6,1))^2;
   
    %Hydrogen
    LP_H(position_test)=mean(collect_point(1:count_collect(position_test,2),position_test,2,2));
    %LP_H_err(position_test)=std(collect_point(1:count_collect(position_test,2),position_test,2,2),1)^2;
    %LP_H_err(position_test)=LP_H_err(position_test)+mean(collect_point(1:count_collect(position_test,2),position_test,3,2))^2;
    %LP_H_err(position_test)=LP_H_err(position_test)^0.5;
    ParcH=ParcH+mean(collect_point(1:count_collect(position_test,2),position_test,1,2));
    %ParcH_err=ParcH_err+std(collect_point(1:count_collect(position_test,2),position_test,1,2),1)^2;
    %ParcH_err=ParcH_err+mean(collect_point(1:count_collect(position_test,2),position_test,6,2))^2;
   
   endfor
   ParcD=ParcD/12;
   ParcH=ParcH/12;
   ParcD_err=(ParcD_err/12)^0.5;
   ParcH_err=(ParcH_err/12)^0.5;
   
      hh=errorbar((1:12)*20-140,LP_D(:,1),LP_D_err(:,1)*.5);
   set(hh,'linestyle','none')
   set(hh,'marker','*')
   set(hh,'color','r')
   
   hold on
   
   hh=errorbar((1:12)*20-140,LP_H(:,1),LP_H_err(:,1)*.5);
   set(hh,'linestyle','none')
   set(hh,'marker','*')
   set(hh,'color','b')
   
   xlabel('Position [mm]','fontsize',14)
   ylabel('n^{-}_{LP} [10^{17}m^{-3}]','fontsize',14)
   title(['LP Density Profile'],'fontsize',14)
   grid on
   legend('cs-Deuterium','cs-Hydrogen','location','northwest')
   

annotation('textbox',[0.4 0.85 0.2 0.2],'edgecolor',[1 1 1],'color','r','fontsize',14,'string',['Arc Power ' num2str( ParcD,'%2.2f'  ) '\pm' num2str( ParcD_err,'%1.2f'  )]);
annotation('textbox',[0.4 0.8 0.2 0.2],'edgecolor',[1 1 1],'color','b','fontsize',14,'string',['Arc Power ' num2str( ParcH,'%2.2f' ) '\pm' num2str( ParcH_err,'%1.2f'  )]);

   
   set(gca,'fontsize',14)
   
   ylim([0 2])
   xlim([-200 200])
   saveas(hh,['LPprofile.png'])