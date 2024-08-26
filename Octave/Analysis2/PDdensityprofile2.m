
  y_data_type_index=3;

  figure('visible','on')
 PD=zeros(14,4);
 PD_std=zeros(14,4);
for plasma_type=3:4

 power_index=4
       count_t=0;
      arc_st=0;
      data_value=zeros(100,1);
      for position_test=6
       for K=1:data_count(position_test,power_index,plasma_type);
         arc_power=(data(K,position_test,2,power_index,plasma_type))*1e-17;
         if(arc_power>2.3&&arc_power<2.7)
          count_t=count_t+1;
          arc_st=arc_st+arc_power;
         end
       endfor
      end
       arc_st=arc_st/count_t;
       
    arc_av=0;
    arc_sq=0;
    arc_i=0;
    for position_test=1:12
      count_t=0;
      arc_avl=0;
      data_value=zeros(100,1);
       for K=1:data_count(position_test,power_index,plasma_type);
         arc_power=(data(K,position_test,2,power_index,plasma_type))*1e-17;
         if(arc_power>2.3&&arc_power<2.7)
          aud=data(K,position_test,y_data_type_index,power_index,plasma_type);
          if(aud>0)
          count_t=count_t+1;
          arc_i=arc_i+1;
          arc_avl=arc_avl+arc_power;
          arc_av=arc_av+arc_power;
          arc_sq=arc_sq+arc_power.^2;
          data_value(count_t)=(aud).*1.259*2/50*1000;
          end
         end
       endfor
       arc_avl=arc_avl/count_t;
       PD(1+position_test,plasma_type)=mean(data_value(1:count_t))/arc_avl*arc_st;
       PD_std(1+position_test,plasma_type)=std(data_value(1:count_t))/arc_avl*arc_st;

  endfor%position_test
     xxx=[-175 (1:12)*20-140 175];
     if(plasma_type==3)
        ax=errorbar(xxx,abs(PD(:,plasma_type))./k_D,PD_std(:,plasma_type)./k_D);
        hold on
     set(ax,'linestyle','-')
     set(ax,'marker','*')
     set(ax,'color','r')
       elseif(plasma_type==4)
        ax=errorbar(xxx,abs(PD(:,plasma_type))./k_H,PD_std(:,plasma_type)./k_H);
        hold on
     set(ax,'linestyle','-')
     set(ax,'marker','*')
     set(ax,'color','b')
     end
     

if(plasma_type==3)
annotation('textbox',[0.15 0.80 0.2 0.2],'edgecolor',[1 1 1],'color','r','fontsize',14,'string',['Deuterium                CRD Density ' num2str( arc_av/arc_i,'%2.2f'  ) '\pm' num2str( ((arc_sq/arc_i)-(arc_av/arc_i)^2)^0.5,'%1.2f') '[10^{17}m^{-3}]']);
elseif(plasma_type==4)
annotation('textbox',[0.15 0.85 0.2 0.2],'edgecolor',[1 1 1],'color','b','fontsize',14,'string',['Hydrogen                 CRD Density ' num2str( arc_av/arc_i,'%2.2f' ) '\pm' num2str( ((arc_sq/arc_i)-(arc_av/arc_i)^2)^0.5,'%1.2f') '[10^{17}m^{-3}]']);
endif

endfor%plasma_type

 ylim([0 3])
 xlim([-200 200])

 ylabel('PD Negative ion Density [10^{17}m^{-3}]','fontsize',14)
 xlabel('Position [mm]','fontsize',14)
 title(['PD integral Profile Peak ' num2str(pd_peak)],'fontsize',14)
 %legend('Cs-Deuterium','Cs-Hydrogen','location','northeast')
 grid on
 set(gca,'fontsize',14)
 saveas(ax,['PDprofileP' num2str(pd_peak) '.png'])
 %close