%close all
color_code1=['*r' ; '*g' ; '*b' ; '*k'];
color_code2=['or' ; 'og' ; 'ob' ; 'ok'];
%----------------------------------------------plot density profile x-position y-parc z-param---------------------------------------------------------

collect_point=zeros(10,12,6,2);
color_indicator=zeros(10,12,1,2);
count_collect=zeros(12,2);

ksum_D=zeros(12,2);
ksum_H=zeros(12,2);

for plasma_type=3:4
disp('here')
  y_data_type_index=3;

 if(plasma_type==1)
  name='HwoCs';
  mi=1.67262192369*10^(-27);
 elseif(plasma_type==2)
  name='DwoCs';
  mi=2.014*1.67262192369*10^(-27);
 elseif(plasma_type==4)
  name='HwCs';
  mi=1.67262192369*10^(-27);
 elseif(plasma_type==3)
  name='DwCs';
  mi=2.014*1.67262192369*10^(-27);
 endif

  for position_test=1:12

   insert_name='Negative ion'

   %figure('name',[name ' ' insert_name], 'visible','off')
   point_count=0;

    power_limit(3)=4;
    power_limit(4)=3;
    power_limit(1)=4;
    power_limit(2)=4;
    for power_index=4%power_limit(plasma_type)
     
     %data_count,position_test,data_type_index,power_index,plasma_type
     if(data_count(position_test,power_index,plasma_type)>3)
       data_start_det=false;
       data_stop_det=false;
       LHJ_min=0;
       LHJ_max=0;
       for data_index=1:data_count(position_test,power_index,plasma_type)
        if(data(data_index,position_test,9,power_index,plasma_type)>min_time&&!data_start_det)
         LHJ_min=data_index;
         data_start_det=true;
        endif
        if(data(data_index,position_test,9,power_index,plasma_type)>max_time&&!data_stop_det)
         LHJ_max=data_index;
         data_start_det=true;
        endif
       endfor
       if(LHJ_min==0||LHJ_max==0)
        break
       endif
       KK=LHJ_min:LHJ_max;

       power_value = data(KK,position_test,1,power_index,plasma_type)*data_scale(1);
       
       %LP
       H1=data(KK,position_test,6,power_index,plasma_type); %np
       H2=abs(data(KK,position_test,7,power_index,plasma_type)); %as
       HHH=10 .^H1.*H2./(H2+1);
       HHH=HHH/10^17;
       disp('not here')
       
       %PD
       if(plasma_type==3)
        k=k_D;
        k2=k_D1;
       elseif(plasma_type==4)
        k=k_H;
        k2=k_H1;
       endif
       GGG=((data(KK,position_test,y_data_type_index,power_index,plasma_type)*(1.259*2/50*1000)))/k+k2/350;
       
       %keep
       Q1=data(KK,position_test,1,power_index,plasma_type)*data_scale(1);
       Q2=HHH;
       Q2_err=0;
       Q4=GGG;
       Q4_err=0;
       Q5=0;
       for iii=1:length(Q1)
         pass_Q=0;
         %if((abs(mean(Q2)-Q2(iii))/mean(Q2))<0.7)
         if(Q2(iii)<1.5*mean(Q2)&&plasma_type==4)
          pass_Q=1;
         endif
         if(Q2(iii)<0.65*mean(Q2)&&plasma_type==3)
          pass_Q=1;
         endif
         if(pass_Q)
         point_count=point_count+1;
         collect_point(point_count,position_test,1,plasma_type-2)=Q1(iii);
         collect_point(point_count,position_test,2,plasma_type-2)=Q2(iii);
         collect_point(point_count,position_test,3,plasma_type-2)=0;%Q2_err(iii);
         collect_point(point_count,position_test,4,plasma_type-2)=Q4(iii);
         collect_point(point_count,position_test,5,plasma_type-2)=0;%Q4_err(iii);
         collect_point(point_count,position_test,6,plasma_type-2)=0;%Q5(iii);
         color_indicator(point_count,position_test,3,plasma_type-2)=power_index;
         count_collect(position_test,plasma_type-2)=point_count;
         endif
       endfor %iii
     
     endif

    endfor %power_index

 endfor %position_test

endfor %plasma_type
%-----------------------------------------------------------------------------------------------------------------------------------------
     
     for position_test=1:12;
     
     %Deuterium
     A_D=squeeze(collect_point(1:count_collect(position_test,1),position_test,1,1));
     PD_D=squeeze(collect_point(1:count_collect(position_test,1),position_test,4,1));
     LP_D=squeeze(collect_point(1:count_collect(position_test,1),position_test,2,1));
     p_D=squeeze(color_indicator(1:count_collect(position_test,1),position_test,3,1));
     %Hydrogen
     A_H=squeeze(collect_point(1:count_collect(position_test,2),position_test,1,2));
     PD_H=squeeze(collect_point(1:count_collect(position_test,2),position_test,4,2));
     LP_H=squeeze(collect_point(1:count_collect(position_test,2),position_test,2,2));
     p_H=squeeze(color_indicator(1:count_collect(position_test,2),position_test,3,2));

     f=@(b,x) b(1).*x+b(2);
     
     kk_D=[0 0];
     kk_H=[0 0];
     
     if(length(PD_D)<3)
      break
     end
     
     [kk_D,R,J,CovB,MSE]=nlinfit(PD_D,LP_D,f,kk_D,opts);
     kk_err_D=sqrt(diag(CovB));
     [kk_H,R,J,CovB,MSE]=nlinfit(PD_H,LP_H,f,kk_H,opts);
     kk_err_H=sqrt(diag(CovB));
     
     ksum_D(position_test,1)=kk_D(1);
     ksum_H(position_test,1)=kk_H(1);
     ksum_D(position_test,2)=kk_err_D(1);
     ksum_H(position_test,2)=kk_err_H(1);
     
     f_x=0:0.1:4;
     fvecD=f(kk_D,f_x);
     fvecH=f(kk_H,f_x);
     
     if(position_test==7)
     %if(1)
     figure('name',[name ' ' insert_name], 'visible','on')
     if(0)
     pmean1=zeros(100,1);pmean2=zeros(100,1);pmean3=zeros(100,1);pmean4=zeros(100,1);
     mean1=zeros(100,1);mean2=zeros(100,1);mean3=zeros(100,1);mean4=zeros(100,1);
     co_p1=0;co_p2=0;co_p3=0;co_p4=0;
     for ie=1:length(PD_D)
       if(p_D(ie)==1)
       co_p1=co_p1+1;
       mean1(co_p1)=LP_D(ie);
       pmean1(co_p1)=PD_D(ie);
       elseif(p_D(ie)==2)
       co_p2=co_p2+1;
       mean2(co_p2)=LP_D(ie);
       pmean2(co_p2)=PD_D(ie);
       elseif(p_D(ie)==3)
       co_p3=co_p3+1;
       mean3(co_p3)=LP_D(ie);
       pmean3(co_p3)=PD_D(ie);
       elseif(p_D(ie)==4)
       co_p4=co_p4+1;
       mean4(co_p4)=LP_D(ie);
       pmean4(co_p4)=PD_D(ie);
       endif
     endfor
     hh=errorbar(mean(pmean1(1:co_p1)),mean(mean1(1:co_p1)),std(mean1(1:co_p1)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean2(1:co_p2)),mean(mean2(1:co_p2)),std(mean2(1:co_p2)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean3(1:co_p3)),mean(mean3(1:co_p3)),std(mean3(1:co_p3)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean4(1:co_p4)),mean(mean4(1:co_p4)),std(mean4(1:co_p4)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
          pmean1=zeros(100,1);pmean2=zeros(100,1);pmean3=zeros(100,1);pmean4=zeros(100,1);
     mean1=zeros(100,1);mean2=zeros(100,1);mean3=zeros(100,1);mean4=zeros(100,1);
     co_p1=0;co_p2=0;co_p3=0;co_p4=0;
     for ie=1:length(PD_H)
       if(p_H(ie)==1)
       co_p1=co_p1+1;
       mean1(co_p1)=LP_H(ie);
       pmean1(co_p1)=PD_H(ie);
       elseif(p_H(ie)==2)
       co_p2=co_p2+1;
       mean2(co_p2)=LP_H(ie);
       pmean2(co_p2)=PD_H(ie);
       elseif(p_H(ie)==3)
       co_p3=co_p3+1;
       mean3(co_p3)=LP_H(ie);
       pmean3(co_p3)=PD_H(ie);
       elseif(p_H(ie)==4)
       co_p4=co_p4+1;
       mean4(co_p4)=LP_H(ie);
       pmean4(co_p4)=PD_H(ie);
       endif
     endfor
     hh=errorbar(mean(pmean1(1:co_p1)),mean(mean1(1:co_p1)),std(mean1(1:co_p1)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean2(1:co_p2)),mean(mean2(1:co_p2)),std(mean2(1:co_p2)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean3(1:co_p3)),mean(mean3(1:co_p3)),std(mean3(1:co_p3)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean4(1:co_p4)),mean(mean4(1:co_p4)),std(mean4(1:co_p4)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
   end  
     hh=plot(PD_D,LP_D,'.r')
     hold on
     hh=plot(PD_H,LP_H,'.b')
     hold on
     hh=plot(f_x,fvecD,'r','linewidth',1.5)
     hold on
     hh=plot(f_x,fvecH,'b','linewidth',1.5)
     
     xlim([0 2])
     if(pd_peak==3)
      xlim([0 1.5])
     endif
     ylim([0 2])
     
     annotation('textbox',[0.15 0.85 0.2 0.2],'edgecolor',[1 1 1],'color','k','fontsize',14,'string',['y = ax']);
     annotation('textbox',[0.15 0.80 0.2 0.2],'edgecolor',[1 1 1],'color','r','fontsize',14,'string',['a = ' num2str(kk_D(1),'%1.2f') '\pm' num2str(kk_err_D(1)^0.5,'%1.2f') ' Deuterium']);
     annotation('textbox',[0.15 0.75 0.2 0.2],'edgecolor',[1 1 1],'color','b','fontsize',14,'string',['a = ' num2str(kk_H(1),'%1.2f') '\pm' num2str(kk_err_H(1)^0.5,'%1.2f') ' Hydrogen']);      
     annotation('textbox',[0.15 0.70 0.2 0.2],'edgecolor',[1 1 1],'color','k','fontsize',14,'string',['position ' num2str(position_test*20-140) ' [mm]']);                                                                               
     xlabel('n^{-}_{PD} [10^{17} m^{-3}]','fontsize',14)
     ylabel('n^{-}_{LP} [10^{17} m^{-3}]','fontsize',14)
     title(['LP vs PD Linear Correlation Peak ' num2str(pd_peak)],'fontsize',14)
     grid on
     %legend('cs-Deuterium','cs-Hydrogen','location','northwest')
     set(gca,'fontsize',14)
     saveas(hh, ['lincorrP' num2str(pd_peak) '.png'])
   endif
   
 if(position_test==7)
     kk_D=[0 0];
     kk_H=[0 0];
f=@(b,x) b(1).*x+b(2);
     [kk_D,R,J,CovB,MSE]=nlinfit(A_D,LP_D,f,kk_D,opts);
     kk_err_D=sqrt(diag(CovB));
     [kk_H,R,J,CovB,MSE]=nlinfit(A_H,LP_H,f,kk_H,opts);
     kk_err_H=sqrt(diag(CovB));
     
     f_x=0:0.1:80;
     fvecD=f(kk_D,f_x);
     fvecH=f(kk_H,f_x);
     
     figure('name',[name ' ' insert_name], 'visible','on')
     if(0)
      pmean1=zeros(100,1);pmean2=zeros(100,1);pmean3=zeros(100,1);pmean4=zeros(100,1);
     mean1=zeros(100,1);mean2=zeros(100,1);mean3=zeros(100,1);mean4=zeros(100,1);
     co_p1=0;co_p2=0;co_p3=0;co_p4=0;
     for ie=1:length(PD_D)
       if(p_D(ie)==1)
       co_p1=co_p1+1;
       mean1(co_p1)=LP_D(ie);
       pmean1(co_p1)=A_D(ie);
       elseif(p_D(ie)==2)
       co_p2=co_p2+1;
       mean2(co_p2)=LP_D(ie);
       pmean2(co_p2)=A_D(ie);
       elseif(p_D(ie)==3)
       co_p3=co_p3+1;
       mean3(co_p3)=LP_D(ie);
       pmean3(co_p3)=A_D(ie);
       elseif(p_D(ie)==4)
       co_p4=co_p4+1;
       mean4(co_p4)=LP_D(ie);
       pmean4(co_p4)=A_D(ie);
       endif
     endfor
     hh=errorbar(mean(pmean1(1:co_p1)),mean(mean1(1:co_p1)),std(mean1(1:co_p1)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean2(1:co_p2)),mean(mean2(1:co_p2)),std(mean2(1:co_p2)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean3(1:co_p3)),mean(mean3(1:co_p3)),std(mean3(1:co_p3)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean4(1:co_p4)),mean(mean4(1:co_p4)),std(mean4(1:co_p4)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
          pmean1=zeros(100,1);pmean2=zeros(100,1);pmean3=zeros(100,1);pmean4=zeros(100,1);
     mean1=zeros(100,1);mean2=zeros(100,1);mean3=zeros(100,1);mean4=zeros(100,1);
     co_p1=0;co_p2=0;co_p3=0;co_p4=0;
     for ie=1:length(PD_H)
       if(p_H(ie)==1)
       co_p1=co_p1+1;
       mean1(co_p1)=LP_H(ie);
       pmean1(co_p1)=A_H(ie);
       elseif(p_H(ie)==2)
       co_p2=co_p2+1;
       mean2(co_p2)=LP_H(ie);
       pmean2(co_p2)=A_H(ie);
       elseif(p_H(ie)==3)
       co_p3=co_p3+1;
       mean3(co_p3)=LP_H(ie);
       pmean3(co_p3)=A_H(ie);
       elseif(p_H(ie)==4)
       co_p4=co_p4+1;
       mean4(co_p4)=LP_H(ie);
       pmean4(co_p4)=A_H(ie);
       endif
     endfor
     hh=errorbar(mean(pmean1(1:co_p1)),mean(mean1(1:co_p1)),std(mean1(1:co_p1)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean2(1:co_p2)),mean(mean2(1:co_p2)),std(mean2(1:co_p2)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean3(1:co_p3)),mean(mean3(1:co_p3)),std(mean3(1:co_p3)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean4(1:co_p4)),mean(mean4(1:co_p4)),std(mean4(1:co_p4)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
    end
     hh=plot(A_D,LP_D,'.r')
     hold on
     hh=plot(A_H,LP_H,'.b')
     hold on
     hh=plot(f_x,fvecD,'r','linewidth',1.5)
     hold on
     hh=plot(f_x,fvecH,'b','linewidth',1.5)
     
     xlim([0 80])
     if(pd_peak==3)
      xlim([0 1.5])
     endif
     ylim([0 2])
     
     annotation('textbox',[0.15 0.85 0.2 0.2],'edgecolor',[1 1 1],'color','k','fontsize',14,'string',['y = ax']);
     annotation('textbox',[0.15 0.75 0.2 0.2],'edgecolor',[1 1 1],'color','r','fontsize',14,'string',['a = ' num2str(kk_D(1)*1000,'%1.2f') '\pm' num2str(kk_err_D(1)^0.5*1000,'%1.2f') ; ' b = ' num2str(kk_D(2)*1000,'%1.2f') '\pm' num2str(kk_err_D(2)^0.5*1000,'%1.2f')  ' Deuterium']);
     annotation('textbox',[0.15 0.65 0.2 0.2],'edgecolor',[1 1 1],'color','b','fontsize',14,'string',['a = ' num2str(kk_H(1)*1000,'%1.2f') '\pm' num2str(kk_err_H(1)^0.5*1000,'%1.2f') ; ' b = ' num2str(kk_H(2)*1000,'%1.2f') '\pm' num2str(kk_err_H(2)^0.5*1000,'%1.2f') ' Hydrogen']);      
     annotation('textbox',[0.15 0.60 0.2 0.2],'edgecolor',[1 1 1],'color','k','fontsize',14,'string',['position ' num2str(position_test*20-140) ' [mm]']);                                                                               
     xlabel('Arc Power [kW]','fontsize',14)
     ylabel('n^{-}_{LP} [10^{17} m^{-3}]','fontsize',14)
     title(['LP vs Arc Power'],'fontsize',14)
     grid on
     %legend('cs-Deuterium','cs-Hydrogen','location','northwest')
     set(gca,'fontsize',14)
     saveas(hh, ['lincorrArc' num2str(pd_peak) '.png'])
     endif
%-------------------------------------------------------------------------------
if(position_test==7)
     kk_D=[0 0];
     kk_H=[0 0];
f=@(b,x) b(1).*x+b(2);
     [kk_D,R,J,CovB,MSE]=nlinfit(A_D,PD_D,f,kk_D,opts);
     kk_err_D=sqrt(diag(CovB));
     [kk_H,R,J,CovB,MSE]=nlinfit(A_H,PD_H,f,kk_H,opts);
     kk_err_H=sqrt(diag(CovB));
     
     f_x=0:0.1:80;
     fvecD=f(kk_D,f_x);
     fvecH=f(kk_H,f_x);
     
     figure('visible','on')
     if(0)
      pmean1=zeros(100,1);pmean2=zeros(100,1);pmean3=zeros(100,1);pmean4=zeros(100,1);
     mean1=zeros(100,1);mean2=zeros(100,1);mean3=zeros(100,1);mean4=zeros(100,1);
     co_p1=0;co_p2=0;co_p3=0;co_p4=0;
     for ie=1:length(PD_D)
       if(p_D(ie)==1)
       co_p1=co_p1+1;
       mean1(co_p1)=PD_D(ie);
       pmean1(co_p1)=A_D(ie);
       elseif(p_D(ie)==2)
       co_p2=co_p2+1;
       mean2(co_p2)=PD_D(ie);
       pmean2(co_p2)=A_D(ie);
       elseif(p_D(ie)==3)
       co_p3=co_p3+1;
       mean3(co_p3)=PD_D(ie);
       pmean3(co_p3)=A_D(ie);
       elseif(p_D(ie)==4)
       co_p4=co_p4+1;
       mean4(co_p4)=PD_D(ie);
       pmean4(co_p4)=A_D(ie);
       endif
     endfor
     hh=errorbar(mean(pmean1(1:co_p1)),mean(mean1(1:co_p1)),std(mean1(1:co_p1)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean2(1:co_p2)),mean(mean2(1:co_p2)),std(mean2(1:co_p2)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean3(1:co_p3)),mean(mean3(1:co_p3)),std(mean3(1:co_p3)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
     hh=errorbar(mean(pmean4(1:co_p4)),mean(mean4(1:co_p4)),std(mean4(1:co_p4)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','r','markerfacecolor','r')
     hold on
          pmean1=zeros(100,1);pmean2=zeros(100,1);pmean3=zeros(100,1);pmean4=zeros(100,1);
     mean1=zeros(100,1);mean2=zeros(100,1);mean3=zeros(100,1);mean4=zeros(100,1);
     co_p1=0;co_p2=0;co_p3=0;co_p4=0;
     for ie=1:length(PD_H)
       if(p_H(ie)==1)
       co_p1=co_p1+1;
       mean1(co_p1)=PD_H(ie);
       pmean1(co_p1)=A_H(ie);
       elseif(p_H(ie)==2)
       co_p2=co_p2+1;
       mean2(co_p2)=PD_H(ie);
       pmean2(co_p2)=A_H(ie);
       elseif(p_H(ie)==3)
       co_p3=co_p3+1;
       mean3(co_p3)=PD_H(ie);
       pmean3(co_p3)=A_H(ie);
       elseif(p_H(ie)==4)
       co_p4=co_p4+1;
       mean4(co_p4)=PD_H(ie);
       pmean4(co_p4)=A_H(ie);
       endif
     endfor
     hh=errorbar(mean(pmean1(1:co_p1)),mean(mean1(1:co_p1)),std(mean1(1:co_p1)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean2(1:co_p2)),mean(mean2(1:co_p2)),std(mean2(1:co_p2)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean3(1:co_p3)),mean(mean3(1:co_p3)),std(mean3(1:co_p3)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
     hold on
     hh=errorbar(mean(pmean4(1:co_p4)),mean(mean4(1:co_p4)),std(mean4(1:co_p4)))
     set(hh,'linestyle','none')
     set(hh,'marker','o','markersize',5)
     set(hh,'color','b','markerfacecolor','b')
    end
     hold on
     hh=plot(A_D,PD_D,'.r')
     hold on
     hh=plot(A_H,PD_H,'.b')
     hold on
     hh=plot(f_x,fvecD,'r','linewidth',1.5)
     hold on
     hh=plot(f_x,fvecH,'b','linewidth',1.5)
     
     %xlim([0 80])
     %if(pd_peak==3)
     % xlim([0 1.5])
     %endif
     %yticks([0 30 60 90 120 150])
     ylim([0 3])
     
     annotation('textbox',[0.15 0.45 0.2 0.2],'edgecolor',[1 1 1],'color','k','fontsize',14,'string',['y = ax+b']);
     annotation('textbox',[0.15 0.35 0.2 0.2],'edgecolor',[1 1 1],'color','r','fontsize',14,'string',['a = ' num2str(kk_D(1)*1000,'%1.2f') '\pm' num2str(kk_err_D(1)^0.5*1000,'%1.2f') ; ' b = ' num2str(kk_D(2),'%1.2f') '\pm' num2str(kk_err_D(2)^0.5,'%1.2f')  ' Deuterium']);
     annotation('textbox',[0.15 0.25 0.2 0.2],'edgecolor',[1 1 1],'color','b','fontsize',14,'string',['a = ' num2str(kk_H(1)*1000,'%1.2f') '\pm' num2str(kk_err_H(1)^0.5*1000,'%1.2f') ; ' b = ' num2str(kk_H(2),'%1.2f') '\pm' num2str(kk_err_H(2)^0.5,'%1.2f') ' Hydrogen']);      
     annotation('textbox',[0.15 0.20 0.2 0.2],'edgecolor',[1 1 1],'color','k','fontsize',14,'string',['position ' num2str(position_test*20-140) ' [mm]']);                                                                               
     xlabel('Arc Power [kW]','fontsize',14)
     ylabel('n^{-}_{PD} [10^{17} m^{-3}]','fontsize',14)
     title(['PD vs Arc Power'],'fontsize',14)
     grid on
     %legend('cs-Deuterium','cs-Hydrogen','location','northwest')
     set(gca,'fontsize',14)
     saveas(hh, ['lincorrArc' num2str(pd_peak) '.png'])
     endif
   endfor
%----------------------------------------------------------------------------------------------------------------------------------
   figure
   
   hh=errorbar((1:12)*20-140,ksum_D(:,1),ksum_D(:,2));
   set(hh,'linestyle','none')
   set(hh,'marker','*')
   set(hh,'color','r')

   hold on
   
   hh=errorbar((1:12)*20-140,ksum_H(:,1),ksum_H(:,2));
   set(hh,'linestyle','none')
   set(hh,'marker','*')
   set(hh,'color','b')
   
   xlabel('Position [mm]','fontsize',14)
   ylabel('n^{-}_{LP}/n^{-}_{PD}','fontsize',14)
   title(['LP vs PD Correlation Factor Profile Peak ' num2str(pd_peak)],'fontsize',14)
   grid on
   legend('cs-Deuterium','cs-Hydrogen','location','northwest')
   set(gca,'fontsize',14)
   
   ylim([0 2.5])
   %if(pd_peak==3)
   % ylim([0 20])
   %endif
   xlim([-150 150])
   saveas(hh,['corrfacP' num2str(pd_peak) '.png'])