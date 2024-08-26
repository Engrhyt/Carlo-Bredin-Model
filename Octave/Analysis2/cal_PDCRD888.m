 %close all
count_arc=zeros(4,1);
 cal_facgg=zeros(200,4);
 CRD_plot=zeros(200,4);

  y_data_type_index=3;

 figure('visible','on')%,'position',[10 10 500 400])
 plasma_type_list=[4 3 2 1];
for plasma_i=1:2
  plasma_type=plasma_type_list(plasma_i);

 %------------------------------------------------------------------------------

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
 %------------------------------------------------------------------------------

 %------------------------------------------------------------------------------
 disp(['plasma type = ' num2str(plasma_type)])
 power_index=4
 disp(['power index = ' num2str(power_index)])

  HH=0;
  HH_err=0;

  len1=55;
  len2=75;

  callintegral
  
  if(plasma_type==1||plasma_type==2) 
   HH=HH./180;
   for ijj=1:length(C0)
    if(!isna(C0(ijj))&&!isna(HH(ijj))&&HH(ijj)>0&&HH(ijj)<50)
     count_arc(plasma_type)=count_arc(plasma_type)+1;
     CRD_plot(count_arc(plasma_type),plasma_type)=C0(ijj);
     cal_facgg(count_arc(plasma_type),plasma_type)=HH(ijj);
    end%checkisnan  
   end%%ijj
  end%plasma_type
  
if(plasma_type==3||plasma_type==4)
  HH=HH./180 .*1.259*2/50*1000;
 for ijj=1:length(C0)-4
  if(!isna(C0(ijj))&&!isna(HH(ijj))&&HH(ijj)>2)
   if(pd_peak==1)
    if(plasma_type==3)
     if(C0(ijj)<2.5)
      if(HH(ijj)>92)
       count_arc(plasma_type)=count_arc(plasma_type)+1;
       CRD_plot(count_arc(plasma_type),plasma_type)=C0(ijj);
       cal_facgg(count_arc(plasma_type),plasma_type)=HH(ijj);
      end%HH92
     else%C02.5
      if(HH(ijj)>130)
       count_arc(plasma_type)=count_arc(plasma_type)+1;
       CRD_plot(count_arc(plasma_type),plasma_type)=C0(ijj);
       cal_facgg(count_arc(plasma_type),plasma_type)=HH(ijj);
      end%HH130
     end%C02.5
    else%plasmatype3
    count_arc(plasma_type)=count_arc(plasma_type)+1;
    CRD_plot(count_arc(plasma_type),plasma_type)=C0(ijj);
    cal_facgg(count_arc(plasma_type),plasma_type)=HH(ijj);
   end%pd_peak1
  end%isnan
 end%ijj
end%plasma_type if
  
end%plasma_type  
 
  hh=plot(CRD_plot(1:count_arc(plasma_type),plasma_type)*180,cal_facgg(1:count_arc(plasma_type),plasma_type)*180/1000);
  set(hh,'linestyle','none')
  set(hh,'marker','.')
  if(plasma_type==3)
  set(hh,'color','r')
  elseif(plasma_type==4)
  set(hh,'color','b')
  elseif(plasma_type==1)
  set(hh,'color',[0 0.5 1])
  elseif(plasma_type==2)
  set(hh,'color','m')
  endif
  hold on

 x=0:0.1:3.5;
 if(plasma_type==2||plasma_type==1)
  f=@(b,x) b(1).*(x);
  k=[0];
  [k,R,J,CovB,MSE]=nlinfit(CRD_plot(1:count_arc(plasma_type),plasma_type),cal_facgg(1:count_arc(plasma_type),plasma_type),f,k,opts);
  %cal_facgg(1:count_arc(plasma_type),plasma_type)
  k
  kerr=sqrt(diag(CovB));
  if(plasma_type==1)
   hh=plot(x*180,f(k,x)*180/1000,'linestyle','-','color',[0 0.5 1]);
  elseif(plasma_type==2)
   hh=plot(x*180,f(k,x)*180/1000,'-m');
  endif
  annotation('textbox',[0.65 0.35 0.2 0.2],'edgecolor','none','color','k','fontsize',14,'string'...
  ,['y = a(x)'])
  if(plasma_type==1)
   annotation('textbox',[0.6 0.15 0.2 0.2],'edgecolor','none','color',[0 0.5 1],'fontsize',14,'string'...
   ,['w/o-Cs Hydrogen' ; 'a= ' num2str(k(1),'%1.2f') '\pm'  num2str(kerr(1),'%1.2f') ' x 10^{-3}']);
  elseif(plasma_type==2)
   annotation('textbox',[0.6 0.25 0.2 0.2],'edgecolor','none','color','m','fontsize',14,'string'...
   ,['w/o-Cs Deuterium' ; 'a= ' num2str(k(1),'%1.2f')  '\pm'  num2str(kerr(1),'%1.2f') ' x 10^{-3}']);
  endif
 elseif(plasma_type==3||plasma_type==4)
  f=@(b,x) b(1).*(x);
  k=[0];
  [k,R,J,CovB,MSE]=nlinfit(CRD_plot(1:count_arc(plasma_type),plasma_type),cal_facgg(1:count_arc(plasma_type),plasma_type),f,k,opts);
  %cal_facgg(1:count_arc(plasma_type),plasma_type)
  k
  kerr=sqrt(diag(CovB));
  if(plasma_type==3)
   hh=plot(x*180,f(k,x)*180/1000,'-r');
  elseif(plasma_type==4)
   hh=plot(x*180,f(k,x)*180/1000,'-b');
  end
  annotation('textbox',[0.15 0.85 0.2 0.2],'edgecolor','none','color','k','fontsize',14,'string'...
  ,['y = a(x)'])
  if(plasma_type==3)
   annotation('textbox',[0.15 0.7 0.2 0.2],'edgecolor','none','color','r','fontsize',14,'string'...
   ,['w-Cs Deuterium' ; 'a= ' num2str(k(1),'%1.2f') '\pm'  num2str(kerr(1),'%1.2f') ' x 10^{-3}' ]);
  elseif(plasma_type==4)
   annotation('textbox',[0.15 0.55 0.2 0.2],'edgecolor','none','color','b','fontsize',14,'string'...
   ,['w-Cs Hydrogen' ; 'a= ' num2str(k(1),'%1.2f')  '\pm'  num2str(kerr(1),'%1.2f') ' x 10^{-3}']);
  end
 end

 ylabel(['Line Integral Current'  ' by PD [Am]'],'fontsize',14)
 xlabel(['Line Integral Density'  ' by CRD [10^{17}m^{-2}]'],'fontsize',14)
 title(['K_{pd} Peak ' num2str(pd_peak)],'fontsize',14)
 ylim([0 40])
 xlim([0 600])
 
 if(plasma_type==3)
   k_D=k(1);
   k_D1=0;
 elseif(plasma_type==4)
   k_H=k(1);
   k_H1=0;
 endif
endfor%plasma_type
 grid on
 set(gca,'fontsize',14)
 saveas(hh,['PDCRDP' num2str(pd_peak) '.png'])





