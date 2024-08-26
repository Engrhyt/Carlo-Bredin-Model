parameter_table = dlmread([workset,'/T',num2str(table_index),'.txt'],',',0,0);
short_index_max=length(parameter_table(:,1));
%-----Plasma Condition-----------------------------------------------------------------------------------------
m=1.6605390666*10^-27;
if(parameter_table(3,2)==1)
%| atom | Caesium | Flow Rate | 
condi=["|H_{2} w/o-CS","|FR.:",num2str(parameter_table(3,1)),"|"];
condi2=['HwoCs'];
arc_rawfileplace=[directory,'/rawdata/Praw/HwoCs'];
CRD_rawfileplace=[directory,'/rawdata/CRD/HwoCs'];
PD_rawfileplace=[directory,'/rawdata/PD/HwoCs'];
elseif(parameter_table(3,2)==2)
condi=["|D_{2} w/o-Cs","|FR.:",num2str(parameter_table(3,1)),"|"];
condi2=['DwoCs'];
m=m*2.014;
arc_rawfileplace=[directory,'/rawdata/Praw/DwoCs'];
CRD_rawfileplace=[directory,'/rawdata/CRD/DwoCs'];
PD_rawfileplace=[directory,'/rawdata/PD/DwoCs'];
elseif(parameter_table(3,2)==3)
condi=["|H_{2} w-Cs","|FR.:",num2str(parameter_table(3,1)),"|"];
condi2=['HwCs'];
arc_rawfileplace=[directory,'/rawdata/Praw/HwCs'];
CRD_rawfileplace=[directory,'/rawdata/CRD/HwCs'];
PD_rawfileplace=[directory,'/rawdata/PD/HwCs'];
elseif(parameter_table(3,2)==4)
condi=["|D_{2} w-Cs","|FR.:",num2str(parameter_table(3,1)),"|"];
condi2=['DwCs'];
m=m*2.014;
arc_rawfileplace=[directory,'/rawdata/Praw/DwCs'];
CRD_rawfileplace=[directory,'/rawdata/CRD/DwCs'];
PD_rawfileplace=[directory,'/rawdata/PD/DwCs'];
else
condi=["|Ar","|FR.:",num2str(parameter_table(3,1)),"|"];
condi2=['Ar'];
endif
%----------Probe Scan--------------------------------------------------------------------------------------------
if(parameter_table(1,1)==0)
scan=["P_{arc}[kW]"];
elseif(parameter_table(1,1)==1)
scan=["[mm]z-axis"];
elseif(parameter_table(1,1)==2)
scan=["[mm]x-axis"];
elseif(parameter_table(1,1)==3)
scan=["Celcius"];
else
scan=["error"];
endif 
%memory allocation----------------------------------------------------------------------------------------------
 position=zeros(short_index_max-4,1);
 P_arc=zeros(short_index_max-4,1);
 np=zeros(short_index_max-4,1);
 as=zeros(short_index_max-4,1);
 CRD=zeros(short_index_max-4,1);
 PD1=zeros(short_index_max-4,1);
 PD2=zeros(short_index_max-4,1);
 PD3=zeros(short_index_max-4,1);
 tt=zeros(short_index_max-4,1);
%plasma short
disp(['New table, perform short call'])
long_IV_curve_index=0;
%short_index_start=5;
for short_index=short_index_start:short_index_max
 call=num2str(parameter_table(short_index,2));
 disp('-------------------------------------------------------------')
 disp([condi num2str(parameter_table(short_index,1)) scan])
 disp(['call: ' call])
 disp(['short ' num2str(short_index-4) ' of ' num2str(short_index_max-5)])
%CRD parameter------------------------------------------------------------------
disp('CRD parameter')
CRD_parameter=dlmread([CRD_rawfileplace,'/',call,'.txt'],'	',0,0);
%PDparameter--------------------------------------------------------------------
disp('PD parameter')
call2=num2str(parameter_table(short_index,2)-1);
PD_pre=dlmread([PD_rawfileplace,'/',call2,'.txt'],'	',0,0);
if(table_index>11)
 sampling_rate=40000;
else
 sampling_rate=14000;
endif
n_max=length(PD_pre(:,1))/sampling_rate;
disp(['PD sampling_raterate = ' num2str(n_max)])
if(n_max>=60)
  n_max=60;
else
  n_max;
endif
%arc parameter------------------------------------------------------------------
 disp(['base I_arc'])
 arc_parameter=dlmread([arc_rawfileplace,'/',call,'.txt'],'	',0,0);
 tshift=0;
 %-2.58
 t_start=-3.1
 t_stop=0.1
 start_arc=findPointA(t_start,arc_parameter(:,1)+tshift,[3/5 1],0.001);
 stop_arc=findPointA(t_stop,arc_parameter(:,1)+tshift,[3/5 1],0.001);
 arc_range=(start_arc:stop_arc);
 prerange=1:100;
 I_arc_base=sum(arc_parameter(prerange,17))/length(prerange);
%IV parameter---------------------------------------------------------------
 disp(['Plasma Start-Stop'])
 IV_parameter=dlmread([savefileplace,'/',call,'.txt'],',',0,0);
 end_length=0;
 plasma_start=findPointA(t_start,real(IV_parameter(1:end,9)),[0 1],0.001);
 plasma_stop=findPointA(t_stop,real(IV_parameter(1:end,9)),[0 1],0.001);
 disp(['Plasma Start-Stop found ! : length = ' num2str(length(plasma_start:plasma_stop))])
 %keep memory allocation----------------------------------------------------------------------------------------------
 P_arc_keep=zeros(short_index_max-4,1);
 np_keep=zeros(short_index_max-4,1);
 as_keep=zeros(short_index_max-4,1);
 CRD_keep=zeros(short_index_max-4,1);
 PD1_keep=zeros(short_index_max-4,1);
 PD2_keep=zeros(short_index_max-4,1);
 PD3_keep=zeros(short_index_max-4,1);
 tt_keep=zeros(short_index_max-4,1);
 
  plasma_parameter_max=0;
  P_arc_buffraw=0;
  np_buffraw=0;
  as_buffraw=0;
  CRD_buffraw=0;
  PD1_buffraw=0;
  PD2_buffraw=0;
  PD3_buffraw=0;
  tt_buffraw=0;
  
  %Each IV curve
  arc_count=0;
  t_arc_index=1;
  IV_curve_index=0;
  PD_index=0;
  disp('--------------taking data--------------')
 for plasma_range_index=plasma_start:plasma_stop
   
   %disp('catch arc')
   Parc_buff_index=0;
   I_arc_buff=zeros(2,1);
   V_arc_buff=zeros(2,1);
   catch_arc=true;
   while(catch_arc) 
    if(IV_parameter(plasma_range_index,9)-1/40>arc_parameter(t_arc_index,1)+tshift)
      t_arc_index=t_arc_index+1;
    elseif(IV_parameter(plasma_range_index,9)-1/40 <= arc_parameter(t_arc_index,1)+tshift && IV_parameter(plasma_range_index,9)+1/40 >= arc_parameter(t_arc_index,1)+tshift)
      Parc_buff_index=Parc_buff_index+1;
      I_arc_buff(Parc_buff_index)=arc_parameter(t_arc_index,17);
      V_arc_buff(Parc_buff_index)=arc_parameter(t_arc_index,16);
      t_arc_index=t_arc_index+1;
      arc_count=arc_count+1;
    else
      t_arc_index=t_arc_index-1;
      catch_arc=false;
    endif
   endwhile
   
    P_arc_raw=30*2*(mean(I_arc_buff)-I_arc_base)*mean(V_arc_buff);

    %disp('catch LP')
    np_raw=real(IV_parameter(plasma_range_index,4));
    as_raw=real(IV_parameter(plasma_range_index,5));
    
    %disp('catch CRD')
    for kk=1:length(CRD_parameter(1,:))
      t_CRD=(kk/20)-15.1;
      if((t_CRD-1/40)<IV_parameter(plasma_range_index,9)&&(t_CRD+1/40)>IV_parameter(plasma_range_index,9))
        break
      endif
    endfor
    CRD_raw=CRD_parameter(1,kk)-mean(CRD_parameter(1,335:350));
    
    %disp('catch PD')
    for kk=1:80
     t_PD=(kk/20)-3.28;
     if((t_PD-1/40)<IV_parameter(plasma_range_index,9)&&(t_PD+1/40)>IV_parameter(plasma_range_index,9))
        break
      endif
    endfor
    PD_index=kk;
    if(parameter_table(3,2)<=2)
     n_max=length(PD_pre(:,1))/14000;
     PDCS7
    elseif(parameter_table(3,2)>=3)
     n_max=length(PD_pre(:,1))/40000;
     PDCS7
     %PDA PDA1 PDA2 PD1 PD2 PDM1 PDM2
    endif
    
    tt_raw=real(IV_parameter(plasma_range_index,9));
    
   if(isnan(P_arc_raw)||isinf(P_arc_raw)) continue endif
    
  plasma_parameter_max=plasma_parameter_max+1; 
  
  %disp(['isnan valid ' num2str(plasma_range_index)])
  P_arc_buffraw=P_arc_buffraw+P_arc_raw;
  np_buffraw=np_buffraw+np_raw;
  as_buffraw=as_buffraw+as_raw;
  CRD_buffraw=CRD_buffraw+CRD_raw;
  PD1_buffraw=PD1_buffraw+PD1_raw; 
  PD2_buffraw=PD2_buffraw+PD2_raw;
  PD3_buffraw=PD2_buffraw+PD2_raw;
  %----------------------------------------------------------------------------------------
 %disp('keep parameter')
 IV_curve_index=IV_curve_index+1;
 P_arc_keep(IV_curve_index)=P_arc_raw;
 np_keep(IV_curve_index)=np_raw;
 as_keep(IV_curve_index)=as_raw;
 CRD_keep(IV_curve_index)=CRD_raw;
 PD1_keep(IV_curve_index)=PD1_raw;
 PD2_keep(IV_curve_index)=PD2_raw;
 PD3_keep(IV_curve_index)=PD3_raw;
 tt_keep(IV_curve_index)=tt_raw;
  
 endfor%plasma_range
 
 %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 disp(['total IV curve : ' num2str(length(plasma_start:plasma_stop))])
 disp(['isnan valid IV curve : ' num2str(IV_curve_index)])
 disp(['arc point per IV curve : ' num2str(floor(arc_count/length(plasma_start:plasma_stop)))])
 
 disp('----------------- amplitude-------------------')
 valid_IV_index=0;
 for IV_curve_index=1:length(P_arc_keep)
   
 P_arc_raw=P_arc_keep(IV_curve_index);
 np_raw=np_keep(IV_curve_index);
 as_raw=as_keep(IV_curve_index);
 CRD_raw=CRD_keep(IV_curve_index);
 PD1_raw=PD1_keep(IV_curve_index);
 PD2_raw=PD2_keep(IV_curve_index);
 PD3_raw=PD3_keep(IV_curve_index);
 tt_raw=tt_keep(IV_curve_index);

 
 valid_IV_index=valid_IV_index+1;
 long_IV_curve_index=long_IV_curve_index+1;
 
 position(long_IV_curve_index)=parameter_table(short_index,1);
 P_arc(long_IV_curve_index)=P_arc_raw;
 np(long_IV_curve_index)=np_raw;
 as(long_IV_curve_index)=as_raw;
 CRD(long_IV_curve_index)=CRD_raw;
 PD1(long_IV_curve_index)=PD1_raw;
 PD2(long_IV_curve_index)=PD2_raw;
 PD3(long_IV_curve_index)=PD3_raw;
 tt(long_IV_curve_index)=tt_raw;
   
 endfor
 
 %plot time sequence
 %valid point count
 disp(['total IV curve : ' num2str(length(plasma_start:plasma_stop))])
 disp(['valid IV curve : ' num2str(valid_IV_index)])
 disp(['arc point per IV curve : ' num2str(floor(arc_count/length(plasma_start:plasma_stop)))])
 
 %Plot time sequence prepare----------------------------------------------------------------------------------------------------------
 A1=IV_parameter(:,4);
 t_plasma=IV_parameter(:,9);
 t_arc=arc_parameter(:,1);
 I_arc=arc_parameter(:,17);
 V_arc=arc_parameter(:,16);
 
 range3=1:length(t_plasma)-10;
 figure('visible','off','Position',[1 1 900 900])
 
 %plot time sequence------------------------------------------------------------ 
 [ax, h1, h2] = plotyy (t_plasma(range3),A1(range3),t_arc(:,1)+tshift,(I_arc(:).-I_arc_base).*V_arc(:).*(30*2));
 hold on
 [ ~, h3, h4] = plotyy (t_plasma(plasma_start:plasma_stop),A1(plasma_start:plasma_stop),t_arc(arc_range,1)+tshift,(I_arc(arc_range).-I_arc_base).*V_arc(arc_range).*(30*2));
 hold on
 set(ax,{'ycolor'},{'b';'r'})
 set(ax(1),'fontsize',15)
 set(ax(2),'fontsize',15)
 set([h1],'color','c','linewidth',2)
 set([h2],'color','g','linewidth',1)
 set([h3],'color','b','linewidth',2)
 set([h4],'color','r','linewidth',1)

 xlabel("time [s]",'fontsize',15);
 xlim([-15 14])
 %ylim(ax(1),[0 0.5])
 ylim(ax(2),[0 120])
 %ylim(ax(1),[0 2.5])
 ylabel(ax(1), "np [log]",'color','b','fontsize',15);
 ylabel(ax(2), "P_{arc} [kW]",'color','r','fontsize',15);
 set(gca,'xminortick','on','yminortick','on','linewidth',1)
 mkdir([save_analyze_place,'/time/',condi2])
 saveas(gcf,[save_analyze_place,'/time/',condi2,'/',call,'.png'])
 close 

endfor%short_index
disp('------------------------------------------------------------------------')
disp('End of table, perform save')
disp('------------------------------------------------------------------------')
%Vf_err
csvwrite([save_analyze_file_place,'/T',num2str(table_index,'%1.0f'),'raw.csv'],[P_arc CRD PD1 PD2 PD3 np as position tt])