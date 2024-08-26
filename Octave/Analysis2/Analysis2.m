close all
pkg load optim
pkg load control
pkg load signal
%constant
q=1.60217663*10^(-19);
me=9.1093837*10^(-31);
A=2*pi*0.00025*0.01+pi*0.00025^2;
%oml coefficient
ie=2*pi/q/A*(pi*me/8/q)^0.5*10^(-17);
directory_name = ['/home/engrhyt/Desktop/carlofitfast/process3'];
color_code=['*r' ; '*g' ; '*b' ; '*k'];
color_code=['*r' ; '*g' ; '*b' ; '*k'];


%1 P_arc 
%2 CRD 
%3 PD1 
%4 PD2 
%5 PD3 
%6 np 
%7 as 
%8 position
%9 tt
data_type=[1 2 3 4 5 6 7 8 9];
data_unit=['Arc Power [kW]' ; 'CRD [10^{17} m^{-3}]' ; 'PD cummulative area [a.u.]' ; 'PD peak1 area [a.u.]' ...
         ; 'PD peak2 area [a.u.]' ; 'PD peak1 amplitude [a.u.]' ; 'PD peak2 amplitude [a.u.]' ; 'PD peak1 max [a.u.]' ...
         ; 'PD peak2 max [a.u.]' ; 'OML plus [mA]' ; 'OML minus [mA]' ; '/beta' ; '\alpha_{+ion}' ; '\alpha_{-ion,e}'...
         ; '+ Charge Saturation [mA]' ; '- Charge Saturation [mA]' ; 'bias [V]' ; 'Stienbrechel minus [mA]' ; 'Stienbrechel plus [mA]' ; '\beta_{stien}'];
data_unit2=['Arc Power [kW]' ; 'CRD [10^{17} m^{-3}]' ; 'PD cummulative area [10^{17} m^{-3}]' ; 'PD peak1 area [10^{17} m^{-3}]' ...
         ; 'PD peak2 area [10^{17} m^{-3}]' ; 'PD peak1 amplitude [10^{17} m^{-3}]' ; 'PD peak2 amplitude [10^{17} m^{-3}]' ; 'PD peak1 max [10^{17} m^{-3}]' ...
         ; 'PD peak2 max [10^{17} m^{-3}]' ; 'OML plus [mA]' ; 'OML minus [mA]' ; '\beta' ; '\alpha_{+ion}' ; '\alpha_{-ion,e}' ...
         ; '+ Charge Saturation [mA]' ; '- Charge Saturation [mA]' ; 'bias [V]' ; 'Stienbrechel minus [mA]' ; 'Stienbrechel plus [mA]' ; '\beta_{stien}'];
data_name=[ 'Arc' ; 'CRD' ; 'PD cummulative area' ; 'PD peak1 area' ...
         ; 'PD peak2 area' ; 'PD peak1 amplitude' ; 'PD peak2 amplitude' ; 'PD peak1 max' ...
         ; 'PD peak2 max' ; 'OML plus' ; 'OML minus' ; 'beta' ; 'alpha_minus' ; 'alpha_plus' ...
         ; '+ Charge Saturation' ; '- Charge Saturation' ; 'bias' ; 'Stienbrechel minus' ; 'Stienbrechel plus' ; 'beta stien'];
data_scale=[1 1e-17 1 1 1 1 1 1 1];

%-----------------------------collectdata---------------------------------------


if(1)
data_count=zeros(12,4,4);
data=zeros(10,12,length(data_type),4,4);
data_err=zeros(10,12,length(data_type),4,4);
for plasma_type=3:4

 if(plasma_type==1)
 name='HwoCs';
 elseif(plasma_type==2)
 name='DwoCs';
 elseif(plasma_type==4)
 name='HwCs';
 elseif(plasma_type==3)
 name='DwCs';
 endif
 disp(['------------' name '------------'])
 l=4;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

 disp('taking data')
 for table_index=1:l
  rawname=['T' num2str(table_index+plasma_type*4-1)];
  disp(rawname)
  H=dlmread([directory_name '/Traw12/' rawname 'raw.csv'],',',0,0);
  disp(['Power max = ' num2str(max(H(:,1))) ' min = ' num2str(min(H(:,1)))])

  %keep data
  for data_index=1:length(H(:,1))

   enable_in=false;
   if(H(data_index,9)<-3.8&&H(data_index,9)>=-3.9)
     power_index=1;
     enable_in=true;
   elseif(H(data_index,9)<-3.7&&H(data_index,9)>=-3.8)
     power_index=2;
     enable_in=true;
   elseif(H(data_index,9)<-3.1&&H(data_index,9)>=-3.7)
     power_index=3;
     enable_in=true;
   elseif(H(data_index,9)<1&&H(data_index,9)>=-3)
     power_index=4;
     enable_in=true;
   endif

   if(enable_in)
    for position_test=1:12

     if(H(data_index,8)<(position_test*20+5)&&H(data_index,8)>(position_test*20-5)&&enable_in)
      reject_point=false;
      data_count(position_test,power_index,plasma_type)=data_count(position_test,power_index,plasma_type)+1;
      for data_type_index=1:length(data_type)
        data(data_count(position_test,power_index,plasma_type),position_test,data_type_index,power_index,plasma_type)=H(data_index,data_type(data_type_index));
      endfor%data_type_index
      if(reject_point)
       data_count(position_test,power_index,plasma_type)=data_count(position_test,power_index,plasma_type)-1;
      endif

     endif%position check

    endfor%position_test

   endif%enable_in

  endfor%data_index
  %clear('H','H_err')
 endfor%table_index

endfor %plasma_type

endif

if(1)
opts = optimset ("MaxIter",70000,'TolFun',1e-20,'TolX',1e-20,"FinDiffType","central")
min_time=-3;
max_time=-0.1;
lin_yes=true;
pd_peak=1;
CRDfit6
pd_peak = 1
  cal_PDCRD888
  PDdensityprofile2
  Negativeiondensitypower7
  LPprofile
endif