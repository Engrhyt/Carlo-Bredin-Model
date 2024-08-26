pkg load optim
HD=2;
j1=1;
%j2=1to7 npmin=15
Te=1;Ti=2.2;np=17;as=10;rss=-3;
Temax=5;Timax=5;npmax=19;asmax=1e5,vsmax=10;rsmax=-2;
Temin=0.1;Timin=0.01;npmin=14;asmin=0.1;vsmin=-5;rssmin=-5;
if(HD==1)
mr=1.5;
elseif(HD==2)
mr=2*1.5;
end
m=3;

AA=dlmread(['/home/engrhyt/Desktop/carlofitfast/work/A' num2str(j1) '.txt']);
Bin=dlmread(['/home/engrhyt/Desktop/carlofitfast/magnetic/Spap.txt']);

for j2=5:length(AA(:,1))
call=AA(j2,HD+1);

for Bi=1:length(Bin(:,1))
 if(Bin(Bi,3)==AA(j2,1))
  posi=Bin(Bi,3);
  Spap=Bin(Bi,1);
  B=Bin(Bi,2);
  disp(['catch x=' num2str(Bin(Bi,3)) '[mm] Spap=' num2str(Spap) '[m^{2}] B=' num2str(B) '[mT]'])
  break
 end
end
mkdir('pa')
delete(['pa/' num2str(call) '.txt']);
for j3=50:198
in=dlmread(['/home/engrhyt/Desktop/carlofitfast/IVraw/' num2str(call) '/' num2str(call) '_' num2str(j3) '.txt'],',');
y=in(1:end-1,2);
x=in(1:end-1,1);
t=in(end,1);
ivs=findPointA(0,y,[0.25 0.75],1e-6);
vs=x(ivs);
disp([ num2str(j3) ':' num2str(call) '_' num2str(j3) ' x=' num2str(posi) '[mm] t=' num2str(t) '[s]'])
%Te=Tn=Tp=np=as=Vs=rss=Spap=B=
u=time;
if(m==1)
rs_mode=1;%rss=rs;
bu=[Temax,Timax,npmax,asmax,vsmax];
bl=[Temin,Timin,npmin,asmin,vsmin];
b=[Te,Ti,np, as, vs];
ff=@(b,x) f([b(1) b(2) b(2) b(3) b(4) b(5) rss Spap B rs_mode mr],x);
elseif(m==2)
rs_mode=2;%rss=1;
bu=[Temax,Timax,npmax,asmax,vsmax];
bl=[Temin,Timin,npmin,asmin,vsmin];
b=[Te,Ti,np, as, vs];
ff=@(b,x) f([b(1) b(2) b(2) b(3) b(4) b(5) rss Spap B rs_mode mr],x);
elseif(m==3)
rs_mode=3;%rss=fit;
bu=[Temax,Timax,npmax,asmax,vsmax,rsmax];
bl=[Temin,Timin,npmin,asmin,vsmin,rssmin];
b=[Te,Ti,np, as, vs, rss];
ff=@(b,x) f([b(1) b(2) b(2) b(3) b(4) b(5) b(6) Spap B rs_mode mr],x);
elseif(m==4)
rs_mode=4;%approx OML;
bu=[Temax,Timax,npmax,asmax,vsmax];
bl=[Temin,Timin,npmin,asmin,vsmin];
b=[Te,Ti,np, as, vs];
ff=@(b,x) f([b(1) b(2) b(2) b(3) b(4) b(5) 1 Spap B rs_mode mr],x);
elseif(m==5)
rs_mode=1;%rss=rs;
bu=[npmax,asmax,vsmax];
bl=[npmin,asmin,vsmin];
b=[np, as, vs];
ff=@(b,x) f([Te Ti Ti b(1) b(2) b(3) rss Spap B rs_mode mr],x);
elseif(m==6)
rs_mode=2;%rss=1; 
bu=[npmax,asmax,vsmax];
bl=[npmin,asmin,vsmin];
b=[np, as, vs];
ff=@(b,x) f([Te Ti Ti b(1) b(2) b(3) rss Spap B rs_mode mr],x);
elseif(m==7)
rs_mode=3;%rss=fit;
bu=[npmax,asmax,vsmax,rssmax];
bl=[npmin,asmin,vsmin,rssmin];
b=[np, as, vs, rss];
ff=@(b,x) f([Te Ti Ti b(1) b(2) b(3) b(4) Spap B rs_mode mr],x);
elseif(m==8)
rs_mode=4;%approx OML;
bu=[npmax,asmax,vsmax];
bl=[npmin,asmin,vsmin];
b=[np, as, vs];
ff=@(b,x) f([Te Ti Ti b(1) b(2) b(3) 1 Spap B rs_mode mr],x);
end
try
 %b=nlinfit(x,y,ff,b);
 b=lsqcurvefit(ff,b,x,y,bl,bu);
catch
 try
  %b=nlinfit(x,y,ff,b);
  b=lsqcurvefit(ff,b,x,y,bl,bu);
 catch
  try
   %b=nlinfit(x,y,ff,b);
   b=lsqcurvefit(ff,b,x,y,bl,bu);
  catch
   continue
  end
 end
end
A=dlmread('out.csv',',');
u=time-u;
disp(['time =' num2str(u) '[s]'])
if(m==1)
rec=[b(1) b(2) b(2) b(3) b(4) b(5) rss posi t];
st1=['Te=' num2str(b(1)) '[eV]' ' Ti=' num2str(b(2)) '[ev]' ' np=' num2str(b(3)) '[m^{-3}]' ; 'as=' num2str(b(4)) ' Vs=' num2str(b(5)) '[V]' ' rss=rs' '[m]'];
elseif(m==2)
rec=[b(1) b(2) b(2) b(3) b(4) b(5) rss posi t];
st1=['Te=' num2str(b(1)) '[eV]' ' Ti=' num2str(b(2)) '[ev]' ' np=' num2str(b(3)) '[m^{-3}]' ; 'as=' num2str(b(4)) ' Vs=' num2str(b(5)) '[V]' ' rss=1' '[m]'];
elseif(m==3)
rec=[b(1) b(2) b(2) b(3) b(4) b(5) b(6) posi t];
st1=['Te=' num2str(b(1)) '[eV]' ' Ti=' num2str(b(2)) '[ev]' ' np=' num2str(b(3)) '[m^{-3}]' ; 'as=' num2str(b(4)) ' Vs=' num2str(b(5)) '[V]' ' rss=' num2str(10^(b(6))) '[m]'];
elseif(m==4)
rec=[b(1) b(2) b(2) b(3) b(4) b(5) 1 posi t];
st1=['Te=' num2str(b(1)) '[eV]' ' Ti=' num2str(b(2)) '[ev]' ' np=' num2str(b(3)) '[m^{-3}]' ; 'as=' num2str(b(4)) ' Vs=' num2str(b(5)) '[V]' ' approx OML'];
elseif(m==5)
rec=[Te Ti Ti b(1) b(2) b(3) rss posi t];
st1=['Te=' num2str(Te)   '[eV]' ' Ti=' num2str(Ti)   '[eV]'  'np=' num2str(b(1)) '[m^{-3}]' ; 'as=' num2str(b(2)) ' Vs=' num2str(b(3)) '[V]' ' rss=rs[m]'];
elseif(m==6)
rec=[Te Ti Ti b(1) b(2) b(3) rss posi t];
st1=['Te=' num2str(Te)   '[eV]' ' Ti=' num2str(Ti)   '[eV]'  'np=' num2str(b(1)) '[m^{-3}]' ; 'as=' num2str(b(2)) ' Vs=' num2str(b(3)) '[V]' ' rss=1[m]'];
elseif(m==7)
rec=[Te Ti Ti b(1) b(2) b(3) b(4) posi t];
st1=['Te=' num2str(Te)   '[eV]' ' Ti=' num2str(Ti)   '[eV]'  'np=' num2str(b(1)) '[m^{-3}]' ; 'as=' num2str(b(2)) ' Vs=' num2str(b(3)) '[V]' ' rss=' num2str(10^(b(4))) '[m]'];
elseif(m==8)
rec=[Te Ti Ti b(1) b(2) b(3) b(4) posi t];
st1=['Te=' num2str(Te)   '[eV]' ' Ti=' num2str(Ti)   '[eV]'  'np=' num2str(b(1)) '[m^{-3}]' ; 'as=' num2str(b(2)) ' Vs=' num2str(b(3)) '[V]' ' approx OML'];
end
dlmwrite(['pa/' num2str(call) '.txt'],rec,',','-append')
close all
figure('visible','off')
h=plot(x,y.*1000,'k');
hold on
h=plot(A(:,1),A(:,2).*1000,'r');
ylabel('Probe Current [mA]','fontsize',14)
xlabel('Probe Voltage [V]','fontsize',14)
xlim([-40 40])
%ylim([-5 30])
grid on
pos=[0.15 0.85 0.3 0.3];
st=['time = ' num2str(u) ' [s]' ' mode=' num2str(m)];
 ht = annotation ("textbox", pos, "string",st, ...
                  "edgecolor", "none", "linewidth", 3, "color", "k", ...
                  "verticalalignment", "bottom", "fontsize", 12);
pos=[0.15 0.75 0.3 0.3];
 ht = annotation ("textbox", pos, "string",st1, ...
                  "edgecolor", "none", "linewidth", 3, "color", "k", ...
                  "verticalalignment", "bottom", "fontsize", 12);
set(gca,'fontsize',14)
mkdir(['pic/' num2str(call) '/Imode_' num2str(m)])
saveas(h,['pic/' num2str(call) '/Imode_' num2str(m) '/' num2str(j3) '.png'])
close

figure('visible','off')
h=plot(x,y.*1000,'k');
hold on
h=plot(A(:,1),A(:,3).*1000,'g');
hold on
h=plot(A(:,1),A(:,4).*1000,'r');
hold on
h=plot(A(:,1),A(:,5).*1000,'b');
ylabel('Probe Current [mA]','fontsize',14)
xlabel('Probe Voltage [V]','fontsize',14)
xlim([-40 40])
%ylim([-5 30])
grid on
pos=[0.15 0.85 0.3 0.3];
st=['time = ' num2str(u) ' [s]' ' mode=' num2str(m)];
 ht = annotation ("textbox", pos, "string",st, ...
                  "edgecolor", "none", "linewidth", 3, "color", "k", ...
                  "verticalalignment", "bottom", "fontsize", 12);
pos=[0.15 0.75 0.3 0.3];
 ht = annotation ("textbox", pos, "string",st1, ...
                  "edgecolor", "none", "linewidth", 3, "color", "k", ...
                  "verticalalignment", "bottom", "fontsize", 12);
set(gca,'fontsize',14)
mkdir(['pic/' num2str(call) '/IeIpInmode_' num2str(m)])
saveas(h,['pic/' num2str(call) '/IeIpInmode_' num2str(m) '/' num2str(j3) '.png'])
close

figure('visible','off')
h=plot(A(:,1),A(:,6).*1000,'k');
ylabel('Probe Sheath Radius [mm]','fontsize',14)
xlabel('Probe Voltage [V]','fontsize',14)
xlim([-40 40])
%ylim([0 1])
grid on
pos=[0.15 0.85 0.3 0.3];
st=['time = ' num2str(u) ' [s]' ' mode=' num2str(m)];
 ht = annotation ("textbox", pos, "string",st, ...
                  "edgecolor", "none", "linewidth", 3, "color", "k", ...
                  "verticalalignment", "bottom", "fontsize", 12);
pos=[0.15 0.75 0.3 0.3];
 ht = annotation ("textbox", pos, "string",st1, ...
                  "edgecolor", "none", "linewidth", 3, "color", "k", ...
                  "verticalalignment", "bottom", "fontsize", 12);
set(gca,'fontsize',14)
mkdir(['pic/' num2str(call) '/Rmode_' num2str(m)])
saveas(h,['pic/' num2str(call) '/Rmode_' num2str(m) '/' num2str(j3) '.png'])
close

end%j3
end%j2