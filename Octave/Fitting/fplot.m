u=time;
if(m==1)
rs_mode=1;%rss=rs;
ff=@(b,x) f([b(1) b(2) b(2) b(3) b(4) b(5) rss Spap B rs_mode mr],x);
elseif(m==2)
rs_mode=2;%rss=1;
ff=@(b,x) f([b(1) b(2) b(2) b(3) b(4) b(5) rss Spap B rs_mode mr],x);
elseif(m==3)
rs_mode=3;%rss=fit;
ff=@(b,x) f([b(1) b(2) b(2) b(3) b(4) b(5) b(6) Spap B rs_mode mr],x);
elseif(m==4)
rs_mode=4;%approx OML;
ff=@(b,x) f([b(1) b(2) b(2) b(3) b(4) b(5) 1 Spap B rs_mode mr],x);
elseif(m==5)
rs_mode=1;%rss=rs;
ff=@(b,x) f([Te Ti Ti b(1) b(2) b(3) rss Spap B rs_mode mr],x);
elseif(m==6)
rs_mode=2;%rss=1; 
ff=@(b,x) f([Te Ti Ti b(1) b(2) b(3) rss Spap B rs_mode mr],x);
elseif(m==7)
rs_mode=3;%rss=fit;
ff=@(b,x) f([Te Ti Ti b(1) b(2) b(3) b(4) Spap B rs_mode mr],x);
elseif(m==8)
rs_mode=4;%approx OML;
ff=@(b,x) f([Te Ti Ti b(1) b(2) b(3) 1 Spap B rs_mode mr],x);
end
%x=-40:1:40;
%x=x';
%m=3;
%b=[1 0.1 0.1 17 10 5 -3 1e-5 0.1 3 1];
%b=[b(1) b(2) b(2) b(3) b(4) b(5) b(6) Spap B rs_mode mr];
f(b,x);
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
A=dlmread('out.csv',',');
u=time-u;
disp(['time =' num2str(u) '[s]'])
figure('visible','on')
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

figure('visible','on')
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

figure('visible','on')
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