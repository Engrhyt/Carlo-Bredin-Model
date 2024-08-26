disp(['fitting_' num2str(PD_index)])
    k1=40000*(PD_index-1)+1;
    k2=40000*PD_index;
    k3=k1+11000;
    k4=k1+17000;
    clear('G')
    G=-PD_pre(k3+500:k3+2500,1);
    G=(G)-mean(G(1:200));
    Gt=(1:length(G));
    Gx=(1:length(G));
    
    %starting parameter
    A0=max(G(535:555));A1=0.5*A0;A2=1*A0;
    ub1=[555     10     A2];
    lb1=[535     0.006  A1];
    b1 =[550     1      A0]; 
    H1=@(b,x) M(b(1:3),x,ub(1:3),lb(1:3));
    opts = optimset ("MaxIter",1000,'TolFun',1e-12,'TolX',1e-12,"FinDiffType","central");
    try
     [b1]=lsqcurvefit(H1,b1,(1:555),G(1:555)',lb1,ub1,opts);
    catch
     try 
      [b1]=nlinfit((1:555),G(1:555)',H1,b1,opts);
      if(b1(3)>A0*3)
       b1(3)=A0;
      end
     catch
      b1(3)=A0;
     end
    end
  
    B0=max(G(560:590));B1=0.5*B0;B2=1.5*B0;
    C0=max(G(590:1100));C1=0.01*C0;C2=1*C0;
    ub = [590    200 500  B2   1100   500   500  C2];
    lb = [560      0   0  B1    590     0     0  C1];
    b  = [564     30  20  B0    650   180    50  C0];
    HH = @(b,x) W1(b(1:4),x,ub(1:4),lb(1:4)) + W1(b(5:8),x,ub(5:8),lb(5:8));
    GH=G'-M(b1(1:3),1:length(G),ub1(1:3),lb1(1:3));

    try
     [b]=lsqcurvefit(HH,b,(1:length(G)),GH,lb,ub);
    catch
     try 
       [b]=nlinfit((1:length(G)),GH,HH,b,opts);
       if(b(4)>B0*3)
        b(4)=B0;
       end
       if(b(8)>C0*3)
        b(8)=C0;
       end
     catch
       b(4)=B0;
       b(8)=C0;
     end
    end
    
    figure('visible','off')
    ff=plot((Gt-1),G,'.k');
 
    hold on
    EH=HH(b,1:0.2:length(G));
    ff=plot((0:0.2:length(G)-1),EH,'k');
    
    hold on
    EH=M(b1(1:3),1:0.2:length(G),ub1(1:3),lb1(1:3));
    ff=plot((0:0.2:length(G)-1),EH,'b');
    
    hold on
    EH=W1(b(1:4),1:0.2:length(G),ub(1:4),lb(1:4));
    ff=plot((0:0.2:length(G)-1),EH,'r');
    
    hold on
    EH=W1(b(5:8),1:0.2:length(G),ub(5:8),lb(5:8));
    ff=plot((0:0.2:length(G)-1),EH,'g');
    ylim([0 3])
    xlim([0 1500])
    xlabel('time [ns]','fontsize',15)
    ylabel('PD amplitude','fontsize',15)
    set(gca,'fontsize',15)
    mkdir([save_analyze_file_place '/T' num2str(table_index,'%1.0f') '/' num2str(call)])
    saveas(ff,[save_analyze_file_place '/T' num2str(table_index,'%1.0f') '/' num2str(call) '/' num2str(call) '_' num2str(PD_index) '.png'])
    close
    
    PD1_raw=b1(3);
    PD2_raw=b(4);
    PD3_raw=b(8);