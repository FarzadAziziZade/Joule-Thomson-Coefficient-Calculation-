clear all
close all
clc
stop=1;
ii=1;
TT=1;
J=[];
ZZ=[];
R=8.314472;
PC_CH4=4640;
PC_C2H6=4880;
TC_CH4=190.6;
TC_C2H6=305.4;
A_Methane=3.3298e+04;
B_Methane=7.9933e+04;
C_Methane=2.0869e+03;
D_Methane=4.1602e+04;
E_Methane=9.9196e+0;
A_Ethane=4.0326e+04;
B_Ethane=1.3422e+05;
C_Ethane=1.6555e+03;
D_Ethane=7.3223e+04;
E_Ethane=7.5287e+02;
e_CH4=0.011;
e_C2H6=0.099;
Error=0.000001;
Z_OLD=1;
b_CH4=(0.077796074*R*(TC_CH4)/(PC_CH4));
b_C2H6=(0.077796074*R*(TC_C2H6)/(PC_C2H6));

w_CH4=input("Enter molar fraction of Methane : ");
w_C2H6=100-w_CH4;

x=(w_CH4*b_CH4)+(w_C2H6*b_C2H6);
m_CH4=0.37464+(1.54226*e_CH4)-(0.26992*e_CH4^2);
m_C2H6=0.37464+(1.54226*e_C2H6)-(0.26992*e_C2H6^2);
ac_CH4=(0.45723553*(R^2)*(TC_CH4^2)/PC_CH4);
ac_C2H6=(0.45723553*(R^2)*(TC_C2H6^2)/PC_C2H6);
for T=200:50:400
    Z_OLD=1;
    d=((ac_CH4)*((1+m_CH4*(1-sqrt(T/TC_CH4)))^2));
    f=((ac_C2H6)*((1+m_C2H6*(1-sqrt(T/TC_C2H6)))^2));
    Y=(w_CH4*w_CH4*sqrt(d*d))+(w_CH4*w_C2H6*sqrt(d*f))+(w_CH4*w_C2H6*sqrt(d*f))+(w_C2H6*w_C2H6*sqrt(f*f));
   for i=100:1000:70100
       P=i;
       E=((Y*P)/((R*T)^2));
       G=((x*P)/(R*T));
       K=G-1;
       L=E-(2*G)-(3*(G^2));
       N=(G^3)+(G^2)-(G*E);
       Z_OLD=1;
      while stop
        Fz=(Z_OLD^3)+(K*(Z_OLD^2))+(L*Z_OLD)+N;
        FZprim=(3*(Z_OLD^2))+(2*K*Z_OLD)+L;
        Z_NEW=Z_OLD-(Fz/FZprim);
        if (abs(Z_NEW-Z_OLD)>Error)
          Z_OLD=Z_NEW;
        else
            stop=0;
        end   
      end
      stop=1;
      ZZ(ii,TT)=Z_NEW;
      v_act=Z_NEW*R*T/P;
      ad_CH4=(-1)*(d*(m_CH4))/((1+(m_CH4)*(1-sqrt(T/TC_CH4)))*(sqrt(T*TC_CH4)));
      ad_C2H6=(-1)*(f*m_C2H6)/((1+(m_C2H6)*(1-sqrt(T/TC_C2H6)))*(sqrt(T*TC_C2H6)));
      ad_mix=(((w_CH4)^2)*(ad_CH4))+(((w_C2H6)^2)*(ad_C2H6))+(w_CH4)*(w_C2H6)*(((sqrt(d/f))*(ad_C2H6))+((sqrt(f/d))*(ad_CH4)));
      ad2_CH4=(ac_CH4)*(m_CH4)*(sqrt((TC_CH4)/T))*(1+(m_CH4))/(2*(T)*(TC_CH4));
      ad2_C2H6=(ac_C2H6)*(m_C2H6)*(sqrt((TC_C2H6)/T))*(1+(m_C2H6))/(2*(T)*(TC_C2H6));
      ad2_mix=0.5*w_CH4^2*(ad_CH4^2/sqrt(d^2)+2*ad2_CH4-ad2_CH4/sqrt(d^2))+0.5*w_C2H6^2*(ad_C2H6^2/sqrt(f^2)+2*ad2_C2H6-ad2_C2H6/sqrt(f^2))+w_CH4*w_C2H6*(ad_CH4*ad_C2H6/sqrt(d*f)+ad2_CH4*sqrt(f)/sqrt(d)+ad2_C2H6*sqrt(d)/sqrt(f)+0.5*(ad_CH4^2*sqrt(f)/sqrt(d^3)+ad_C2H6^2*sqrt(d)/sqrt(f^3)));
      CV_R=T*(ad2_mix)/(x*sqrt(8.))*(log((Z_NEW+G*(1+sqrt(2.)))/(Z_NEW+G*(1-sqrt(2.)))));
      pT_v=(R/((v_act)-x))-((ad_mix)/((v_act)*((v_act)+x)+(x*((v_act)-x))));
      aT_p=(P*((ad_mix)-(2*Y/T)))/((R*T)^2);
      bT_P=((-1)*x*P)/(R*(T^2));
      ZT_P1=(aT_p)*(G-(Z_NEW))+(bT_P)*(6*G*(Z_NEW)+2*(Z_NEW)-(3*(G^2))-(2*G)+(E)-((Z_NEW)^2));
      ZT_P=(ZT_P1)/((3*((Z_NEW)^2))+(2*(G-1)*(Z_NEW))+(E-(2*G)-(3*(G^2))));
      VT_P=(R/P)*(T*(ZT_P)+Z_NEW);
      CP_R=(CV_R)+(T*(pT_v)*(VT_P))-R;
      CP_CH4=A_Methane+(B_Methane*((C_Methane/T)/sinh(C_Methane/T))^2)+(D_Methane*((E_Methane/T)/cosh(E_Methane/T))^2);
      CP_C2H6=A_Ethane+(B_Ethane*((C_Ethane/T)/sinh(C_Ethane/T))^2)+(D_Ethane*((E_Ethane/T)/cosh(E_Ethane/T))^2);
      C_P=(CP_R)+(((w_CH4)*(CP_CH4))+((w_C2H6)*(CP_C2H6)));
      J(ii,TT)=-((((T*(VT_P)-(v_act))/(C_P)))*10^6);
      ii=ii+1;
   end
   ii=1;
   TT=TT+1;
end
writematrix(J,'J_raw.xlsx');
writematrix(ZZ,'Z_raw.xlsx');
disp('Done!')