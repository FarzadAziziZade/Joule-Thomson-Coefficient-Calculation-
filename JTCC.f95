!in the name of Allah, the entirely merciful, the easpecially merciful
Program JouleThomson
implicit none
real::R,PC_CH4,PC_C2H6,TC_CH4,TC_C2H6,b_CH4,b_C2H6,w_CH4,w_C2H6,x
real::m_CH4,m_C2H6,ac_CH4,ac_C2H6,e_CH4,e_C2H6
real::d,f,y,E,G,K,L,A_Methane,B_Methane,C_Methane,D_Methane,E_Methane
real::A_Ethane,B_Ethane,C_Ethane,D_Ethane,E_Ethane,CP_CH4,CP_C2H6
real::N,Z_NEW,Z_OLD,Fz,Fzprim,Error,ad_CH4,ad_C2H6,ad_mix,ad2_CH4,ad2_C2H6
real::ad2_mix,CV_R,v_act,pT_v,aT_p,bT_p,zT_P,zT_P1,VT_P,CP_R
real::T,P,C_P,J_T,i
R=8.314472
PC_CH4=4640
PC_C2H6=4880
TC_CH4=190.6
TC_C2H6=305.4
A_Methane=3.3298e+04
B_Methane=7.9933e+04
C_Methane=2.0869e+03
D_Methane=4.1602e+04
E_Methane=9.9196e+0
A_Ethane=4.0326e+04
B_Ethane=1.3422e+05
C_Ethane=1.6555e+03
D_Ethane=7.3223e+04
E_Ethane=7.5287e+02
e_CH4=0.011
e_C2H6=0.099
Error=0.000001
Z_old=1
b_CH4=(0.077796074*R*(TC_CH4)/(PC_CH4))
b_C2H6=(0.077796074*R*(TC_C2H6)/(PC_C2H6))
write(*,*)'Enter molar fraction of Methane'
read(*,*)w_CH4
write(*,*)'Enter molar fraction of Ethane'
read(*,*)w_C2H6
x=(w_CH4*b_CH4)+(w_C2H6*b_C2H6)
m_CH4=0.37464+(1.54226*e_CH4)-(0.26992*e_CH4**2)
m_C2H6=0.37464+(1.54226*e_C2H6)-(0.26992*e_C2H6**2)
ac_CH4=(0.45723553*(R**2)*(TC_CH4**2)/PC_CH4)
ac_C2H6=(0.45723553*(R**2)*(TC_C2H6**2)/PC_C2H6)
Do T=200,400,50
   Z_old=1
   d=((ac_CH4)*((1+m_CH4*(1-sqrt(T/TC_CH4)))**2))
   f=((ac_C2H6)*((1+m_C2H6*(1-sqrt(T/TC_C2H6)))**2))
   Y=(w_CH4*w_CH4*sqrt(d*d))+(w_CH4*w_C2H6*sqrt(d*f))+(w_CH4*w_C2H6*sqrt(d*f))+(w_C2H6*w_C2H6*sqrt(f*f))
   do i=1,70000,100  
      P=i    
      E=((y*P)/((R*T)**2))
      G=((x*P)/(R*T))
      K=G-1
      L=E-(2*G)-(3*(G**2))
      N=(G**3)+(G**2)-(G*E)
      Z_old=1
      do
        Fz=(Z_OLD**3)+(K*(Z_OLD**2))+(L*Z_OLD)+N
        Fzprim=(3*(Z_OLD**2))+(2*K*Z_OLD)+L
        Z_NEW=Z_OLD-(Fz/FZprim)
        if(abs(Z_NEW-Z_OLD)>Error)then
          Z_OLD=Z_NEW
        else
           write(*,*)Z_NEW
           exit
        end if   
      end do
      v_act=z_new*R*T/P
      ad_CH4=(-1)*(d*(m_CH4))/((1+(m_ch4)*(1-sqrt(T/TC_CH4)))*(sqrt(T*TC_CH4)))
      ad_C2H6=(-1)*(f*m_C2H6)/((1+(m_C2H6)*(1-sqrt(T/TC_C2H6)))*(sqrt(T*TC_C2H6)))
      ad_mix=(((w_CH4)**2)*(ad_CH4))+(((w_C2H6)**2)*(ad_C2H6))&
      +(w_CH4)*(w_C2H6)*(((sqrt(d/f))*(ad_C2H6))+((sqrt(f/d))*(ad_CH4)))
      ad2_CH4=(ac_CH4)*(m_CH4)*(sqrt((TC_CH4)/T))*(1+(m_CH4))/(2*(T)*(TC_CH4))
      ad2_C2H6=(ac_C2H6)*(m_C2H6)*(sqrt((TC_C2H6)/T))*(1+(m_C2H6))/(2*(T)*(TC_C2H6))
      ad2_mix=0.5*w_ch4**2*(ad_ch4**2/sqrt(d**2)+2*ad2_ch4-ad2_ch4/sqrt(d**2))&
      +0.5*w_c2h6**2*(ad_c2h6**2/sqrt(f**2)+2*ad2_c2h6-ad2_c2h6/sqrt(f**2))&
      +w_ch4*w_c2h6*(ad_ch4*ad_c2h6/sqrt(d*f)+ad2_ch4*sqrt(f)/sqrt(d)+ad2_c2h6*sqrt(d)/sqrt(f)&
      +0.5*(ad_ch4**2*sqrt(f)/sqrt(d**3)+ad_c2h6**2*sqrt(d)/sqrt(f**3)))
      CV_R=T*(ad2_mix)/(x*sqrt(8.))*(log((z_NEW+G*(1+sqrt(2.)))/(z_NEW+G*(1-sqrt(2.)))))
      pT_v=(R/((v_act)-x))-((ad_mix)/((v_act)*((v_act)+x)+(x*((v_act)-x))))
      aT_p=(P*((ad_mix)-(2*Y/T)))/((R*T)**2)
      bT_P=((-1)*x*P)/(R*(T**2))
      zT_P1=(aT_p)*(G-(Z_NEW))+(bT_P)*(6*G*(Z_NEW)+2*(Z_NEW)-(3*(G**2))-(2*G)+(E)-((Z_NEW)**2))
      ZT_P=(zT_P1)/((3*((Z_NEW)**2))+(2*(G-1)*(Z_NEW))+(E-(2*G)-(3*(G**2))))
      VT_P=(R/P)*(T*(zT_P)+Z_New)
      CP_R=(CV_R)+(T*(pT_v)*(VT_P))-R
      CP_CH4=A_Methane+(B_Methane*((C_Methane/T)/sinh(C_Methane/T))**2)+(D_Methane*((E_Methane/T)/cosh(E_Methane/T))**2)
      CP_C2H6=A_Ethane+(B_Ethane*((C_Ethane/T)/sinh(C_Ethane/T))**2)+(D_Ethane*((E_Ethane/T)/cosh(E_Ethane/T))**2)
      C_P=(CP_R)+(((W_CH4)*(CP_ch4))+((W_C2H6)*(CP_c2h6)))
      J_T=(((T*(VT_P)-(v_act))/(C_P)))*10**6
      write(*,*)J_T
  end do
end do      
            
End program        
        