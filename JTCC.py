import numpy as np

JJ=[]
ZZ=[]
R=float(8.314472)
PC_CH4=float(4640)
PC_C2H6=float(4880)
TC_CH4=float(190.6)
TC_C2H6=float(305.4)
A_Methane=float(3.3298e+04)
B_Methane=float(7.9933e+04)
C_Methane=float(2.0869e+03)
D_Methane=float(4.1602e+04)
E_Methane=float(9.9196e+0)
A_Ethane=float(4.0326e+04)
B_Ethane=float(1.3422e+05)
C_Ethane=float(1.6555e+03)
D_Ethane=float(7.3223e+04)
E_Ethane=float(7.5287e+02)
e_CH4=float(0.011)
e_C2H6=float(0.099)
Error=float(0.1)
Z_OLD=float(1)
stop=True
b_CH4=float((0.077796074*R*(TC_CH4)/(PC_CH4)))
b_C2H6=float((0.077796074*R*(TC_C2H6)/(PC_C2H6)))

w_CH4=float(input("Enter molar fraction of Methane : "))
w_C2H6=100-w_CH4

x=float((w_CH4*b_CH4)+(w_C2H6*b_C2H6))
m_CH4=float(0.37464+(1.54226*e_CH4)-(0.26992*e_CH4**2))
m_C2H6=float(0.37464+(1.54226*e_C2H6)-(0.26992*e_C2H6**2))
ac_CH4=float((0.45723553*(R**2)*(TC_CH4**2)/PC_CH4))
ac_C2H6=float((0.45723553*(R**2)*(TC_C2H6**2)/PC_C2H6))
JJ_old=np.array(range(1,70000,1000))
for T in range(200,400,50):
    Z_OLD=float(1)
    d=float((ac_CH4)*((1+m_CH4*(1-np.sqrt(T/TC_CH4)))**2))
    f=float((ac_C2H6)*((1+m_C2H6*(1-np.sqrt(T/TC_C2H6)))**2))
    Y=float(w_CH4*w_CH4*np.sqrt(d*d))+(w_CH4*w_C2H6*np.sqrt(d*f))+(w_CH4*w_C2H6*np.sqrt(d*f))+(w_C2H6*w_C2H6*np.sqrt(f*f))
    for i in range(100,70000,1000):
        P=float(i)
        E=float((Y*P)/((R*T)**2))
        G=float((x*P)/(R*T))
        K=float(G-1)
        L=float(E-(2*G)-(3*(G**2)))
        N=float((G**3)+(G**2)-(G*E))
        while stop:
            Fz=float((Z_OLD**3)+(K*(Z_OLD**2))+(L*Z_OLD)+N)
            Fzprim=float((3*(Z_OLD**2))+(2*K*Z_OLD)+L)
            Z_NEW=float(Z_OLD)-float(Fz/Fzprim)
            if (np.abs(Z_NEW-Z_OLD)>=Error):
                Z_OLD=float(Z_NEW)
            else:
                stop=False
                
        Z_OLD=float(1)
        stop=True
        v_act=float(Z_NEW*R*T/P)
        ad_CH4=float((-1)*(d*(m_CH4))/((1+(m_CH4)*(1-np.sqrt(T/TC_CH4)))*(np.sqrt(T*TC_CH4))))
        ad_C2H6=float((-1)*(f*m_C2H6)/((1+(m_C2H6)*(1-np.sqrt(T/TC_C2H6)))*(np.sqrt(T*TC_C2H6))))
        ad_mix=float((((w_CH4)**2)*(ad_CH4))+(((w_C2H6)**2)*(ad_C2H6))+(w_CH4)*(w_C2H6)*(((np.sqrt(d/f))*(ad_C2H6))+((np.sqrt(f/d))*(ad_CH4))))
        ad2_CH4=float((ac_CH4)*(m_CH4)*(np.sqrt((TC_CH4)/T))*(1+(m_CH4))/(2*(T)*(TC_CH4)))
        ad2_C2H6=float((ac_C2H6)*(m_C2H6)*(np.sqrt((TC_C2H6)/T))*(1+(m_C2H6))/(2*(T)*(TC_C2H6)))
        ad2_mix=float(0.5*w_CH4**2*(ad_CH4**2/np.sqrt(d**2)+2*ad2_CH4-ad2_CH4/np.sqrt(d**2))+0.5*w_C2H6**2*(ad_C2H6**2/np.sqrt(f**2)+2*ad2_C2H6-ad2_C2H6/np.sqrt(f**2))+w_CH4*w_C2H6*(ad_CH4*ad_C2H6/np.sqrt(d*f)+ad2_CH4*np.sqrt(f)/np.sqrt(d)+ad2_C2H6*np.sqrt(d)/np.sqrt(f)+0.5*(ad_CH4**2*np.sqrt(f)/np.sqrt(d**3)+ad_C2H6**2*np.sqrt(d)/np.sqrt(f**3))))
        CV_R=float(T*(ad2_mix)/(x*np.sqrt(8.))*(np.log(Z_NEW+G*(1+np.sqrt(2.)))-np.log(Z_NEW+G*(1-np.sqrt(2.)))))
        pT_v=float((R/((v_act)-x))-((ad_mix)/((v_act)*((v_act)+x)+(x*((v_act)-x)))))
        aT_p=float((P*((ad_mix)-(2*Y/T)))/((R*T)**2))
        bT_P=float(((-1)*x*P)/(R*(T**2)))
        zT_P1=float((aT_p)*(G-(Z_NEW))+(bT_P)*(6*G*(Z_NEW)+2*(Z_NEW)-(3*(G**2))-(2*G)+(E)-((Z_NEW)**2)))
        ZT_P=float((zT_P1)/((3*((Z_NEW)**2))+(2*(G-1)*(Z_NEW))+(E-(2*G)-(3*(G**2)))))
        VT_P=float((R/P)*(T*(ZT_P)+Z_NEW))
        CP_R=float((CV_R)+(T*(pT_v)*(VT_P))-R)
        CP_CH4=float(A_Methane+(B_Methane*((C_Methane/T)/np.sinh(C_Methane/T))**2)+(D_Methane*((E_Methane/T)/np.cosh(E_Methane/T))**2))
        CP_C2H6=float(A_Ethane+(B_Ethane*((C_Ethane/T)/np.sinh(C_Ethane/T))**2)+(D_Ethane*((E_Ethane/T)/np.cosh(E_Ethane/T))**2))
        C_P=float((CP_R)+(((w_CH4)*(CP_CH4))+((w_C2H6)*(CP_C2H6))))
        J_T=float((((T*(VT_P)-(v_act))/(C_P)))*10**6)
        JJ.append(J_T)
        ZZ.append(Z_NEW)

    np.savetxt("J_at_{0}.csv".format(T), JJ, delimiter=",")
    np.savetxt("Z_at_{0}.csv".format(T), ZZ, delimiter=",")
