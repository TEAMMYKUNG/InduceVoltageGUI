from icecream import ic

# !Define (x,y) phase A,B,C
Xa = 0; # phase A
Ya = 15.1
Xb = 4  # phase B
Yb = 15.1
Xc = 4  # phase C
Yc = 12.6
Xp = -2.3  # กำหนดจุด P
Yp = 15.1
r = 12.825  # radius of cable or GMR
# !Defind Phase A,B,C
r_a = 115
theta_a = 0
r_b = 115
theta_b = -120
r_c = 115
theta_c = 120

# ===============================================
# ============= Define Function ================
# ==============================================
import numpy as np

E0=8.854*pow(10,-12)

def pol2cart(r,theta):
    z = r * np.exp(1j * theta)
    x, y = z.real, z.imag
    return x, y
def cart2pol(x, y):
    z = x + y * 1j
    r,theta = np.abs(z), np.angle(z)*57.2958
    # 1rad = 57.2958 deg 
    return r,theta

# ================================================

Xva,Yva = pol2cart(r_a,np.radians(theta_a))
Xvb,Yvb = pol2cart(r_b,np.radians(theta_b))
Xvc,Yvc = pol2cart(r_c,np.radians(theta_c))
V=np.array([[complex(Xva,Yva)/np.sqrt(3)],[complex(Xvb,Yvb)/np.sqrt(3)],[complex(Xvc,Yvc)/np.sqrt(3)]])
dab = np.sqrt((Xa-Xb)**2+(Ya-Yb)**2)
dac = np.sqrt((Xa-Xc)**2+(Ya-Yc)**2)
dbc = np.sqrt((Xb-Xc)**2+(Yb-Yc)**2)
dab_p = np.sqrt((Xa-Xb)**2+(Ya+Yb)**2)
dac_p = np.sqrt((Xa-Xc)**2+(Ya+Yc)**2)
dbc_p = np.sqrt((Xb-Xc)**2+(Yb+Yc)**2)
Matrix=np.array([[np.log(2*Ya/r) , np.log(dab_p/dab) , np.log(dac_p/dac)],
                [np.log(dab_p/dab) , np.log(2*Yb/r) , np.log(dbc_p/dbc)],
                [np.log(dac_p/dac) , np.log(dbc_p/dbc) , np.log(2*Yc/r)]])
Inv_Matrix = np.linalg.inv(Matrix)
Matrix_multi =np.matmul(Inv_Matrix, V)
q_cart=2*np.pi*E0*np.matmul(Inv_Matrix, V)
(Qa,Qra)=cart2pol(np.real(q_cart[0]), np.imag(q_cart[0]))
(Qb,Qrb)=cart2pol(np.real(q_cart[1]), np.imag(q_cart[1]))
(Qc,Qrc)=cart2pol(np.real(q_cart[2]), np.imag(q_cart[2]))
np.set_printoptions(precision=3)


ic(Qa,Qra,Qb,Qrb,Qc,Qrc)
ic(q_cart)

dpa = np.sqrt((Xp-Xa)**2+(Yp-Ya)**2)
dpb = np.sqrt((Xp-Xb)**2+(Yp-Yb)**2)
dpc = np.sqrt((Xp-Xc)**2+(Yp-Yc)**2)
dpa_p = np.sqrt((Xp-Xa)**2+(Yp+Ya)**2)
dpb_p = np.sqrt((Xp-Xb)**2+(Yp+Yb)**2)
dpc_p = np.sqrt((Xp-Xc)**2+(Yp+Yc)**2)

ic(dpa,dpb,dpc,dpa_p,dpb_p,dpc_p)

Vp=((q_cart[0]*np.log(dpa_p/dpa))+(q_cart[1]*np.log(dpb_p/dpb))+(q_cart[2]*np.log(dpc_p/dpc)))/(2*np.pi*E0)
(VP,VI)=cart2pol(np.real(Vp), np.imag(Vp))
ic(Vp,VP,VI)

#====== MAGNETIC FIELD ======
Dpa = np.sqrt((Xp-Xa)**2+(0-Ya)**2)
Dpb = np.sqrt((Xp-Xb)**2+(0-Yb)**2)
Dpc = np.sqrt((Xp-Xc)**2+(0-Yc)**2)
ic(Dpa,Dpb,Dpc)

I = 1606.539
Ia,Ima = pol2cart(I,np.radians(theta_a))
Ib,Imb = pol2cart(I,np.radians(theta_b))
Ic,Imc = pol2cart(I,np.radians(theta_c))
ic(Ia,Ima,Ib,Imb,Ic,Imc)
Iphase = np.array([[complex(Ia,Ima)],[complex(Ib,Imb)],[complex(Ic,Imc)]])
ic(Iphase)
SuperPosition = np.array([[np.log(Dpa/dpa) ,np.log(Dpb/dpb) ,np.log(Dpc/dpc)]])
ic(SuperPosition)
matrix2 = np.matmul(SuperPosition, Iphase)
ic(matrix2)
Ep_c = 2*(10**-7)*100*np.pi*matrix2
ic(Ep_c)
(Ep,Ei)=cart2pol(np.real(Ep_c), np.imag(Ep_c))
ic(Ep,Ei)