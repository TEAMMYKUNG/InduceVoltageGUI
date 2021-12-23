
# Define (x,y) phase A,B,C
Xa = 0; # phase A
Ya = 30
Xb = 20  # phase B
Yb = 15
Xc = 0  # phase C
Yc = 15
Xp = 0  # กำหนดจุด P
Yp = 0
r = 0.25  # radius of cable or GMR
# กำหนด V
r_a = 500
theta_a = 0
r_b = 500
theta_b = -120
r_c = 500
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
    r,theta = np.abs(z), np.angle(z)
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
(Va,Ia)=cart2pol(np.real(q_cart[0]), np.imag(q_cart[0]))
(Vb,Ib)=cart2pol(np.real(q_cart[1]), np.imag(q_cart[1]))
(Vc,Ic)=cart2pol(np.real(q_cart[2]), np.imag(q_cart[2]))
print("ระยะห่างระหว่างจุด")
print(f"dab={dab:.3f} , dac={dac:.3f} , dbc={dbc:.3f} , dab_p={dab_p:.3f} , dac_p={dac_p:.3f} , dbc_p={dbc_p:.3f}")
print("ระยะห่างระหว่างจุด")
np.set_printoptions(precision=3)
print(f"Va={Va}∠{Ia},Vb={Vb}∠{Ib},Vc={Vc}∠{Ic}")