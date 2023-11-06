import numpy as np
import control as cnt
import blockbeamParam as P

class ctrlStateFeedback:
    def __init__(self):
        tr_th = 0.15
        zeta_th = 0.707
        wn_th = 2.2 / tr_th

        M = 10
        tr_z = M * tr_th
        zeta_z = 0.707
        wn_z = 2.2 / tr_z

        z_e = P.length/2.0

        A = np.array([
            [0.0,0.0,1.0,0.0],
            [0.0,0.0,0.0,1.0],
            [0.0,-1*P.g,0.0,0.0],
            [-1*(P.m1*P.g)/(((P.m2*P.length**2)/3)+P.m1*z_e**2),0.0,0.0,0.0]
            ])
        B = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [P.length/(((P.m2*P.length**2)/3)+P.m1*z_e**2)]
                    ])

        C = np.array([[1.0,0.0,0.0,0.0],
                      [0.0,1.0,0.0,0.0]])
        
        des_char_poly = np.convolve(
            [1,2*zeta_th*wn_th,wn_th**2],
            [1,2*zeta_z*wn_z,wn_z**2]
            )
        
        des_poles = np.roots(des_char_poly)

        print(A)
        print(cnt.ctrb(A, B))

        
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print("The system is not controllable")
        else:
            self.K = cnt.acker(A, B, des_poles)
            Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
            self.kr = -1.0 / ((Cr @ np.linalg.inv(A - (B @ self.K))) @ B)

        print('K: ', self.K)
        print('kr: ', self.kr)

    def update(self,z_r,x):
        z = x[0][0]
        F_e = (P.m1*P.g / P.length) * z + P.m2*P.g/2.0
        F_tilde = -self.K @ x + self.kr * z_r
        F_unsat = F_e + F_tilde[0][0]
        F = saturate(F_unsat,P.Fmax)
        return F

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u