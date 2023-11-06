import numpy as np
import control as cnt
import massParam as P

class ctrlStateFeedback:
    def __init__(self):
        tr = 2
        zeta = 0.707
        w_n = 2.2/tr

        A = np.array([[0.0,1.0],
                      [-P.k/P.m,-1.0*P.b/P.m]])
        B = np.array([[0.0],
                      [1.0/P.m]])
        C = np.array([[1.0,0.0]])

        # Do the gain calculation
        des_char_poly = [1,2*zeta*w_n,w_n**2]
        des_poles = np.roots(des_char_poly)
        # Find the gains, but only if the system in controllable
        if np.linalg.matrix_rank(cnt.ctrb(A,B)) != 2:
            print("The system is not controllable")
        else:
            self.K = (cnt.acker(A,B,des_poles))
            self.kr = -1.0 / (C @ np.linalg.inv(A - B @ self.K) @ B)
        print('K: ', self.K)
        print('kr: ', self.kr)
        print('Desired poles: ', des_poles)

    def update(self,z_r,x):
        z = x[0][0]
        # Find the F_fl with feedback linearization
        F_fl = P.k * z
        # Compute the state feedback controller
        F_tilde = -self.K @ x + self.kr * z_r
        # Find the total force
        F_unsat = F_fl + F_tilde[0][0]
        # Saturate and return
        F = saturate(F_unsat,P.F_max)
        return F

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u