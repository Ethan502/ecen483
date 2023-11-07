import numpy as np
import control as cnt
import VTOLParam as P

class ctrlStateFeedback:
    def __init__(self):
        # Tuning parameters
        tr_h = 2.5
        zeta_h = 0.707
        wn_h = 2.2 / tr_h

        tr_th = 0.18
        zeta_th = 0.707
        wn_th = 2.2 / tr_th

        M = 10
        tr_z = 10*tr_th
        zeta_z = 0.707
        wn_z = 2.2 / tr_z

        j = (P.mc + 2*P.mr)

        A_lat = np.array([
            [0.0,0.0,1.0,0.0],
            [0.0,0.0,0.0,1.0],
            [0.0,-P.Fe/j,-P.mu/j,0.0],
            [0.0,0.0,0.0,0.0]
                      ])
        
        B_lat = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [1/(P.Jc+2*P.mr*P.d**2)]
                      ])
        
        C_lat = np.array([[1.0,0.0,0.0,0.0],
                      [0.0,1.0,0.0,0.0]])
        
        # Longitude matrices

        A_long = np.array([[0.0,1.0],
                           [0.0,0.0]])
        B_long = np.array([[0.0],[1/j]])
        C_long = np.array([[1.0,0.0]])

        des_char_poly_long = [1, 2 * zeta_h * wn_h,wn_h**2]
        des_poles_long = np.roots(des_char_poly_long)

        
        des_char_poly_lat = np.convolve([1, 2 * zeta_th * wn_th, wn_th**2],
                                    [1, 2 * zeta_z * wn_z, wn_z**2])
                                    
        des_poles_lat = np.roots(des_char_poly_lat)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A_lat, B_lat)) != 4:
            print("The system is not controllable")
        else:
            self.K_lat = cnt.acker(A_lat, B_lat, des_poles_lat)
            Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
            self.kr_lat = -1.0 / (Cr @ np.linalg.inv(A_lat - B_lat @ self.K_lat) @ B_lat)
        # print gains to terminal
        print('K_lat: ', self.K_lat)
        print('kr_lat: ', self.kr_lat)

        if np.linalg.matrix_rank(cnt.ctrb(A_long, B_long)) != 2:
            print("The system is not controllable")
        else:
            self.K_long = cnt.acker(A_long, B_long, des_poles_long)
            Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
            self.kr_long = -1.0 / (C_long @ np.linalg.inv(A_long - B_long @ self.K_long) @ B_long)
        # print gains to terminal
        print('K_long: ', self.K_long)
        print('kr_long: ', self.kr_long)

    def update(self,h_r,z_r,x):
        z = x[0][0]
        h = x[1][0]
        theta = x[2][0]
        zdot = x[3][0]
        hdot = x[4][0]
        thetadot = x[5][0]

        array_long = np.array([[h],[hdot]])
        array_lat = np.array([[z],[theta],[zdot],[thetadot]])

        # Update the vertical control
        F_tilde = -self.K_long @ array_long + self.kr_long * h_r 
        F_unsat = P.Fe + F_tilde[0][0]
        F = saturate(F_unsat, 2*P.fmax)

        # Update the horizontal control
        tau_unsat = -self.K_lat @ array_lat + self.kr_lat * z_r 
        tau = saturate(tau_unsat[0][0], 2*P.fmax*P.d)
        motor_thrusts = P.mixing @ np.array([[F],[tau]])
        return motor_thrusts


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u