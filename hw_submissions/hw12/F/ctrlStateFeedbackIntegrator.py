import numpy as np
import control as cnt
import VTOLParam as P

class ctrlStateFeedbackIntegrator:
    def __init__(self):
        # Tuning parameters
        tr_h = 3.0
        zeta_h = 0.707
        wn_h = 2.2 / tr_h

        tr_th = 0.1
        zeta_th = 0.707
        wn_th = 2.2 / tr_th

        M = 10
        tr_z = 10*tr_th
        zeta_z = 0.707
        wn_z = 2.2 / tr_z

        j = (P.mc + 2*P.mr)

        pI_lat = np.array([-2.0])
        pI_long = np.array([-1.7])

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
        
        Cr_lat = np.array([[1.0, 0.0, 0.0, 0.0]])

        A1_lat = np.array([[0.0,0.0,1.0,0.0,0.0],
                           [0.0,0.0,0.0,1.0,0.0],
                           [0.0,-P.Fe/j,-P.mu/j,0.0,0.0],
                           [0.0,0.0,0.0,0.0,0.0],
                           [-1.0, 0.0, 0.0, 0.0,0.0]])
        
        B1_lat = np.array([[0.0],
                           [0.0],
                           [0.0],
                           [1/(P.Jc+2*P.mr*P.d**2)],
                           [0.0]
                      ])
        
        # print(f"A1_lat: {A1_lat}")
        # print(f"A1_lat: {B1_lat}")

        
        # Longitude matrices

        A_long = np.array([[0.0,1.0],
                           [0.0,0.0]])
        B_long = np.array([[0.0],[1/j]])
        C_long = np.array([[1.0,0.0]])
        Cr_long = np.array([[1.0, 0.0, 0.0, 0.0]])

        A1_long = np.array([[0.0,1.0,0.0],
                            [0.0,0.0,0.0],
                            [-1.0,0.0,0.0]])
        
        B1_long = np.array([[0.0],
                            [1/j],
                            [0.0]])
        
        # print(f"A1_long: {A1_long}")
        # print(f"B1_long: {B1_long}")


        des_char_poly_long = np.convolve([1, 2 * zeta_h * wn_h,wn_h**2],np.poly(pI_long))
        des_poles_long = np.roots(des_char_poly_long)

        
        des_char_poly_lat = np.convolve(np.convolve([1, 2 * zeta_th * wn_th, wn_th**2],
                                    [1, 2 * zeta_z * wn_z, wn_z**2]),np.poly(pI_lat))
                                    
        des_poles_lat = np.roots(des_char_poly_lat)

        # Compute the gains if the system is controllable for lat control
        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != 5:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1_lat,B1_lat,des_poles_lat)
            self.K_lat = K1[0][0:4]
            self.ki_lat = K1[0][4]
        # print gains to terminal
        print('K_lat: ', self.K_lat)
        print('ki_lat: ', self.ki_lat)

        # Compute the gains of system if controllable for long control
        if np.linalg.matrix_rank(cnt.ctrb(A1_long, B1_long)) != 3:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1_long, B1_long, des_poles_long)
            self.K_long = K1[0][0:2]
            self.ki_long = K1[0][2]
        # print gains to terminal
        print('K_long: ', self.K_long)
        print('ki_long: ', self.ki_long)

        # Variables for integration
        self.integrator_z = 0.0
        self.integrator_h = 0.0
        self.error_z_d1 = 0.0
        self.error_h_d1 = 0.0

    def update(self,h_r,z_r,x):
        z = x[0][0]
        h = x[1][0]
        theta = x[2][0]
        zdot = x[3][0]
        hdot = x[4][0]
        thetadot = x[5][0]

        array_long = np.array([[h],[hdot]])
        array_lat = np.array([[z],[theta],[zdot],[thetadot]])
        # Make the errors
        error_z = z_r - z
        error_h = h_r - h
        # Integrate the errors
        self.integrator_h = self.integrator_h + (P.Ts / 2.0) * (error_h + self.error_h_d1)
        self.integrator_z = self.integrator_z + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_h_d1 = error_h
        self.error_z_d1 = error_z

        # Update the vertical control
        F_tilde = -self.K_long @ array_long - self.ki_long * self.integrator_h 
        F_unsat = P.Fe + F_tilde
        F = saturate(F_unsat[0], 2*P.fmax)

        # Update the horizontal control
        tau_unsat = -self.K_lat @ array_lat - self.ki_lat * self.integrator_z
        tau = saturate(tau_unsat[0], 2*P.fmax*P.d)
        print(tau_unsat)
        print(tau)
        motor_thrusts = P.mixing @ np.array([[F],[tau]])
        return motor_thrusts


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u