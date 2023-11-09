import numpy as np 
import VTOLParam as P


class VTOLDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],  # initial base angle
            [P.h0],  # initial panel angle
            [P.theta0],  # initial angular velocity of base
            [P.zdot0],
            [P.hdot0],
            [P.thetadot0]  # initial angular velocity of panel
        ])
        # simulation time step
        self.Ts = P.Ts
        # mass of engine
        self.mr = P.mr * (1.+alpha*(2.*np.random.rand()-1.))
        # mass of center
        self.mc = P.mc * (1.+alpha*(2.*np.random.rand()-1.))
        # distance from center to engine
        self.d = P.d * (1.+alpha*(2.*np.random.rand()-1.))
        # mew drag
        self.mu = P.mu * (1.+alpha*(2.*np.random.rand()-1.))  
        # gravity
        self.g = P.g
        # intertia of the center
        self.Jc = P.Jc * (1.+alpha*(2.*np.random.rand()-1.))

    def update(self, u,d = 0.0):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input torque
        # u = saturate(u, self.torque_limit)
        self.rk4_step(u,d)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y

    def f(self, state, u, d = 0.0):
        # Return xdot = f(x,u)
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        zdot = state[3][0]
        hdot = state[4][0]
        thetadot = state[5][0]

        fr = u[0][0]
        fl = u[1][0]
        # The equations of motion.
        M = np.array([[self.mc + 2*self.mr,0,0],
                      [0,self.mc + 2*self.mr,0],
                      [0,0,self.Jc + 2*self.mr*self.d**2]])
        C = np.array([[-1*(fr+fl)*np.sin(theta) - self.mu * zdot + d],
                      [-1*(self.mc + 2*self.mr)*self.g + (fr+fl)*np.cos(theta)],
                      [self.d*(fr - fl)]])
        tmp = np.linalg.inv(M) @ C
        zddot = tmp[0][0]
        hddot = tmp[1][0]
        thetaddot = tmp[2][0]

        # build xdot and return
        xdot = np.array([[zdot], [hdot], [thetadot], [zddot], [hddot], [thetaddot]])
        return xdot

    def h(self):
        # return y = h(x)
        theta = self.state[0][0]
        phi = self.state[1][0]
        y = np.array([[theta], [phi]])
        return y

    def rk4_step(self, u, d=0.0):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u, d)
        F2 = self.f(self.state + self.Ts / 2 * F1, u, d)
        F3 = self.f(self.state + self.Ts / 2 * F2, u, d)
        F4 = self.f(self.state + self.Ts * F3, u, d)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

        
def saturate(u, limit):
    if abs(u) > limit:
        u = limit*np.sign(u)
    return u
