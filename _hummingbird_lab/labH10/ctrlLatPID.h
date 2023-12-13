
#ifndef CONTROL_H
#define CONTROL_H

#include <math.h>

struct {
  float kp_phi = 0.010164; //0.091476;
  float kd_phi = 0.002059804; //0.020479412;
  float kp_psi = 0.308; // 0.7370847;
  float kd_psi = 0.238; // 0.473744;
  float ki_psi = 0.01;
  float kp_theta = 0.686381121089803;
  float ki_theta = 0.02;
  float kd_theta = 0.2205779330047685;
  float km = 0.365;
} gains;

#include "tuning_utilities.h"

// physical parameters of the system
static struct {
  float m1=0.108862;
  float ell1=0.247;
  float m2=0.4717;
  float ell2=-0.039;
  float m3=.1905;
  float g = 9.81;
  float ellT = 0.29;
  float ell3x=-.007;
  float ell3y=.018;
  float J1x = 0.000189;
  float J1y = 0.001953;
  float J1z = 0.001894;
  float J2x = 0.00231;
  float J2y = 0.003274;
  float J2z = 0.003416;
  float J3x = 0.0002222;
  float J3y = 0.0001956;
  float J3z = 0.000027;
  float d = 0.12;
  float fe = (m1*ell1+m2*ell2)*g/ellT;  
  float force_max = 2.0;
} P;

// reference structure the reference signals for psi and theta
struct Reference {
  float theta = 0.0;
  float psi = 0.0;
  float phi = 0.0;  
};

// Lateral controller for hummingbird
class CtrlLatPID {
  public:
    float phi_d1;
    float phi_dot_d2;
    float phi_dot_d3;
    float psi_d1;
    float psi_dot_d2;
    float psi_dot_d3;    
    float integrator_psi;
    float error_psi_d1; 
    float integrator_theta;
    float theta_d2;
    float theta_d3;
    float theta_dot_d2;
    float theta_dot_d3;
    float error_theta_d1;   
    
    CtrlLatPID() {  
    }

    void init() {
      // persistent variables
      phi_d1 = 0.0;
      phi_dot_d2 = 0.0;
      phi_dot_d3 = 0.0;
      psi_d1 = 0.0;
      psi_dot_d2 = 0.0;
      psi_dot_d3 = 0.0; 
      integrator_psi = 0.0;     
      error_psi_d1 = 0.0;     

      integrator_theta = 0.0;
      theta_d2 = 0.0;
      theta_d3 = 0.0;
      theta_dot_d2 = 0.0;
      theta_dot_d3 = 0.0;
      error_theta_d1 = 0.0; 
    }

    void update(float theta_ref,
                float psi_ref, 
                SensorUtilities &sensors, 
                MotorUtilities &rotors, 
                float Ts) {

      // tune gains
      tuneGains();

      // compute phi_dot/psi_dot (with quadratic prediction)
      float phi = sensors.roll;
      float phi_dot_d1 = (phi - phi_d1) / Ts;
      float phi_dot = 3*phi_dot_d1 - 3*phi_dot_d2 + phi_dot_d3;      
      float psi = sensors.yaw;
      float psi_dot_d1 = (psi - psi_d1) / Ts;
      float psi_dot = 3*psi_dot_d1 - 3*psi_dot_d2 + psi_dot_d3; 

      float theta_d1 = sensors.pitch;
      float theta = 3*theta_d1 - 3*theta_d2 + theta_d3;
      float theta_dot_d1 = (sensors.pitch - theta_d2) / Ts;
      float theta_dot = 3*theta_dot_d1 - 3*theta_dot_d2 + theta_dot_d3;     

      // compute feedback linearized loop
      float force_fl = P.fe*cos(theta);

      //------------------------------------------------------------------
      float error_theta = theta_ref - theta;
      
      // update integrator 
      integrator_theta = integrator_theta + (Ts/2.0) * (error_theta + error_theta_d1);
      
      // pitch control
      float f_tilde = gains.kp_theta * error_theta - gains.kd_theta * theta_dot + gains.ki_theta * integrator_theta;         
      float force = force_fl + f_tilde;
      //float force = saturate(unsat_force, -P.force_max, P.force_max);
      Serial.print("Force_fl:");
      Serial.print(force_fl);
      Serial.print("Theta");
      Serial.println(theta);

      //------------------------------------------------------------------


      // Outer loop
      float error_psi = psi_ref - psi;
      float integrator_psi = integrator_psi + (Ts/2.0) * (error_psi + error_psi_d1);
      float phi_ref = gains.kp_psi * error_psi - gains.kd_psi * psi_dot + integrator_psi * gains.ki_psi;

      // Inner Loop
      float error_phi = phi_ref - phi;
      float torque = gains.kp_phi * error_phi - gains.kd_phi * phi_dot;
      //float torque = 0.0;
      
      
      // convert force and torque to pwm and send to motors
      float left_pwm = (force+torque/P.d)/(1.85*gains.km);
      float right_pwm = (force-torque/P.d)/(2.0*gains.km);
      rotors.update(left_pwm, right_pwm); 

      // update all delayed variables
      phi_d1 = phi;
      phi_dot_d3 = phi_dot_d2;
      phi_dot_d2 = phi_dot_d1;
      psi_d1 = psi;
      psi_dot_d3 = psi_dot_d2;
      psi_dot_d2 = psi_dot_d1;      
      error_psi_d1 = error_psi;

      error_psi_d1 = error_psi;

      theta_d3 = theta_d2;
      theta_d2 = theta_d1;
      theta_dot_d3 = theta_dot_d2;
      theta_dot_d2 = theta_dot_d1;
      error_theta_d1 = error_theta;


      
      // print commanded values
      //Serial.print("Psi_ref:");
      //Serial.print(psi_ref*180/PI);
      //Serial.print(",");
      //Serial.print("Psi:");
      //Serial.println(psi*180/PI);
//      Serial.print(",");      
//      Serial.print("Phi_ref:");
//      Serial.println(phi_ref*180/PI);
//      Serial.print(",");
//      Serial.print("Phi:");
//      Serial.println(phi*180/PI);      
    }

    float saturate(float value, float min_value, float max_value) {
      // Implements the saturation function
      return min(max(min_value, value), max_value);
    }  
};

#endif 
