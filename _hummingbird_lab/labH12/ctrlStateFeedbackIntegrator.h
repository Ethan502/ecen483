/**
 * class to implement controller
 */

#ifndef CONTROL_H
#define CONTROL_H

#include <math.h>

struct {
  float wn_th = 2.5;
  float zeta_th = 0.707;
  float pi_lon = 4.0;
  float wn_phi = 10.0;
  float zeta_phi = 0.707;
  float wn_psi = 1.6;
  float zeta_psi = 0.707;
  float pi_lat = 1.0;
  float km = 0.3300;
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
  float JT = m1*pow(ell1,2)+m2*pow(ell2,2)+J2z+m3*(pow(ell3x,2)+pow(ell3y,2));
  float fe = (m1*ell1+m2*ell2)*g/ellT;  
  float force_max = 2.0;
  float b_theta = ellT/(m1*pow(ell1,2)+m2*pow(ell2,2)+J1y+J2y);
} P;

// Lateral controller for hummingbird
class CtrlStateFeedbackIntegrator {
  private:
    float theta;
    float theta_d2;
    float theta_d3;
    float theta_dot_d2;
    float theta_dot_d3;
    float integrator_lon;
    float error_lon_d1;
    float phi_d1;
    float phi_d2;
    float phi_d3;
    float phi_dot_d2;
    float phi_dot_d3;
    float psi;
    float psi_d1;
    float psi_d2;
    float psi_d3;
    float psi_dot_d2;
    float psi_dot_d3;    
    float integrator_lat;
    float error_lat_d1;    
  
  public:   
    ctrlStateFeedbackIntegrator() {  
    }

    void init() {
      // persistent variables
      theta_d2 = 0.0;
      theta_d3 = 0.0;
      theta_dot_d2 = 0.0;
      theta_dot_d3 = 0.0;
      integrator_lon = 0.0;     
      error_lon_d1 = 0.0;
      phi_d1 = 0.0;
      phi_d2 = 0.0;
      phi_d3 = 0.0;
      phi_dot_d2 = 0.0;
      phi_dot_d3 = 0.0;
      psi_d1 = 0.0;
      psi_d2 = 0.0;
      psi_d3 = 0.0;
      psi_dot_d2 = 0.0;
      psi_dot_d3 = 0.0; 
      integrator_lat = 0.0;     
      error_lat_d1 = 0.0;      
    }

    void update(Reference const &reference, 
                SensorUtilities const &sensors, 
                MotorUtilities &rotors, 
                float Ts) {

      // tune gains
      tuneGains();
      
      float alpha1_lon = gains.pi_lon+2*gains.zeta_th*gains.wn_th;
      float alpha2_lon = 2*gains.zeta_th*gains.wn_th*gains.pi_lon+pow(gains.wn_th,2);
      float alpha3_lon = pow(gains.wn_th,2)*gains.pi_lon;
      float k_theta = alpha2_lon/P.b_theta;
      float k_thetadot = alpha1_lon/P.b_theta;
      float ki_lon = -alpha3_lon/P.b_theta;
      float alpha1_lat = 2*gains.zeta_phi*gains.wn_phi+2*gains.zeta_psi*gains.wn_psi + gains.pi_lat;
      float alpha2_lat = gains.pi_lat*(2*gains.zeta_phi*gains.wn_phi+2*gains.zeta_psi*gains.wn_psi)
                          + pow(gains.wn_phi,2)+pow(gains.wn_psi,2)
                          + 4*gains.zeta_phi*gains.zeta_psi*gains.wn_phi*gains.wn_psi;
      float alpha3_lat = gains.pi_lat*(pow(gains.wn_phi,2)+pow(gains.wn_psi,2)
                        + 4*gains.zeta_phi*gains.zeta_psi*gains.wn_phi*gains.wn_psi)
                        + 2*gains.zeta_phi*gains.wn_phi*pow(gains.wn_psi,2)
                        + 2*gains.zeta_psi*gains.wn_psi*pow(gains.wn_phi,2);
      float alpha4_lat = gains.pi_lat*(2*gains.zeta_phi*gains.wn_phi*pow(gains.wn_psi,2)
                        + 2*gains.zeta_psi*gains.wn_psi*pow(gains.wn_phi,2))
                        + pow(gains.wn_phi,2)*pow(gains.wn_psi,2);
      float alpha5_lat = gains.pi_lat*pow(gains.wn_phi,2)*(gains.wn_psi,2);
      float b1 = 1/P.J1x;
      float a1 = P.ellT*P.fe/(P.JT+P.J1z);
      float k_phi = alpha2_lat/b1;
      float k_psi = alpha4_lat/(a1*b1);
      float k_phidot = alpha1_lat/b1;
      float k_psidot = alpha3_lat/(a1*b1);
      float ki_lat = -alpha5_lat/(a1*b1);
      
      // compute theta and theta_dot (with quadratic prediction)
      float theta_d1 = sensors.pitch;
      float theta = 3*theta_d1 - 3*theta_d2 + theta_d3;
      float theta_dot_d1 = (sensors.pitch - theta_d2) / Ts;
      float theta_dot = 3*theta_dot_d1 - 3*theta_dot_d2 + theta_dot_d3;
      // compute phi/psi and phi_dot/psi_dot (with quadratic prediction)
      float phi_d1 = sensors.roll;
      float phi = 3*phi_d1 - 3*phi_d2 + phi_d3;
      float phi_dot_d1 = (phi_d1 - phi_d2) / Ts;
      phi_dot_d1 = (phi-phi_d1)/Ts;
      float phi_dot = 3*phi_dot_d1 - 3*phi_dot_d2 + phi_dot_d3;
      float psi = sensors.yaw;
      float psi_dot_d1 = (psi - psi_d1) / Ts;
      float psi_dot = 3*psi_dot_d1 - 3*psi_dot_d2 + psi_dot_d3;      

      theta = sensors.pitch;
      // compute feedback linearized force      
      float force_fl = P.g * (P.m1 * P.ell1 + P.m2 * P.ell2) / P.ellT;
      // compute error
      float error_lon = reference.theta - theta;      
      float error_lat = reference.psi - psi; 
      // update integrator 
      integrator_lon += (Ts / 2) * (error_lon + error_lon_d1);
      integrator_lat += (Ts / 2) * (error_lat + error_lat_d1);
      // longitudinal control
      float force_tilde_unsat = -k_theta * theta - k_thetadot * theta_dot - ki_lon * integrator_lon;
      float force_tilde = saturate(force_tilde_unsat,-P.force_max, P.force_max);
      float force = force_tilde + force_fl;
      // lateral controltheta_dot
      float torque = - k_phi * phi - k_psi * psi - k_phidot * phi_dot - k_psidot * psi_dot - ki_lat * integrator_lat;
      // convert force and torque to pwm and send to motors
      float left_pwm = (force+torque/P.d)/(2.0*gains.km);
      float right_pwm = (force-torque/P.d)/(2.0*gains.km);
      rotors.update(left_pwm, right_pwm); 

      // update all delayed variables
      theta_d3 = theta_d2;
      theta_d2 = theta_d1;
      theta_d1 = theta;
      theta_dot_d3 = theta_dot_d2;
      theta_dot_d2 = theta_dot_d1;
      theta_dot_d1 = theta_dot;
      phi_d3 = phi_d2;
      phi_d2 = phi_d1;
      phi_d1 = phi;
      phi_dot_d3 = phi_dot_d2;
      phi_dot_d2 = phi_dot_d1;
      psi_d3 = psi_d2;
      psi_d2 = psi_d1;
      psi_d1 = psi;
      psi_dot_d3 = psi_dot_d2;
      psi_dot_d2 = psi_dot_d1;      
      error_lon_d1 = error_lon;
      error_lat_d1 = error_lat;
      printToSerial(theta, reference.theta, psi, reference.psi);      
    }

    float saturate(float value, float min_value, float max_value) {
      // Implements the saturation function
      return min(max(min_value, value), max_value);
    }  

    void printToSerial(float theta, float theta_ref, float psi, float psi_ref) {
      // print commanded values
      Serial.print("Theta_ref:");
      Serial.print(theta_ref*180/PI);
      Serial.print(",");
      Serial.print("Theta:");
      Serial.print(theta*180/PI);
      Serial.print(",");
      Serial.print("Psi_ref:");
      Serial.print(psi_ref*180/PI);
      Serial.print(",");
      Serial.print("Psi:");
      Serial.print(psi*180/PI);
      Serial.println(",");
    }
};

#endif 
