/**
 * \file reference.h
 * \author Randy Beard <beard@byu.edu>
 *
 * class to manage references
 */

#ifndef REFERENCE_H
#define REFERENCE_H

#include <math.h>
#include <Prandom.h>


// reference structure the reference signals for psi and theta
struct Reference {
  float theta = 0.0;
  float psi = 0.0;
  float phi = 0.0;  
};


class ReferenceUtilities  {
  public:
    float amplitude;
    float frequency;
    float y_offset;

    ReferenceUtilities() {
    }

    void init(float _amplitude=1.0, float _frequency=0.0, float _y_offset=0.0) {
      amplitude = _amplitude;
      frequency = _frequency;
      y_offset = _y_offset;
    }

    float square_signal(float t) {
      if (fmod(t, 1/frequency)<=0.5/frequency) {
        return amplitude + y_offset;
      }
      else {
        return -amplitude + y_offset;
      }
    }

    float sawtooth_signal(float t) {
      float tmp = fmod(t, (0.5/frequency));
      return 4 * amplitude * frequency * tmp - amplitude + y_offset;
    }

    float step_signal(float t) {
      if (t >= 0.0) {
        return amplitude + y_offset;
      }
      else {
        return y_offset;
      }
    }

    // float random_signal(float time) {
    //   //return normalvariate(y_offset, amplitude);
    //   return gauss(y_offset, amplitude);
    // }

    float sin_signal(float  t) {
      float pi = 3.14159;
      return amplitude * sin(2*pi*frequency*t) + y_offset;
    }

    float joystick_pitch() {
      // get reference angle for pitch from joystick
      // when the stick is all up, return 30 degrees, 
      // when the stick is all down, return -30 degrees
      float theta_ref = ((float) analogRead(JOYSTICK_UPDOWN) ) / 1024.0;
        // theta_ref in [0, 1]
      theta_ref = 60 * (PI/180) * (theta_ref - 0.5);
        // theta_ref in [-30, 30] * PI/180
      return theta_ref;
    }

    float joystick_yaw() {
      // get reference angle for yaw from joystick
      // when the stick is full right, return 120 degrees, 
      // when the stick is full left, return -120 degrees
      float psi_ref = ((float) analogRead(JOYSTICK_SIDESIDE) ) / 1024.0;
        // psi_ref in [0, 1]
      psi_ref = 240.0 * (PI/180) * (psi_ref - 0.5);
        // psi_ref in [-120, 120] * PI/180
      return psi_ref;
    }    
};

#endif 
