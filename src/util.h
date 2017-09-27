#ifndef __PP_UTIL__
#define __PP_UTIL__

#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include <cmath>
#include <cassert>

using namespace std;
using vectord = vector<double>;

inline double sqr(double x){return x*x;};

inline double distance(double x1, double y1, double x2, double y2) {
  return sqrt(sqr(x1-x2) + sqr(y1-y2));
}

inline double distance2(pair<double,double> p, pair<double,double> q) {
  return distance(p.first, p.second, q.first, q.second);
}

template <typename T> inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// template<typename T> inline void min(T a, T b) { return (a<b) ? a:b;}
// template<typename T> inline void max(T a, T b) { return (a<b) ? b:a;}

template<typename T> inline void updatemin(T &a, T b) { if(a>b) a=b;}
template<typename T> inline void updatemax(T &a, T b) { if(a<b) a=b;}

inline double eval_normal_cdf(double x){
  return 0.5 + erfc(-x + M_SQRT1_2);
}

struct NormalCDF{
  double mu, sigma;
  NormalCDF(double u, double s): mu(u), sigma(s){}
  inline double operator()(double x) const {
    return eval_normal_cdf((x-mu)/sigma);
  }
};

struct CoordTransformation {
  double px, py, yaw, sin_yaw, cos_yaw;

  CoordTransformation(double x_, double y_, double yaw_): px(x_), py(y_), yaw(yaw_){
    sin_yaw = sin(yaw);
    cos_yaw = cos(yaw);
  }

  // y  /      
  // | /
  // |/ \yaw
  // ----------x

  void operator() (double &x, double &y) {
    // cout<<"Old : ("<<x<<","<<y<<") --> (";
    // translate
    x-=px; 
    y-=py;

    // rotate
    double new_x = x*cos_yaw + y*sin_yaw;
    double new_y = y*cos_yaw - x*sin_yaw;

    x = new_x;
    y = new_y;

    // cout<<new_x<<","<<new_y<<")"<<endl;
  }

  void revert(double &x, double &y) {
    // rotate -yaw
    double old_x = x*cos_yaw - y*sin_yaw;
    double old_y = y*cos_yaw + x*sin_yaw;

    x = old_x;
    y = old_y;

    // translate
    x += px;
    y += py;
  }

};


template<typename T> inline void erase(vector<T> &v, const T &t){
  v.erase(v.remove(v.begin(), v.end(), t), v.end());
}

class JMT {
private:
  vectord solve(const vectord &start, const vectord &end, double T){
    double T2 =  T*T;
    double T3 = T2*T;
    double T4 = T3*T;
    double T5 = T4*T;

    Eigen::Matrix3d A;
    A <<    T3,   T4,    T5,
          3*T2, 4*T3,  5*T4,
          6*T, 12*T2, 20*T3;
      
    Eigen::Vector3d b;     
    b <<  end[0]-(start[0] + start[1]*T + start[2]*T2/2),
          end[1]-(start[1] + start[2]*T),
          end[2]-start[2];
            
    Eigen::Vector3d C = A.inverse()*b;

    return {start[0], start[1], start[2]/2, C(0), C(1), C(2)};
  }

  vectord alpha;

public:

  JMT(){}

  JMT(const vectord &start, const vectord &end, double horizon) {
    alpha = solve(start, end, horizon);
    cout<<"JMT coeffs: "; for(auto a: alpha) cout<<a<<" "; cout<<endl;
    cout<<"input end = "<<end[0]<<", JMT solve end = "<<(*this)(horizon)<<endl;
  }

  JMT(const JMT &jmt): alpha(jmt.alpha){}

  double operator() (double t) const {
    double res = 0;
    for(int i=5; i>=0; i--){
      res = res*t  + alpha[i];
    }
    return res;
  }

  // s(t) = sum(i=0..5) ai*t^i
  // v(t) = sum(i=1..5) i*ai*t^(i-1)
  inline double velocity(double t) const {
    double res;
    for(int i=5; i>=1; i--){
      res = res*t  + i*alpha[i];
    }
    return res;
  }

  // s(t) = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5
  // v(t) = a1 + 2*a2*t + 3*a3*t^2 + 4*a4*t^3 + 5*a5*t^4
  // a(t) = 2*a2 + 6*a3*t + 12*a4*t^2 + 20*a5*t^3
  inline double acceleration(double t) const {
    return 2*alpha[2] + 6*alpha[3]*t + 12*alpha[4]*t*t +20*alpha[5]*t*t*t;
  }

  // j(t) = a3*3*2*1 + a4*4*3*2*t + a5*5*4*3*t^2
  // j(t) = 6*a3 + 24*a4*t + 60*a5*t^2
  inline double jerk(double t) const {
    return 6*alpha[3] + 24*alpha[4]*t + 60*alpha[5]*sqr(t);
  }
  
  friend ostream& operator<<(ostream &os, const JMT &jmt){
    cout<<"ostream JMT"<<endl;
    for(auto a: jmt.alpha) os<<a<<" ";
    cout<<"ostream JMT end"<<endl;
  }
};

class State{
public:
  double p, v, a;
  bool init_p, init_v, init_a;
  State(): init_p(false), init_v(false), init_a(false){
    p=v=a=0;
  }

  void update(double new_p, double dt=0.02) {
    double new_v = (new_p - p)/dt;
    double new_a = (new_v - v)/dt;

    p = new_p;
    if(init_p) {
      v = new_v;
      init_v = true;

      if(init_v) {
        a = new_a;
      }
    }

    init_p = true;
  }

  void set_v(double new_v) {
    if(!init_v) {
      v = new_v;
      init_v = true;
    }
  }

  void advance(double jerk, double dt=0.02) {
    p += v*dt + a/2*dt*dt + jerk/6*dt*dt*dt;
    v += a*dt + jerk/2*dt*dt;
    a += jerk*dt;
  }
};


#endif