#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <algorithm>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"


using namespace std;

// for convenience
using json = nlohmann::json;
using vectord = vector<double>;
using spline = tk::spline;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

inline double sqr(double x) { return x*x; }

inline double distance(double x1, double y1, double x2, double y2) {
	return sqrt(sqr(x2-x1) + sqr(y2-y1));
}

template <typename T> inline T sum(const vector<T> &a){
  T accum = 0;
  for(T aa: a) accum += aa;
  return accum;
}

class map_wp {
public:

  vector<double> x, y, s, dx, dy, bearing;
  double center_x, center_y, initial_bearing;
  int n;

  // the widthe of the lane
  static constexpr double lane_width = 4; 

  double normalized_bearing(double alpha){
    alpha -= initial_bearing;
    while(alpha < 0) alpha += 2*M_PI;
    while(alpha > 2*M_PI) alpha -= 2*M_PI;
    return alpha;
  }

  void intialize(){
    n = x.size();

    // the center of the map
    center_x = sum(x)/n;
    center_y = sum(y)/n;

    // bearing from each way point to the center of the map
    for(int i=0; i<n; i++){
      bearing.push_back(atan2(y[i]-center_y, x[i]-center_x));
    }

    // normalize bearings with respect to the first way point
    // so that the bearings are sorted
    initial_bearing = bearing[0];

    for(double &alpha: bearing) alpha = normalized_bearing(alpha);
  }

  // increase and decrease wp index, wrapping around at the end/start
  inline int inc_wp(int j){ return j==n-1? 0: j+1;}
  inline int dec_wp(int j){ return j==0? n-1: j-1;}


  map_wp(string map_file): x(), y(), s(), dx(), dy(), bearing(){
    ifstream is(map_file.c_str(), ifstream::in);
    for(;;){
      double xi, yi, si, dxi, dyi;
      is>>xi>>yi>>si>>dxi>>dyi;
      if(is.eof()) break;
      x.push_back(xi);
      y.push_back(yi);
      s.push_back(si);
      dx.push_back(dxi);
      dy.push_back(dyi);
    }

    intialize();
  }

  // find the coordinate of the point in the center of each of the 6 lanes
  // that corresponds to the wp-th way point
  void lane_wp(int wp, int lane_id, double &wp_x, double &wp_y) {
    double dist_to_center = lane_id*lane_width + lane_width/2;
    wp_x = x[wp] + dist_to_center*dx[wp];
    wp_y = y[wp] + dist_to_center*dy[wp];
  }

  int nearest_wp(double px, double py) {
    double alpha = normalized_bearing(atan2(py-center_y, px-center_x));

    // find the fist bearing that is greater than alpha 
    // using binary search function upper_bound from C++ STL <algorithm>
    // (this is making use of the fact that the bearing vector is sorted)
    int j = (upper_bound(bearing.begin(), bearing.end(), alpha) - bearing.begin()) % n;

    double nearest_dist = 1e30, dist;

    // try moving forward to find better way point
    // because of the near convexity of the map, break immediately 
    // if we cannot improve the nearest distance
    for(; nearest_dist < (dist=distance(px,py,x[j],y[j])); j = inc_wp(j)){
      nearest_dist = dist;
    }

    // try moving backward to find better way point
    for(; nearest_dist < (dist=distance(px,py,x[j],y[j])); j = dec_wp(j)){
      nearest_dist = dist;
    }

    return j;
  }

  int next_wp(double px, double py, double yaw) {
    int wp = nearest_wp(px, py);

    // the second nearest point can be either wp+1 or wp-1
    int wp2 = inc_wp(wp), wp3 = dec_wp(wp);
    double dist2 = distance(px, py, x[wp2], y[wp2]);
    double dist3 = distance(px, py, x[wp3], y[wp3]);

    // if the second nearest point is wp2 (aka wp+1)
    // then the 2 nearest points will be wp and wp+1
    // we have to increase wp so that the 2 nearest points 
    // will alwasy be wp-1 and wp
    if(dist2 < dist3) wp = wp2;


    // (px,py) lies somewhat between wp-1 and wp
    // presumbly the car is heading forward (wp-1 to wp, wp+1, etc.)
    // the next way point must be wp
    // assuming perfect control,
    // we skip the checking of the yaw angle

    return wp;
  }

  void frenet(double px, double py, double theta, double &ps, double &pd){
    // index of the two nearest way points
    int idx_a = next_wp(px,py,theta);
    int idx_b = dec_wp(idx_a);

    // coordinate of the two nearest way points A, B and the car P = (px,py)
    Eigen::Vector2d A(x[idx_a], y[idx_a]), B(x[idx_b], y[idx_b]), P(px,py);


    // illustratioon
    //    P
    //   /|
    //  / |
    // B--H-----A

    // given P, find ps = BH and pd = PH

    // vector indicate road direction and car position (P) from point B
    Eigen::Vector2d segment_direction = A-B, car_direction = P-B;

    // length of AB
    double segment_length = segment_direction.norm();

    // distance from project H to point B (the length of BH)
    // BH = BP*cos(B) = BP*AB*cos(B)/AB = vec(BA) dot vec(BP) / AB
    double dist_to_projection = segment_direction.dot(car_direction)/segment_length;

    // vector BH is in the same directionas vector BA, so we only to rescale the length
    Eigen::Vector2d direction_to_projection = dist_to_projection/segment_length*segment_direction;

    // coordiate of H = coordiate of B + vector BH
    Eigen::Vector2d projection = B + direction_to_projection, center(center_x, center_y);

    ps = s[idx_b] + dist_to_projection;
    pd = (P-projection).norm();

    // to determine if the point is inside or outside of the lane (negative or postive d)
    // compare the distance to the center of the car and of the projection
    if((center-P).norm() < (center-projection).norm()) pd *= -1;
  }


  void cartesian(double ps, double pd, double &px, double &py) {
    double map_length = s[n-1] + distance(x[0],y[0],x[n-1],y[n-1]);
    ps = fmod(ps, map_length);

    // find the segment somewhat contain the point
    // by searching for the first s that exceeds ps
    // given that s is sorted, binary search can be employed
    int upper_end = (upper_bound(s.begin(), s.end(), ps) - s.begin())%n; 
    int lower_end = dec_wp(upper_end);

    ps -= s[lower_end];

    // illustratioon
    //    P
    //   /|
    //  / |
    // B--H-----A

    // given BH = ps and PH = pd, find P

    Eigen::Vector2d A(x[upper_end], y[upper_end]), B(x[lower_end], y[lower_end]);
    Eigen::Vector2d segment_direction = A-B;
    double segment_length = segment_direction.norm();

    // coordinate of the projection H
    Eigen::Vector2d projection = ps/segment_length*segment_direction + B;

    // if pd = PH is negligible, return the lane point (aka the projection H)
    if(fabs(pd)<1e-6) {
      px = projection(0); py = projection(1);
      return;
    }

    // normal vector of AB (aka unit vector perpendicular to AB)
    Eigen::Vector2d normal(-segment_direction(1), segment_direction(0)), center(center_x, center_y);
    normal /= normal.norm();

    // there are 2 possibilites for P
    Eigen::Vector2d P1 = projection + pd*normal, P2 = projection - pd*normal, P;

    // to choose between thsese 2 possibilities, look at the distance to map center
    double dist_center_to_P1 = (P1-center).norm(), dist_center_to_P2 = (P2-center).norm();

    // if d is positive and P1 is further away from center than P2, use P1
    // there are 2 more cases, but the following logical statement capture them all
    if( (pd>0) == (dist_center_to_P1 > dist_center_to_P2)) P = P1;
    else P = P2;

    px = P(0); py = P(1);

  }

  void smooth_cartersian(double ps, double pd, double &px, double &py) {
    double map_length = s[n-1] + distance(x[0],y[0],x[n-1],y[n-1]);
    ps = fmod(ps, map_length);

    // find a point close enougth to the car
    int j = (upper_bound(s.begin(), s.end(), ps) - s.begin())%n; 

    // we will use 5 way points before j and 5 way points after j 
    // for spline interpolation
    j-=5;
    if(j<0) j+=n;

    vector<double> local_x, local_y, local_s;
    for(int i=0; i<10; i++, j=inc_wp(j)){
      local_x.push_back(x[j]);
      local_y.push_back(y[j]);
      local_s.push_back(s[j]);
    }

    // fitting 2D splines
    spline spline_x, spline_y;
    spline_x.set_points(local_s, local_x);
    spline_y.set_points(local_s, local_y);

    // interpolate to obtain the feet of the projection of car onto lane
    Eigen::Vector2d projection(spline_x(ps), spline_y(ps));

    // tangent of the lane at the projection point is the first derivative at ps
    // aka, tangent = (dspline_x/ds, dspline_y/ds)
    // the normal vector is perpendicular to the tangent vector
    // thus normal = (-dspline_y/ds, dspline_x/ds)
    Eigen::Vector2d normal(-spline_y.deriv(1,ps), spline_x.deriv(1,ps));

    // normalize the normal vector
    normal /= normal.norm();

    // there are 2 possibilites for the car position (P)
    Eigen::Vector2d P1 = projection + pd*normal, P2 = projection - pd*normal, P, center(center_x, center_y);

    // to choose between thsese 2 possibilities, look at the distance to map center
    double dist_center_to_P1 = (P1-center).norm(), dist_center_to_P2 = (P2-center).norm();

    // if d is positive and P1 is further away from center than P2, use P1
    // there are 2 more cases, but the following logical statement capture them all
    if( (pd>0) == (dist_center_to_P1 > dist_center_to_P2)) P = P1;
    else P = P2;

    px = P(0); py = P(1);
  }
};



void test_map() {
  string map_file_ = "../data/highway_map.csv";

  map_wp Map(map_file_);

  int j = Map.next_wp(2134.1, 1632.3, 1);
  cout<<"Nearest way point: "<<j-1<<" ("<<Map.x[j-1]<<","<<Map.y[j-1]<<")";
  cout<<", "<<j<<" ("<<Map.x[j]<<","<<Map.y[j]<<")"<<endl;


  double s, d;
  Map.frenet(2138.1, 1632.3, 1, s, d);
  cout<<"frenet s="<<s<<", d="<<d<<endl;


  double x, y;
  Map.cartesian(3000, 2, x, y);
  cout<<"Cartesian x="<<x<<", y="<<y<<endl;

  Map.smooth_cartersian(3000, 2, x, y);
  cout<<"Smooth cartesian x="<<x<<", y="<<y<<endl;
}

int main() {
  uWS::Hub h;

  string map_file_ = "../data/highway_map.csv";
  map_wp Map(map_file_);

  h.onMessage([&Map](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
        	double car_x = j[1]["x"];
        	double car_y = j[1]["y"];
        	double car_s = j[1]["s"];
        	double car_d = j[1]["d"];
        	double car_yaw = j[1]["yaw"];
        	double car_speed = j[1]["speed"];

        	// Previous path data given to the Planner
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];
        	// Previous path's end s and d values 
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];

        	// Sensor Fusion Data, a list of all other cars on the same side of the road.
        	auto sensor_fusion = j[1]["sensor_fusion"];

        	json msgJson;
          
          int path_size = previous_path_x.size();
          double pos_x, pos_y, angle;

          if(path_size <= 1) {
            pos_x = car_x;
            pos_y = car_y;
            angle = deg2rad(car_yaw);
          }
          else {
            pos_x = previous_path_x[path_size-1];
            pos_y = previous_path_y[path_size-1];

            double pos_x2 = previous_path_x[path_size-2];
            double pos_y2 = previous_path_y[path_size-2];
            angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
          }

        	vector<double> next_x_vals, next_y_vals;

          for(int i = 0; i < path_size; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
          }

          double dist_inc = 0.4, pos_s, pos_d;
          Map.frenet(pos_x, pos_y, angle, pos_s, pos_d);

          for(int i=0; i<200; i++) {
            pos_s += dist_inc;
            double x, y;
            Map.smooth_cartersian(pos_s, 2, pos_x, pos_y);
            next_x_vals.push_back(pos_x);
            next_y_vals.push_back(pos_y);
          }

          cout<<"planed trajectory:"<<endl;
          for(int i=0; i<next_x_vals.size(); i++){
            cout<<"   ("<<next_x_vals[i]<<","<<next_y_vals[i]<<")"<<endl;
          }

          cout<<"*****************"<<endl;


        	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
        	msgJson["next_x"] = next_x_vals;
        	msgJson["next_y"] = next_y_vals;

        	auto msg = "42[\"control\","+ msgJson.dump()+"]";

        	//this_thread::sleep_for(chrono::milliseconds(1000));
        	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































