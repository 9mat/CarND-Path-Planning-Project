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
#include "map_wp.h"
#include "util.h"

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

int main() {
  // vectord aa0 = {0,2,3}, aa1={10,1,1};
  // JMT mmm(aa0, aa1, 2);
  // cout<<mmm(2)<<endl;
  // return 0;
  uWS::Hub h;

  string map_file_ = "../data/highway_map.csv";
  map_wp Map(map_file_);

  double timestamp = 0, target_v=0;

  h.onMessage([&Map, &timestamp, &target_v](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {

        cout<<endl<<endl<<"**************************************************"<<endl;

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

          double heading = 0, x, y;

          Map.smooth_cartersian(car_s, car_d, x, y, heading);


          car_speed *= 0.44704;

        	// Previous path data given to the Planner
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];

        	// Previous path's end s and d values 
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];

          double dt = 0.02; // 20 ms
          timestamp += dt;

          int path_size = previous_path_x.size();
          if(path_size > 10) path_size=10;

        	// Sensor Fusion Data, a list of all other cars on the same side of the road.
        	auto sensor_fusion = j[1]["sensor_fusion"];

        	json msgJson;
          
          double pos_x, pos_y, angle, pos_s, pos_d, vel =0, acc = 0, jerk = 0;

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

          for(int i=1; i<path_size; i++){
            double new_vel = distance(previous_path_x[i-1], previous_path_y[i-1], previous_path_x[i], previous_path_y[i])/dt;
            double new_acc = (new_vel-vel)/dt;

            jerk = (new_acc-acc)/dt;
            acc = new_acc;
            vel = new_vel;
          }

          CoordTransformation local_coord(pos_x, pos_y, angle);

          Map.frenet(pos_x, pos_y, pos_s, pos_d, heading);

          int cur_lane = int(floor(pos_d/4));

          cout<<"Car position: s = "<<car_s<<", d ="<<car_d;
          cout<<" in lane "<<cur_lane<<endl;
          cout<<" Car speed: "<<car_speed<<" target_v: "<<target_v<<endl;
          cout<<" Car yaw  : "<<car_yaw<<" car heading: "<<heading*180/M_PI<<endl;

          vectord nearest_behind(3,1e300), nearest_front(3, 1e300);
          vectord nearest_behind_v(3, 0), nearest_front_v(3, 1e300);

          double buffer_s = 10;
          double buffer_s_lane_change = 30;
          double lane_width = 4;


          for(auto& other: sensor_fusion) {
            double other_x = other[1], other_y = other[2], other_vx = other[3], other_vy = other[4], other_s = other[5], other_d = other[6];
            int lane = int(floor(other_d/lane_width));
            local_coord(other_x, other_y);

            cout<<"Car #"<<other[0]<<" in lane "<<lane<<" at local coord = ("<<other_x<<","<<other_y<<")"<<endl;

            if(lane<0 || lane >= 3) {
              continue;
            }

            // double other_dist = sqrt(sqr(other_x) + sqr(other_y));
            double other_dist = other_s - car_s;
            double other_speed = sqrt(sqr(other_vx) + sqr(other_vy));

            if(other_dist > Map.total_length/2) other_dist -= Map.total_length;
            if(other_dist < -Map.total_length/2) other_dist += Map.total_length;

            if(other_dist < 0 && fabs(other_dist) < nearest_behind[lane]) {
              nearest_behind[lane] = fabs(other_dist);
              nearest_behind_v[lane] = other_speed;
            }
            else if(other_dist > 0 && fabs(other_dist) < nearest_front[lane]) {
              nearest_front[lane] = fabs(other_dist);
              nearest_front_v[lane] = other_speed;
            }
          }

          for(int i=0; i<3; i++){
            cout<<"##### Lane "<<i<<"  nearest_front  = "<<nearest_front[i]<<" v = "<<nearest_front_v[i]<<endl;
            cout<<"              nearest_behind = "<<nearest_behind[i]<<" v = "<<nearest_behind_v[i]<<endl;
          }

          int best_lane = cur_lane;
          double best_lane_v = nearest_front_v[cur_lane];

          for(int lane = cur_lane-1; lane <= cur_lane+1; lane++) {
            if(lane < 0 || lane > 2) continue;
            double keep_lane_preference = (lane == cur_lane) ? 1.2 : 1;
            double safe_distance = (lane==cur_lane) ? buffer_s : buffer_s_lane_change;
            bool feasible = (nearest_front[lane] > safe_distance) && (nearest_behind[lane] > safe_distance);
            if(feasible && nearest_front_v[lane]*keep_lane_preference > best_lane_v) {
              best_lane = lane;
              best_lane_v = nearest_front_v[lane];
            }
          }


          double ref_s = nearest_front[cur_lane] - buffer_s*2;
          double ref_v = nearest_front_v[cur_lane];
          double ref_a = sqr(car_speed - ref_v)/(2*ref_s);

          cout<<"reference s = "<<ref_s<<endl;
          cout<<"reference v = "<<ref_v<<endl;
          cout<<"reference a = "<<ref_a<<endl;

          bool speed_up = true;
          bool slow_down = false;
          bool match_speed = false;

          if(ref_s<0) {
            if(ref_s < -buffer_s || target_v > ref_v+0.2){
              speed_up = false;
              slow_down = true;              
            } 
            else {
              speed_up = false;
              slow_down = false;
              match_speed = true;
            }
          }
          else if(target_v < ref_v) {
            speed_up = target_v < 22;
          }
          else {
            if(ref_a > 10) {
              speed_up = false;
              slow_down = true;
            }
            else if(ref_a < 8) {
              speed_up = target_v < 22;
            }
          }

          if(slow_down && target_v >= 0.2) target_v -= 0.2;
          if(speed_up && target_v < 22) target_v += 0.2;
          if(match_speed) target_v = ref_v;

          if(slow_down && target_v >= 0.2) cout<<"*** Slow down"<<endl;
          else if(speed_up && target_v < 22) cout<<"*** Speed up"<<endl;
          else cout<<"Keep speed unchanged"<<endl;


          // if(nearest_front[cur_lane] < 10) {
          //   target_v -= 0.2;
          // }
          // else if (nearest_front[cur_lane] < 20 && (target_v > nearest_front_v[cur_lane])) {
          //   target_v -= 0.2;
          // }
          // else if(target_v < 22) {
          //   target_v += 0.2;

          // }

          double change_lane_interval = 30;

          vectord trajectory_x, trajectory_y;

          // convert previous postition to local coordinate 
          for(int i=0; i<path_size; i++) {
            trajectory_x.push_back(previous_path_x[i]);
            trajectory_y.push_back(previous_path_y[i]);
            local_coord(trajectory_x[i], trajectory_y[i]);
            // if (trajectory_x[i] > 0) break;
          }

          // // current position, in local coordinate is (0,0)
          // trajectory_x.push_back(0);
          // trajectory_y.push_back(0);


          // cout<<"Previous path x: ";
          // for(int i=0; i<path_size; i++) cout<<previous_path_x[i]<<" "; cout<<endl;

          // cout<<"trajectory path x: ";
          // for(int i=0; i<trajectory_x.size(); i++) cout<<trajectory_x[i]<<" "; cout<<endl;

          // double x, y;
          for(double ds = 0; ds < 50; ds += 5){
            Map.smooth_cartersian(pos_s + change_lane_interval + ds, best_lane*4+2, x, y, heading);
            local_coord(x,y);
            trajectory_x.push_back(x);
            trajectory_y.push_back(y);
          }


          spline spl;
          spl.set_points(trajectory_x, trajectory_y);


        	vector<double> next_x_vals, next_y_vals;

          double cur_x = 0, cur_y = 0;

          for(int i = 0; i < path_size; i++) {     
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          double dist_inc = 0.8;

          // cout<<"planed trajectory:"<<endl;

          while(next_x_vals.size() < 50) {
            double slope = spl.deriv(1, cur_x);

            if(speed_up && 22 - vel > sqr(acc)/(2*45)) {
              if(acc + 45*dt < 9 && vel + acc*dt + 45*dt*dt < 22) acc += 45*dt;
              if(vel < 22 - acc*dt) vel += acc*dt;
            }
            else if(speed_up && 22 - vel < sqr(acc)/(2*45)) {
              if(acc - 45*dt > -9) acc -= 45*dt;
              if(vel > -acc*dt) vel += acc*dt;
            }

            if(slow_down && vel >= 0.2) vel -= 0.2;

            double dx = vel/sqrt(1+sqr(slope))*dt;

            cur_x += dx;
            cur_y = spl(cur_x);

            pos_x = cur_x, pos_y = cur_y;
            local_coord.revert(pos_x, pos_y);

            next_x_vals.push_back(pos_x);
            next_y_vals.push_back(pos_y);

            // cout<<"Slope = "<<slope<<", Before transform: ("<<x<<","<<y<<"), after : ("<<pos_x<<","<<pos_y<<")"<<endl;
          }


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
















































































