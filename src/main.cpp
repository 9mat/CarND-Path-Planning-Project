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
#include "planner.h"

using namespace std;

// for convenience
using json = nlohmann::json;
using vectord = vector<double>;
using spline = tk::spline;


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
  BehaviorPlanner bplanner;
  TrajectoryPlanner tplanner;

  // State of the car, persistent across multiple planning iterations
  State state_x, state_y;

  double timestamp = 0, target_v=0;
  int count = 0;

  h.onMessage([&count, &Map, &timestamp, &target_v, &bplanner, &tplanner, &state_x, &state_y](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

          // convert from mph to m/s
          car_speed *= 0.44704;

          // update the states of the car
          state_x.update(car_x);
          state_y.update(car_y);

          state_x.set_v(car_speed*cos(car_yaw));
          state_y.set_v(car_speed*sin(car_yaw));

          double heading = 0, x, y;

          Map.smooth_cartersian(car_s, car_d, x, y, heading);

        	// Previous path data given to the Planner
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];

        	// Previous path's end s and d values 
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];

          double dt = 0.02; // 20 ms
          timestamp += dt;

          int path_size = previous_path_x.size();

          // To avoid over-dependence on previous path, I will use only 10 points from the previous path
          if(path_size > 10) path_size=10;

        	// Sensor Fusion Data, a list of all other cars on the same side of the road.
        	auto sensor_fusion = j[1]["sensor_fusion"];

          vector<vectord> vehicles(sensor_fusion.size());
          for(int i=0; i<sensor_fusion.size(); i++) {
            for(int k=1; k<7; k++){
              vehicles[i].push_back(double(sensor_fusion[i][k]));
            }
          }

          vectord pre_x, pre_y;
          for(int i=0; i<path_size; i++){
            pre_x.push_back(previous_path_x[i]);
            pre_y.push_back(previous_path_y[i]);
          }


          // Report car position and previous path information
          printf("Car @(%7.2f,%7.2f), s=%7.2f, speed = %7.2f\n", car_x, car_y, car_s, car_speed);
          printf("receive path: \n");
          for(int i=0; i<pre_x.size(); i++){
            if(i>0){  
              double v = distance(pre_x[i], pre_y[i], pre_x[i-1], pre_y[i-1])/dt;
              printf(" @(%7.2f,%7.2f), v = %7.2f\n", pre_x[i], pre_y[i], v);
            }
            else{
              printf(" @(%7.2f,%7.2f)\n", pre_x[i], pre_y[i]);              
            }
          }

          // If no previous, use the car car position to indicate the starting point of the future path
          if(path_size < 1){
            path_size = 1;
            pre_x.push_back(car_x);
            pre_y.push_back(car_y);
          }

          // call behaviour planner to choose the best lane
          int best_lane = bplanner.choose_lane(car_d, car_s, vehicles, Map);


          // call the trajectory planner to find a trajectory towards the best lane
          pair<vectord, vectord> next_path = tplanner.generate_trajectory(best_lane, state_x, state_y, pre_x, pre_y, vehicles, Map);


          json msgJson;
        	msgJson["next_x"] = next_path.first;
        	msgJson["next_y"] = next_path.second;

        	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          printf("Finish at timestamp %7.2g\n", timestamp);

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
















































































