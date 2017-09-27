#ifndef __BPLANNER__
#define __BPLANNER__

#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include <cmath>
#include <cassert>
#include "util.h"
#include "map_wp.h"
#include "spline.h"
using namespace std;

class BehaviorPlanner {
public:
  static constexpr double lane_width = 4;
  static constexpr double safe_distance_for_change_lane = 20;
  static constexpr double safe_time_for_change_lane = 5;
  static constexpr double weight_distance = 1/100;
  static constexpr double weight_speed = 1/100;

  enum ACTION {LANE_KEEP, LANE_CHANGE_LEFT, LANE_CHANGE_RIGHT};

  double lane_cost(bool same_lane, double pos_s, const vector<vectord> &vehicles, const map_wp &Map) const {
    double nearest_behind_s = -1e300, nearest_front_s = 1e300;
    vectord nearest_behind_veh, nearest_front_veh;
    bool has_vehicle_front = false, has_vehicle_behind = false;
    for(auto& veh: vehicles){
      double veh_s = veh[4];
      Map.standardize(pos_s, veh_s);
      if(veh_s < pos_s && veh_s > nearest_behind_s) {
        nearest_behind_s = veh_s;
        nearest_behind_veh = veh;
        has_vehicle_behind = true;
      }
      else if(veh_s >= pos_s && veh_s < nearest_front_s) {
        nearest_front_s = veh_s;
        nearest_front_veh = veh;
        has_vehicle_front = true;
      }
    }

    double nearest_front_dist = fabs(nearest_front_s-pos_s);
    double nearest_behind_dist = fabs(nearest_behind_s-pos_s);

    double nearest_front_speed = has_vehicle_front ? sqrt(sqr(nearest_front_veh[3]) + sqr(nearest_front_veh[4])) : 25;
    double nearest_behind_speed = has_vehicle_behind ? sqrt(sqr(nearest_behind_veh[3]) + sqr(nearest_behind_veh[4])) : 0;

    double cost = -nearest_front_dist*nearest_front_speed;
    if(same_lane) cost *= 1.1;

    printf("  frnt: dist = %7.2e, speed = %7.2e\n", nearest_front_dist, nearest_front_speed);
    printf("  back: dist = %7.2e, speed = %7.2e\n", nearest_behind_dist, nearest_behind_speed);
    printf("----------------------------------\n");

    if(!same_lane) {
      if(nearest_front_dist < safe_distance_for_change_lane) return 1e100;
      if(nearest_behind_dist < safe_distance_for_change_lane) return 1e100;
    }

    return cost;
  }

  int choose_lane(int cur_lane, double pos_s, const vector<vectord> &vehicles, const map_wp &Map) const {
    vector< vector<vectord> > vehicles_in_lane(3, vector<vectord>());

    int check_point = 0;

    for(auto& veh: vehicles){
      int lane = int(floor(veh[5]/lane_width));
      if(lane<0||lane>2) continue;
      vehicles_in_lane[lane].push_back(veh);
    }

    double min_cost = 1e300;
    int best_lane = cur_lane;

    for(int i=0; i<=2; i++){
      printf("Consider lane %d...\n", i);
      double cost = lane_cost(i==cur_lane, pos_s, vehicles_in_lane[i], Map);
      
      if(i==1) cost*=1.05;

      if(abs(i-cur_lane)<=1 && cost < min_cost){
        min_cost = cost;
        best_lane = i;
      }
    }

    printf("*** Best lane: %d, cost = %7.2e\n", best_lane, min_cost);

    return best_lane;
    // return 1;
  }

};

class TrajectoryPlanner {
public:
  static constexpr double change_lane_interval = 30;
  static constexpr double lane_width = 4;
  static constexpr double spline_knot_spacing = 5;

  pair<vectord, vectord> generate_trajectory(int lane, State state_x, State state_y, const vectord &pre_x, const vectord &pre_y, const vector<vectord> &vehicles, const map_wp &Map) {
    int path_size = min(int(pre_x.size()), 10);

    vectord trj_x(pre_x.begin(), pre_x.begin()+path_size);
    vectord trj_y(pre_y.begin(), pre_y.begin()+path_size);

    double cur_s, cur_d, heading, x, y;
    Map.frenet(trj_x.back(), trj_y.back(), cur_s, cur_d, heading);

    cur_s += change_lane_interval;
    cur_d = lane*lane_width + lane_width/2 - (lane==2? 0.3: 0);

    for(int i=0; i<10; i++, cur_s += spline_knot_spacing) {
      Map.smooth_cartersian(cur_s, cur_d, x, y, heading);
      trj_x.push_back(x);
      trj_y.push_back(y);
    }

    printf("Finish generating trajectory\n");

    double dt = 0.02;
    double v = sqrt(sqr(state_x.v) + sqr(state_y.v));
    double angle = atan2(state_y.v, state_x.v);

    CoordTransformation local_coord(state_x.p, state_y.p, angle);

    // printf("pre path length = %d, path length = %d\n", path_size, int(trj_x.size()));

    for(int i=0; i<trj_x.size(); i++){
      local_coord(trj_x[i], trj_y[i]);
      // printf("Point %d: %7.2f, %7.2f\n", i, trj_x[i], trj_y[i]);
    }

    printf("Planned trajectory: ");
    for(auto&x: trj_x) printf("%7.2f ", x); printf("\n");


    spline spl;
    spl.set_points(trj_x, trj_y);

    vectord path_x(pre_x.begin(), pre_x.begin()+path_size);
    vectord path_y(pre_y.begin(), pre_y.begin()+path_size);

    for(int i=0; i<path_size; i++) {
      state_x.update(pre_x[i]);
      state_y.update(pre_y[i]);
    }

    double Ox = 0, Oy = 0;
    local_coord(Ox, Oy);
    local_coord(state_x.p, state_y.p);
    local_coord(state_x.v, state_y.v);
    local_coord(state_x.a, state_y.a);

    state_x.v -= Ox;
    state_x.a -= Ox;
    state_y.v -= Oy;
    state_y.a -= Oy;

    double cur_x = trj_x[path_size-1], cur_y = trj_y[path_size-1], car_s, car_d;
    Map.frenet(pre_x[path_size-1], pre_y[path_size-1], car_s, car_d, heading);

    int car_lane = int(floor(car_d/lane_width));
    double nearest_front_dist = 1e300, nearest_front_speed = 22;
    for(auto& veh: vehicles){
      double veh_v = sqrt(sqr(veh[2]) + sqr(veh[3]));
      double veh_s = veh[4] + veh_v*path_size*dt;
      int veh_lane = int(floor(veh[5]/4));

      Map.standardize(car_s, veh_s);

      if(fabs(veh_s-car_s) > 100) continue;

      if(car_lane == veh_lane && veh_s>car_s) {
        if(veh_s-car_s < nearest_front_dist) {
          nearest_front_dist = veh_s - car_s;
          nearest_front_speed = veh_v;
        }
      }

      // printf("Veh: lane = %d, dist=%7.2f, speed=%7.2f\n", veh_lane, fabs(veh_s-car_s), veh_v);
    }

    printf("TRAJECTORY PLANNER: \n");
    printf("  current state: lane = %d, s = %7.2e, v = %7.2e\n", car_lane, car_s, sqrt(sqr(state_x.v) + sqr(state_y.v)));
    printf("  lane = %d, front distance = %7.2e, front speed = %7.2e\n", car_lane, nearest_front_dist, nearest_front_speed);

    double v_limit = 22;
    double a_limit = 9;
    double j_limit = 50;
    double target_j = 30;

    while(path_x.size() < 20) {
      double slope = spl.deriv(1, state_x.p);
      double speed = state_x.v*sqrt(1+sqr(slope));
      double curvature = spl.deriv(2,state_x.p)/pow(1+sqr(slope),1.5);
      double normal_acc = curvature*sqr(speed);
      double target_acc = abs(normal_acc) < 5 ? sqrt(25 - sqr(normal_acc)) : 0.0;
      // printf("Normal acc = %7.2f, target acc = %7.2f\n", normal_acc, target_acc);

      assert(speed >= -0.01);

      bool speed_up = false, slow_down = false;

      if (nearest_front_dist > 20) speed_up = true, slow_down = false;
      else if(nearest_front_dist > 15 && speed < nearest_front_speed) speed_up = true, slow_down = false;
      else if(nearest_front_dist < 15 && speed > nearest_front_speed) slow_down = true, speed_up = false;
      else if(nearest_front_dist < 10) slow_down = true, speed_up = false;

      if(speed_up) speed += target_acc*dt;
      else if(slow_down) speed -= target_acc*dt;
      // else speed += min(fabs(speed - nearest_front_speed),0.02)*(sgn(nearest_front_speed-speed));
      else {
        if(speed < nearest_front_speed) speed = min(speed+target_acc*dt, nearest_front_speed);
        else speed = max(speed-target_acc*dt, nearest_front_speed);
      }

      // printf("  step %d: ", path_x.size());
      // if(speed_up) printf("speed up\n");
      // else if(slow_down) printf("slow down\n");
      // else printf("keep speed\n");

      updatemin(speed, v_limit);
      updatemax(speed, 0.0);

      double x = state_x.p + speed*dt/sqrt(1+sqr(slope));
      double y = spl(x);

      state_x.update(x);


      local_coord.revert(x,y);

      nearest_front_speed += (nearest_front_speed - speed)*dt;





      // bool speed_up = true, slow_down = false;
      // double target_v = 22.0;
      // double curvature = spl.deriv(2,state_x.p)/pow(1+sqr(slope),1.5);
      // double tangent_x = 1/sqrt(1+sqr(slope));
      // double tangent_y = slope*tangent_x;

      // nearest_front_dist += nearest_front_speed*dt - speed*dt;

      // Map.frenet(cur_x, cur_y, cur_s, cur_d, heading);

      // if(nearest_front_dist < 15 && v > nearest_front_speed) {
      //   speed_up = false; slow_down = true;
      // }
      // else if(nearest_front_dist > 20 && v < nearest_front_speed) {
      //   slow_down = false;
      //   speed_up = (v+0.2) < nearest_front_speed;
      // }
      // else {
      //   speed_up = true;
      //   slow_down = false;
      // }

      // //     jerk
      // // --------- new acc
      // // |\`.    /
      // // | \  `./acc pojection
      // // |acc  /
      // // |   \/
      // // |   /tangent
      // // |  /
      // // | /
      // // |/


      // double horizon = 3*dt;

      // double target_x = state_x.p + state_x.v*horizon + state_x.a*sqr(horizon)/2;

      // double target_vx = state_x.v + state_x.v*horizon;
      // updatemin(target_vx, v_limit/sqrt(1+sqr(slope)));
      // updatemax(target_vx, 0.0);

      // double speed = sqrt(1+sqr(slope))*state_x.v;
      // double acc_y = spl.deriv(2, state_x.p)*sqr(state_x.v);
      // double acc = acc_y*tangent_y + state_x.a*tangent_x;
      // double normal_acc = acc_y*(-tangent_x) + state_x.a*tangent_y;


      // double target_acc = acc + horizon*target_j*(speed_up? 1 : (slow_down? -1 : (acc>target_j ? -1 : (acc<-target_j? 1 : 0))));
      // updatemin(target_acc, sqrt(sqr(a_limit) - sqr(normal_acc)));
      // updatemax(target_acc, -sqrt(sqr(a_limit) - sqr(normal_acc)));

      // double target_slope = spl.deriv(1, target_x);
      // double target_curvature = spl.deriv(2,target_x)/pow(1+sqr(slope),1.5);
      // double target_normal_acc = target_curvature*sqr(target_vx); // TODO: double check the sign


      // double target_acc_x = target_acc*tangent_x + target_normal_acc*tangent_y;


      // vectord start = {state_x.p, state_x.v, state_x.a};
      // vectord end = {target_x, target_vx, target_acc_x};

      // JMT jmt(start, end, horizon);

      // double x = jmt(dt);
      // double y = spl(x);

      // state_x.update(x, dt);

      // local_coord.revert(x,y);

      path_x.push_back(x);
      path_y.push_back(y);

    }

    // printf("Send path:\n");
    // for(int i=0; i<path_x.size(); i++){
    //   printf("(%7.2f,%7.2f)\n", path_x[i], path_y[i]);
    // }
    // printf("Finish TrajectoryPlanner\n");
    return make_pair(path_x, path_y);

  }

};


#endif