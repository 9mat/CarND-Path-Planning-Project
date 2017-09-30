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

inline double target_d(int lane, double lane_width) {
  // The outer-most lane sometimes gives "out-of-lane" error
  // This may be to mismatch between spline interpolation and straight interpolation
  // To remedy for this, I offset the target d for the outer-most lane by 30cm
  return lane*lane_width + lane_width/2 - (lane==2? 0.3: 0);
}


inline int get_lane(double d, double lane_width) {
  return int(floor(d/lane_width));
}

inline double veh_speed(const vectord &veh) {
  // veh[2] and veh[3] are the x and y component of the vehicle's velocity
  return sqrt(sqr(veh[2]) + sqr(veh[3]));
}


class BehaviorPlanner {
public:
  static constexpr double lane_width = 4;
  static constexpr double safe_distance_for_change_lane = 20;
  static constexpr double safe_time_for_change_lane = 5;

  // Current targeted lane
  int target_lane;

  // To avoid sudden lane changes, I impose a minimum time between 2 consecutive lane change
  // This varibale will count down the next time that a lane change is allowed
  int change_lane_count_down;

  BehaviorPlanner() {
    target_lane = -1;
    change_lane_count_down = 0;
  }


  // calculate the cost of changing to a lane
  double lane_cost(bool same_lane, double pos_s, const vector<vectord> &vehicles, const map_wp &Map, double discount) const {
    double nearest_behind_s = -1e300, nearest_front_s = 1e300;
    vectord nearest_behind_veh, nearest_front_veh;
    bool has_vehicle_front = false, has_vehicle_behind = false;

    //  find the nearest car in front and behine
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

    double nearest_front_speed = has_vehicle_front ? veh_speed(nearest_front_veh) : 25;
    double nearest_behind_speed = has_vehicle_behind ? veh_speed(nearest_behind_veh) : 0;

    // I prefer the nearest car in front to be as far away and as fast as possible
    // Thus the cost will be decreasing in distance and speed of the nearest front car
    double cost = -min(nearest_front_dist, 300.0)*nearest_front_speed;

    // If a lane change is required, for safety purpose, 
    // make sure there is no car within a certain distance
    // otherwise, return a very high cost
    if(!same_lane) {
      if(nearest_front_dist < safe_distance_for_change_lane) cost = 1e100;
      if(nearest_behind_dist < safe_distance_for_change_lane) cost = 1e100;
    }

    cost -= discount*abs(cost);

    printf("%10.2e %11.2e %9.2e %10.2e %9.2e\n", nearest_front_dist, nearest_behind_speed, nearest_behind_dist, nearest_behind_speed, cost);


    return cost;
  }

  
  // realize a lane change
  int change_lane(int lane) {
    if(lane != target_lane) {
      printf("CHANGE LANE\n============\n");
      // count down to the next time that another lane change may be allowed
      change_lane_count_down = 50;

      // change the targeted lane
      target_lane = lane;
    } else {
      printf("KEEP LANE\n============\n");
    }

    return lane;
  }



  int choose_lane(double cur_d, double pos_s, const vector<vectord> &vehicles, const map_wp &Map) {

    // printf("All vehicles\n");
    // for(auto&veh: vehicles) {
    //   for(auto &d: veh) {
    //     printf("%7.2e ", d);
    //   }
    //   printf("\n");
    // }

    // Still too close to the last lane change, continue to carry the lane change and settle to the new lane
    if(change_lane_count_down > 0) {
      printf("CHANGING LANE. Count down %3d...\n==================\n", change_lane_count_down);
      change_lane_count_down--;
      return target_lane;
    }

    // Another check, if the horizontal to the target lane is too large
    // indicating lane changing has not finished
    // continue the lane change before looking for another lane change
    if (target_lane >= 0) {
      double distance_to_lane = abs(cur_d - target_d(target_lane, lane_width));
      if (distance_to_lane > lane_width/3) {
        printf("CHANGING LANE. Distance to target lane = %7.2f\n==================\n", distance_to_lane);
        return target_lane;
      } 
    }


    // Maximum curvature in the next 50m
    double future_max_curvature = 0.0;
    for(double future_s = pos_s; future_s <= pos_s + 50; future_s += 2) {
      updatemax(future_max_curvature, Map.curvature_frenet(future_s));
    }

    printf("Max curvature in 50m = %7.4f\n", future_max_curvature);

    // normal acceleration = curvature*speed^2
    // curvature = normal acceleration/speed^2 < acceleration limit/speed^2
    if (target_lane>=0 && future_max_curvature > 10.0/sqr(22.0)/2) {
      printf("KEEP LANE due to high curvature road segment\n==================\n");
      return target_lane;
    }


    int cur_lane = get_lane(cur_d, lane_width);


    // find the vehicles in each lane
    vector< vector<vectord> > vehicles_in_lane(3, vector<vectord>());
    for(auto& veh: vehicles){
      int lane = get_lane(veh[5], lane_width);
      if(lane>=0 && lane<=2) {
        vehicles_in_lane[lane].push_back(veh);
      }
    }

    double min_cost = 1e300;
    int best_lane = cur_lane;


    // I will prioritize the middle lane (since it will give more choices of lane changing in the future)
    // I will also prioritize the current lane
    // This vector will indicate the degree of priority of each lane
    vectord discount(3);
    discount[1] = 0.2;
    if(cur_lane != 1) discount[cur_lane] = 0.1;

    vectord cost(3);

    // header of the reporting table
    printf("Lane front/dist front/speed back/dist back/speed      cost\n");
    printf("---- ---------- ----------- --------- ---------- ---------\n");

    for(int i=0; i<=2; i++) {

      // reporting lane
      printf("%4d ", i);

      cost[i] = lane_cost(i==cur_lane, pos_s, vehicles_in_lane[i], Map, discount[i]);

      // find the lane with minimum cost
      if(min_cost > cost[i]) {
        min_cost = cost[i];
        best_lane = i;
      }
    }

    // reporting table
    printf("---- ---------- ----------- --------- ---------- ---------\n");


    // If the best lane is the current lane or adjacent to the current lane
    // then change to the best lane
    if(abs(best_lane - cur_lane) <= 1) {
      printf("Best lane: %d\n", best_lane);
      return change_lane(best_lane);
    }

    // If the best lane is 2 lanes apart
    // We will consider a lane change in the direction towards to the best lane
    int temp_lane = cur_lane + (best_lane > cur_lane? 1: -1);

    // Adjust the cost of the intermediate lane to capture the potential benefit of changing to the best lane
    // The new cost will be weighted average of the old cost and the cost of the best lane
    // Note, this adjustment is only done when a lane change is possible (the cost is not too large)
    if(cost[temp_lane]<5e99){
      cost[temp_lane] = (min_cost*2 + cost[temp_lane])/3;
    }

    // Try to find the lane with the minimum cost again, 
    // but this time only consider lanes adjacent to the current lane
    min_cost = cost[cur_lane];
    best_lane = cur_lane;
    for(int i=0; i<=2; i++) {
      if(abs(i-cur_lane)<=1 && cost[i] < min_cost) {
        min_cost = cost[i];
        best_lane = i;
      }
    }

    // Report Cost
    printf("Best lane: %d\n", best_lane);

    return change_lane(best_lane);
  }

};

class TrajectoryPlanner {
public:
  static constexpr double change_lane_interval = 40;
  static constexpr double lane_width = 4;
  static constexpr double spline_knot_spacing = 5;

  pair<vectord, vectord> generate_trajectory(int lane, State state_x, State state_y, const vectord &pre_x, const vectord &pre_y, const vector<vectord> &vehicles, const map_wp &Map) {
    int path_size = min(int(pre_x.size()), 10);

    // including at most 10 points of the old path to the new trajectory
    vectord trj_x(pre_x.begin(), pre_x.begin()+path_size);
    vectord trj_y(pre_y.begin(), pre_y.begin()+path_size);

    // convert to frenet
    double cur_s, cur_d, x, y, alpha;
    Map.frenet(trj_x.back(), trj_y.back(), cur_s, cur_d);

    // forward direction along the way points
    double heading = Map.heading_frenet(cur_s);

    // global to local coordinate transformation
    CoordTransformation local_coord(state_x.p, state_y.p, heading);

    double dt = 0.02;
    double v = sqrt(sqr(state_x.v) + sqr(state_y.v));

    // lane change: after a distance of change_lane_interval,
    // the vehicle is expected to be on the new lane, and continue so for 20 more iterations
    cur_s += change_lane_interval;
    cur_d = target_d(lane, lane_width);

    for(int i=0; i<20; i++, cur_s += spline_knot_spacing) {
      Map.smooth_cartersian(cur_s, cur_d, x, y);
      trj_x.push_back(x);
      trj_y.push_back(y);
    }

    // convert the trajectory to local coordinate
    // so that x-coordinates are likely to be sorted
    // so that they can be used for spline fitting
    for(int i=0; i<trj_x.size(); i++){
      local_coord(trj_x[i], trj_y[i]);
    }

    // Spline fitting for smooth lane change
    spline spl;
    spl.set_points(trj_x, trj_y);


    // Output path includes the 10 points on the old path
    vectord path_x(pre_x.begin(), pre_x.begin()+path_size);
    vectord path_y(pre_y.begin(), pre_y.begin()+path_size);

    // Simulate vehicle moving according the first few points in the trajectory
    // Update the states accordingly
    for(int i=0; i<path_size; i++) {
      state_x.update(pre_x[i]);
      state_y.update(pre_y[i]);
    }

    // Convert the states to the local coordinates
    local_coord(state_x.p, state_y.p);
    local_coord(state_x.v, state_y.v);
    local_coord(state_x.a, state_y.a);

    // Note that there is no translation for velocity and acceleration, only rotation
    // But local_coord will do both translation and rotation
    // Thus, for velocation and acceleration, we need to translated them back to the original origin
    // after being transformaed by local_coord
    double Ox = 0, Oy = 0;
    local_coord(Ox, Oy);
    state_x.v -= Ox;
    state_x.a -= Ox;
    state_y.v -= Oy;
    state_y.a -= Oy;

    double cur_x = trj_x[path_size-1], cur_y = trj_y[path_size-1], car_s, car_d;
    Map.frenet(pre_x[path_size-1], pre_y[path_size-1], car_s, car_d);

    // Find the neares car in front in the new lane
    // This vehicle will act as the reference vehicle for the trajectory

    int car_lane = int(floor(car_d/lane_width));
    double nearest_front_dist = 1e300, nearest_front_speed = 22;
    for(auto& veh: vehicles){
      double veh_v = veh_speed(veh);
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
    }

    // printf("TRAJECTORY PLANNER: \n");
    // printf("  current lane = %d, s = %7.2e, v = %7.2e\n", car_lane, car_s, sqrt(sqr(state_x.v) + sqr(state_y.v)));
    // printf("  target lane = %d, front distance = %7.2e, front speed = %7.2e\n", car_lane, nearest_front_dist, nearest_front_speed);

    double v_limit = 22;
    double a_limit = 9;
    double j_limit = 50;
    double target_j = 30;

    while(path_x.size() < 150) {
      double slope = spl.deriv(1, state_x.p);
      double speed = state_x.v*sqrt(1+sqr(slope));
      double curvature = spl.deriv(2,state_x.p)/pow(1+sqr(slope),1.5);

      // Accerelation can be decomposed into 2 componenent
      // The normal component depends on the curvature and the current velocity, thus cannot be controlled
      //    normal acceleration = curvature*velicity^2
      // The tangent component is parellel to the velocity
      // The total acceleraton = sqrt( sqr(normal acceleration) + sqr(tangent acceleration))
      // To make sure that the accerlation is within the limit, we need to adject the tangent acceraleration
      // to be less than sqrt(sqr(limit) - sqr(normal acceleration))
      double normal_acc = curvature*sqr(speed);
      double target_acc = abs(normal_acc) < 4 ? sqrt(20 - sqr(normal_acc)) : 0.0;




      bool speed_up = false, slow_down = false;

      double ref_speed = 22;

      double future_max_curvature = 0.0;
      for(double next_x = state_x.p; next_x < state_x.p+100; next_x += 1) {
        double next_slope = spl.deriv(1, next_x);
        double next_curvature = spl.deriv(2,next_x)/pow(1+sqr(next_slope),1.5);
        updatemax(future_max_curvature, abs(next_curvature));
      }


      // printf("Curvature = %7.4f, Future curv = %7.4f, Normal acc = %7.2f\n", curvature, future_max_curvature, normal_acc);

      if(future_max_curvature > 1e-10) {
        updatemin(ref_speed, sqrt(a_limit/future_max_curvature/1.5));
      }


      assert(normal_acc <= a_limit);

      // The logic that determines whether we should speed up, slow down or keep the speed constant
      // Speed up when the front vehicle is really far away, or when the front vehicle is relatively far and faster than the ego vehicle
      // Slow down when the front vehicle is really close, or when the front vehicle is relatively clse and run slower than the ego vehicle
      // Otherwise, try to match the speed of the front vehicle
      if (nearest_front_dist > 20) speed_up = true, slow_down = false;
      else if(nearest_front_dist > 15 && speed < nearest_front_speed) speed_up = true, slow_down = false;
      else if(nearest_front_dist < 15 && speed > nearest_front_speed) slow_down = true, speed_up = false;
      else if(nearest_front_dist < 10) slow_down = true, speed_up = false;

      if(speed_up && speed < ref_speed) speed += target_acc*dt;
      else if(slow_down) speed -= target_acc*dt;
      else {
        if(speed < nearest_front_speed && speed < ref_speed) speed = min(speed+target_acc*dt, nearest_front_speed);
        else speed = max(speed-target_acc*dt, nearest_front_speed);
      }

      if(speed > ref_speed) speed -= target_acc*dt;




      // Make sure the speed is within the limit
      updatemin(speed, v_limit);
      updatemax(speed, 0.0);

      // New positon (in local coordinate)
      double x = state_x.p + speed*dt/sqrt(1+sqr(slope));
      double y = spl(x);

      // Update the states of the car
      state_x.update(x);

      // Conver from local to global coordinate
      local_coord.revert(x,y);

      // Update the distance to the nearest front vehicle (assume the front vehicle travelling at a constant speed)
      nearest_front_speed += (nearest_front_speed - speed)*dt;

      // Output the new global coordinates
      path_x.push_back(x);
      path_y.push_back(y);
    }

    return make_pair(path_x, path_y);
  }

};


#endif