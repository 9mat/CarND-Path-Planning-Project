# Overview
Path planning usually consists of three main components: _Prediction_, _Behavior planning_ and _Trajectory generation_. For this project, due to simple simulated environment (highway driving) I will simplify away most of the _Prediction_ component, employ a simple heuristics for the _Behavior planning_ and focus on the _Trajectory generation_ component.

# Prediction
For simplicity, I assume all the other vehicles as travelling at a constant speed with no lane changing. Thus the predicted location of a vehicle will simply be the current position plus the advancement due to the constant velocity of that vehicle. In frenet coordinate, the state evolution is as follows:

d(t) = d(0)
s(t) = s(0) + v*t

Due to the simplicity of this prediction model, it will not be implemented explicitly as an independent module or class, but it will be rather incorporated into the calculation of other module.

# Behavior planner
The Behavior planner is implemented in the class `BehaviorPlanner` in the file `planner.h`.

The responsibity of the _Behavior planner_ is to choose a target lane  for the vehicle to move into in each planning iteration.

The _Behavior planner_ will assign a cost to each lane. The calculation of the cost is implemented in the member function `lane_cost` of the class `BehaviorPlanner`. The cost is calculated as follow:

- If it is not possible to move to the lane, return a very high cost (`1e100`); this is case for a lane change to lane with a vehicle within a distance of `20m` (specified by the const `safe_distance_for_change_lane`)
- Otherwise, the cost will be the negative of the product between the distance to the nearest car in front and the speed of the car in front; the cost is defined this way to prioritize lane with nearest car in front which is far away and travelling at high speed
- To prioritize the middle lane, there will be a discount of 20% if the lane is the middle lane
- To discourage frequent lane changing, there will be a discount of 10% for the current lane

After calculating the costs, the function `choose_lane` will determine the target lane for the vehicle to move into. The function works in the following manner:

- If the lane with lowest cost is the current lane or adjacent to current lane, return that lane as the target lane
- Otherwise (i.e. the lane with the lowest cost is at least 2 lane apart from the current lane), I will consider the possibility of moving to an intermiate lane by reducing the cost of the intermediate lane
    + The cost of the intermediate lane will be set ot the weighted average of the min cost and the its old cost, unless it is not possible to move to the intermediate lane, in which case the cost remain at the high cost (`1e100`)
    + Choose the lane with the minimum cost among the current lane and those lanes adjacent to the current lane, return as the target lane

To avoid sudden change of lanes, especially the case where the best lane keep changing in middle of making a lane change, I implement a logic that forbidden chaning of the target lane for a certain number of iterations. The logic is as follows:

- I maintain 2 persistent states in the `BehaviorPlanner` class: `target_lane` and `lane_change_count_down`
- Whenever initiating a lane change, I reset the `lane_change_count_down` to `50`, and set the `target_lane` to the best lane identified by the logic above
- During the next 50 iterations, decrease the `lane_change_count_down` by 1 in each iteration until it reaches 0
- In any planning iteration, if the `lane_change_count_down` is greater than zero, don't do any lane change; i.e. return the `target_lane` of the previous iteration as the new target lane


Moreover, to avoid chaning lane on highly curved road segment, which will lead to a high normal acceleration which is not avoidable, I will restrict lane change not to be done if curvature in the next `50m` exceeds a certain threshold. Given that `normal acceleration = curvature * (speed)^2`, I choose the threshold to be `(acceleration limit)/(speed limit)^2*0.5` with 0.5 a safety margin.

# Trajectory Planning
The _Trajectory Planning_ is implemented in the class `TrajectoryPlanner` from the file `planner.h`

The responsility of the `TrajectoryPlanner` is to generate a trajectory connecting the current postition of the vehicle to the target lane that was identified by the `BehaviorPlanner`.

The Trajectory planning involves 2 stages:

- Stage 1: Generating a smooth path from the current location to the target plane (line `213` to line `249` of the file `planner.h`)
- Stage 2: Space the vehicle along the path over time so that the velocity, acceleration and jerk conform to the specified limits (line `253` to line `369` of the file `planner.h`)


## Generating a smooth trajectory
- Selection anchor points: I select points along the potential trajectory as follows:
    + Up to 10 points from the previous path
    + 20 points along the target lane, the first point is placed _30m_ (specified by the `constexpr` `change_lane_interval`) from the last point of the previous path; these points are placed _5m_ apart from each other (specified by the `constexpr` `spline_knot_spacing`)
- Convert these points to the local coordinate of the vehicle, this is done so that the `x`-coordinate is likely to be sorted, which is a requirement for 1D spline fitting
- Fit these point to a cubic spline using the `spline.h` library

## Space the vehicle over time along the path
- Maintain the variables `state_x` and `state_y` that record the states (position, velocity and acceleration) of the vehicle along the `x` and `y` coordinate 
- Maintain the variable `nearest_front_dist` that records the distant to the nearest other vehicle in front in the target lane
- Incrementally add 150 points for each of the next 150 time steps to the trajectory in the following manner:
    + Find `normal acceleration = curvature * (speed)^2`
    + Find the target tangent acceleration so that he total acceleration is within the limit, i.e. `target accelerattion = sqrt((limit)^2 - (normal acceleration)^2)`
    + Determine if the vehicle should speed up, slow down or keep constant speed 
        * Speed up if either `nearest_front_dist` is greater than 20m, or  `nearest_front_dist` is greater than 15m and the speed of the front vehicle is greater
        * Slow down if either `nearest_front_dist` is less than 10m, or `nearest_front_dist` is less than 15m and the speed of the front vehicle is slower
        * Otherwise, try to match the speed of the vehicle in front
    + Adjust the speed:
        * If speed up: increase the speed according to the target acceleration `v = v + (target acceleration)*dt`
        * If slow down: decrease the speed according to the target acceleration `v = v - (target acceleration)*dt`
        * If match speed: set speed to match the speed of the vehicle in front if such speed is within the target acceleration, otherwise, adjust the speed 
    + Move the vehicle to a new coordinate according to the new speed `x = x + vx * dt` with `vx = v/sqrt(1 + slope^2)`
    + Update the distance to the front vehicle `nearest_front_dist = nearest_front_dist + (front speed - ego speed)*dt`
    + Move to the next time step
        
# Misc

- I implemented the helper class `map` in `map_wp.h` to manage all the way-points operations, including Cartersian to Frenet conversion and vice versa
- I also implemented a smooth version of the Frenet to Cartersian conversion in the function `smooth_cartersian`, as I found the version using straight interpolation is rather jerky. The function `smooth_cartersian` find 10 way points around the point specified by the frenet coordinate `(s,d)`, fit 2 splines using `(x,s)` and `(y,s)` coordinates of those 10 points, and finally use these 2 splines to find the `x` and `y` coordinate corresponding to the the given frenet.
- I implemented a simple deterministic state (position-velocity-acceleration) updating mechanism in the class `State` in the file `util.h`. This can be improved by employing probabilistic state model.

# Results

[![CarND Path-Planning Result Video](http://img.youtube.com/vi/YOUTUBE_VIDEO_ID_HERE/0.jpg)](https://www.youtube.com/watch?v=5LsI5J_MGlI&feature=youtu.be "Result Video")

