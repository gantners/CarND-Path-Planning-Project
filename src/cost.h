#ifndef COST_H
#define COST_H
#include "vehicle.h"

using namespace std;

float calculate_cost(Vehicle vehicle, map<int, vector<Vehicle>> predictions, vector<Vehicle> trajectory);

float goal_distance_cost(Vehicle vehicle, vector<Vehicle> trajectory, map<int, vector<Vehicle>> predictions, map<string, float> data);

float inefficiency_cost(Vehicle vehicle, vector<Vehicle> trajectory, map<int, vector<Vehicle>> predictions, map<string, float> data);

float lane_speed(map<int, vector<Vehicle>> predictions, int lane);

map<string, float> get_helper_data(Vehicle vehicle, vector<Vehicle> trajectory, map<int, vector<Vehicle>> predictions);

#endif