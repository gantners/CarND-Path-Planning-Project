#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>


#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "GNB.h"

using namespace std;

// for convenience
using json = nlohmann::json;

struct val_lc {
    double d_min, d_max;
    double s_min, s_max;
};

bool is_within(double d, double s, val_lc other) {
    return d > other.d_min && d < other.d_max && s > other.s_min && s < other.s_max;
}

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


double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y) {

    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for (int i = 0; i < maps_x.size(); i++) {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x, y, map_x, map_y);
        if (dist < closestLen) {
            closestLen = dist;
            closestWaypoint = i;
        }

    }

    return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y) {

    int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];
    double heading = atan2((map_y - y), (map_x - x));

    double angle = fabs(theta-heading);
    angle = min(2*pi() - angle, angle); // XXX bug fix

    if(angle > pi()/4) {
        closestWaypoint++;
        if (closestWaypoint == maps_x.size()) {
            closestWaypoint = 0; // XXX bug fix
        }
    }
   /* double angle = abs(theta - heading);

    if (angle > pi() / 4) {
        closestWaypoint++;
    }*/

    return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y) {
    int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

    int prev_wp;
    prev_wp = next_wp - 1;
    if (next_wp == 0) {
        prev_wp = maps_x.size() - 1;
    }

    double n_x = maps_x[next_wp] - maps_x[prev_wp];
    double n_y = maps_y[next_wp] - maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
    double proj_x = proj_norm * n_x;
    double proj_y = proj_norm * n_y;

    double frenet_d = distance(x_x, x_y, proj_x, proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000 - maps_x[prev_wp];
    double center_y = 2000 - maps_y[prev_wp];
    double centerToPos = distance(center_x, center_y, x_x, x_y);
    double centerToRef = distance(center_x, center_y, proj_x, proj_y);

    if (centerToPos <= centerToRef) {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for (int i = 0; i < prev_wp; i++) {
        frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
    }

    frenet_s += distance(0, 0, proj_x, proj_y);

    return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double>
getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y) {
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();

    double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s - maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
    double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

    double perp_heading = heading - pi() / 2;

    double x = seg_x + d * cos(perp_heading);
    double y = seg_y + d * sin(perp_heading);

    return {x, y};

}

int get_lane(double d){
    if(d < (2 + 4 * 0 + 2) && d > (2 + 4 * 0 - 2))
        return 0;
    else if(d < (2 + 4 * 1 + 2) && d > (2 + 4 * 1 - 2))
        return 1;
    else if(d < (2 + 4 * 2 + 2) && d > (2 + 4 * 2 - 2))
        return 2;
    else
        return -1;
}

vector<vector<double> > Load_State(string file_name)
{
    ifstream in_state_(file_name.c_str(), ifstream::in);
    vector< vector<double >> state_out;
    string line;


    while (getline(in_state_, line))
    {
        istringstream iss(line);
        vector<double> x_coord;

        string token;
        while( getline(iss,token,','))
        {
            x_coord.push_back(stod(token));
        }
        state_out.push_back(x_coord);
    }
    return state_out;
}
vector<string> Load_Label(string file_name)
{
    ifstream in_label_(file_name.c_str(), ifstream::in);
    vector< string > label_out;
    string line;
    while (getline(in_label_, line))
    {
        istringstream iss(line);
        string label;
        iss >> label;

        label_out.push_back(label);
    }
    return label_out;

}

void train_classifier(GNB &gnb){
    cout << "training classifier" << endl;
    vector< vector<double> > X_train = Load_State("../data/train_states.txt");
    vector< vector<double> > X_test  = Load_State("../data/test_states.txt");
    vector< string > Y_train  = Load_Label("../data/train_labels.txt");
    vector< string > Y_test   = Load_Label("../data/test_labels.txt");

   /* cout << "X_train number of elements " << X_train.size() << endl;
    cout << "X_train element size " << X_train[0].size() << endl;
    cout << "Y_train number of elements " << Y_train.size() << endl;*/

    gnb.train(X_train, Y_train);

   /* cout << "X_test number of elements " << X_test.size() << endl;
    cout << "X_test element size " << X_test[0].size() << endl;
    cout << "Y_test number of elements " << Y_test.size() << endl;*/

    /*int score = 0;
    for(int i = 0; i < X_test.size(); i++)
    {
        vector<double> coords = X_test[i];
        string predicted = gnb.predict(coords);
        if(predicted.compare(Y_test[i]) == 0)
        {
            score += 1;
        }
    }

    float fraction_correct = float(score) / Y_test.size();
    cout << "You got " << (100*fraction_correct) << " correct" << endl;*/

    cout << "Classifier ready." << endl;
}

int main() {
    uWS::Hub h;

    GNB gnb = GNB();

    train_classifier(gnb);

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;


    ifstream in_map_(map_file_.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line)) {
        istringstream iss(line);
        double x;
        double y;
        float s;
        float d_x;
        float d_y;
        iss >> x;
        iss >> y;
        iss >> s;
        iss >> d_x;
        iss >> d_y;
        map_waypoints_x.push_back(x);
        map_waypoints_y.push_back(y);
        map_waypoints_s.push_back(s);
        map_waypoints_dx.push_back(d_x);
        map_waypoints_dy.push_back(d_y);
    }

    double ref_vel = 0.0;
    int lane = 1;

    h.onMessage(
            [&gnb, &ref_vel, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &lane](
                    uWS::WebSocket <uWS::SERVER> ws, char *data, size_t length,
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

                            //cout << car_d << endl;
                            /* cout << "My car: (" << car_x << "," << car_y << ")"
                                     "(" << car_yaw << "," << car_speed << ")"
                                     "(" << car_s << "," << car_d << ")"
                                     << endl;*/

                            // Sensor Fusion Data, a list of all other cars on the same side of the road.
                            auto sensor_fusion = j[1]["sensor_fusion"];

                            int prev_size = previous_path_x.size();

                            if (prev_size > 0) {
                                car_s = end_path_s;
                            }

                            //logic

                            bool too_close = false;
                            int next_lane = lane;
                            bool f_occupying_left = false;
                            bool f_occupying_right = false;
                            double dist_back = 20.0;
                            double dist_front = 20.0;
                            double target_vel = 49.5;
                            double speed_step = .224;
                            double m_per_s = .02;

                            double min_dist = 100.0;

                            for (int i = 0; i < sensor_fusion.size(); i++) {
                                auto f_car = sensor_fusion[i];
                                double vx = f_car[3];
                                double vy = f_car[4];
                                double check_car_s = f_car[5];
                                double check_car_d = f_car[6];
                                double check_speed = sqrt(vx * vx + vy * vy);
                                double orientation = atan2(vx, vy)/(M_PI/180);

                                //predict s into the future for this car
                                check_car_s += ((double) prev_size * m_per_s * check_speed);

                                //Check cars movement speed
                                double speed_diff = check_car_s - (double)f_car[5];
                                double distance_to_my_car = check_car_s - car_s;

                                min_dist = min(min_dist,distance_to_my_car);

                                bool in_front = distance_to_my_car >= 0;
                                bool my_lane = lane == get_lane(check_car_d);

                                //check s values greater than mine and s gap
                                if (in_front) {
                                    if(my_lane){
                                        if(distance_to_my_car < dist_front && distance_to_my_car){
                                           /* cout << "Approaching car id " << f_car[0] << " - Collision in " << distance_to_car << endl;*/
                                            too_close = true;
                                            //cout << check_speed << ","  << vx <<  "," << vy << ", " << orientation << endl;
                                        }
                                    }
                                }

                                //if lange change possible change lanes, otherwise slow down
                                for (int j = 0; j < sensor_fusion.size(); j++) {
                                    auto f2_car = sensor_fusion[j];
                                    float f_s = f2_car[5];
                                    float f_d = f2_car[6];

                                    //Calculate area of safe lane change and check if any car is occupying this space
                                    val_lc l;
                                    l.s_min = car_s - dist_back;
                                    l.s_max = car_s + dist_front;

                                    switch (lane) {
                                        case 0:
                                            f_occupying_left |= true;
                                            l.d_min = 4.0;
                                            l.d_max = 8.0;
                                            f_occupying_right |= is_within(f_d, f_s, l);
                                            break;
                                        case 1:
                                            l.d_min = 0.0;
                                            l.d_max = 4.0;
                                            f_occupying_left |= is_within(f_d, f_s, l);
                                            l.d_min = 8.0;
                                            l.d_max = 12.0;
                                            f_occupying_right |= is_within(f_d, f_s, l);
                                            break;
                                        case 2:
                                            l.d_min = 4.0;
                                            l.d_max = 8.0;
                                            f_occupying_left |= is_within(f_d, f_s, l);
                                            f_occupying_right |= true;
                                            break;
                                    }
                                }
                            }

                            /**
                             * Notice, however, that the car went from 0 MPH to 56 MPH in a single 20 ms frame, causing a spike in acceleration.
                             * Acceleration is calculated by comparing the rate of change of average speed over .2 second intervals.
                             * In this case total acceleration at one point was as high as 75 m/s^2. Jerk was also very high.
                             * The jerk is calculated as the average acceleration over 1 second intervals.
                             * In order for the passenger to have an enjoyable ride both the jerk and the total acceleration should not exceed 10 m/s^2
                             *
                             * Part of the total acceleration is the normal component, AccN which measures the centripetal acceleration from turning.
                             * The tighter and faster a turn is made, the higher the AccN value will be.
                             */

                            if (too_close) {
                                if (f_occupying_left && f_occupying_right) {
                                    cout << "Lane change not possible." << endl;
                                    if(min_dist < 2.0) {
                                        cout << "Emergency brake!" << endl;
                                        ref_vel /= 2;
                                    }
                                    else{
                                        ref_vel -= speed_step;
                                    }

                                } else {
                                    lane = f_occupying_left ? lane + 1 : lane - 1;
                                    cout << "safe to change lane to " << lane << endl;
                                }
                            }
                                //Otherwise accelerate
                            else if (ref_vel < target_vel) {
                                ref_vel += speed_step;
                            }

                            //cout << "target_vel: " << target_vel << ", ref_vel: " << ref_vel << endl;

                            //Create a list of widely spaced (x,y) waypoints, evenly spaces at 30m
                            //later we will interpolate these waypoints with a spline and fill it in with more points that control speed.
                            vector<double> ptsx, ptsy;

                            //Reference states
                            //Either we will reference the starting point as where the car is or at the previous paths end point
                            double ref_x = car_x;
                            double ref_y = car_y;
                            double ref_yaw = deg2rad(car_yaw);

                            //if previous size is almost empty, use the car as a starting point
                            if (prev_size < 2) {
                                //use two points that make the path tangent to the car
                                double prev_car_x = car_x - cos(car_yaw);
                                double prev_car_y = car_y - sin(car_yaw);

                                ptsx.push_back(prev_car_x);
                                ptsx.push_back(car_x);

                                ptsy.push_back(prev_car_y);
                                ptsy.push_back(car_y);
                            } else {
                                //Redefine reference state as previous path end point
                                ref_x = previous_path_x[prev_size - 1];
                                ref_y = previous_path_y[prev_size - 1];

                                double ref_x_prev = previous_path_x[prev_size - 2];
                                double ref_y_prev = previous_path_y[prev_size - 2];
                                ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

                                //use two points that make the path tangent to the previous path's end point
                                ptsx.push_back(ref_x_prev);
                                ptsx.push_back(ref_x);

                                ptsy.push_back(ref_y_prev);
                                ptsy.push_back(ref_y);
                            }

                            //In Frenet add evenly m spaced points ahead of the starting reference
                            int waypoints = 3;
                            int waypoint_steps = 30;

                            for(int i = 0; i < waypoints; i++){
                                int step_s = (i + 1) * waypoint_steps;
                                vector<double> next_wp = getXY(car_s + step_s, (2 + 4 * lane), map_waypoints_s,
                                                               map_waypoints_x,
                                                               map_waypoints_y);
                                //cout << next_wp[0] << " , " << next_wp[1] << endl;
                                ptsx.push_back(next_wp[0]);
                                ptsy.push_back(next_wp[1]);
                            }

                            for (int i = 0; i < ptsx.size(); i++) {
                                //shift car reference angle to 0 degrees
                                double shift_x = ptsx[i] - ref_x;
                                double shift_y = ptsy[i] - ref_y;

                                ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
                                ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
                            }

                            tk::spline s;
                            s.set_points(ptsx, ptsy);

                            json msgJson;

                            vector<double> next_x_vals;
                            vector<double> next_y_vals;

                            //path size for tracking the past waypoints for smooth trajectory
                            int path_size = previous_path_x.size();
                            for (int i = 0; i < path_size; i++) {
                                next_x_vals.push_back(previous_path_x[i]);
                                next_y_vals.push_back(previous_path_y[i]);
                            }

                            //Calculate how to break up spline points so that we travel at our desired reference velocity
                            double target_x = 30.0;
                            double target_y = s(target_x);
                            double target_dist = sqrt(target_x * target_x + target_y * target_y);

                            double x_add_on = 0;

                            //Fill up the rest of our path planner after filling it with previous points, here we will always output 50 points
                            int desired_waypoints = 50;
                            int new_waypoints = desired_waypoints - path_size;

                            for (int i = 0; i < new_waypoints; i++) {
                                double N = (target_dist / (m_per_s * ref_vel / 2.24));
                                double x_point = x_add_on + (target_x / N);
                                double y_point = s(x_point);

                                x_add_on = x_point;

                                double x_ref = x_point;
                                double y_ref = y_point;

                                //rotate back to normal after rotating earlier
                                x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
                                y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

                                x_point += ref_x;
                                y_point += ref_y;

                                next_x_vals.push_back(x_point);
                                next_y_vals.push_back(y_point);
                            }

                            msgJson["next_x"] = next_x_vals;
                            msgJson["next_y"] = next_y_vals;

                            auto msg = "42[\"control\"," + msgJson.dump() + "]";

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

    h.onConnection([&h](uWS::WebSocket <uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket <uWS::SERVER> ws, int code,
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
