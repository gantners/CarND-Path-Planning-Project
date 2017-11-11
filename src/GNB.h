#ifndef CLASSIFIER_H
#define CLASSIFIER_H
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;
#include "Eigen-3.3/Eigen/Core"
using Eigen::ArrayXd;

class GNB {
public:

    vector<string> possible_labels = {"left","keep","right"};

    ArrayXd left_means;
    ArrayXd left_sds;

    ArrayXd keep_means;
    ArrayXd keep_sds;

    ArrayXd right_means;
    ArrayXd right_sds;


    /**
      * Constructor
      */
    GNB();

    /**
     * Destructor
     */
    virtual ~GNB();

    void train(vector<vector<double> > data, vector<string>  labels);

    string predict(vector<double>);

};

#endif