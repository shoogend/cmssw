#ifndef BEND_CALC
#define BEND_CALC

#include <cmath>
#include <string>
#include <vector>

using namespace std;

double phi(double del_r, double invPt);

int module_type(int layer, double r);

double int_bend(double bend, int module);

double meter_bend(double bend, int module);

double get_theta(double x, double y, double z);

double radius(double x, double y);

double get_radian(double deg);

double get_tilt(double z, int layer);

double get_separation(double r, double z, int layer);

double pt_bend(double x, double y, double z, int layer, double invPt, double charge);

double stub_bend_calc(double pt, double charge, int module, double r, double z, int layer);

double new_bend(double inner, double outer, double x, double y, double z, int layer, bool isFlipped);

double bend_isFlipped_variations(double inner, double outer, double x, double y, double z, int layer, bool isFlipped, int variation);

double r2(double r, double z, double tilt);

double z2(double r, double z, double tilt);

#endif
