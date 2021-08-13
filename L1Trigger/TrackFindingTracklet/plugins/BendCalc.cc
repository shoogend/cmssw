#include <cmath>
#include <vector>
#include <string>

using namespace std;

// function to return separation angle between two hits
double phi(double del_r, double invPt) {
  return 0.3*del_r*3.8*0.5*invPt;
}

//function that returns 0 for PS modules, 1 for 2S modules
int module_type(int layer, double r) {
  if (layer==0||layer==1||layer==2) {
    return 0; }
  else if (layer==3||layer==4||layer==5) {
    return 1; }
  else if (layer==6||layer==7){
    if (r < 0.662) {
      return 0; }
    else {
      return 1; }
  }
  else {
    if (r<0.63) {
      return 0; }
    else {
      return 1; }
  }
}

// convert bend from int representation to meters and vice versa (works for all layers)
double int_bend(double bend, int module) {
  if (module == 0) {
    return bend/0.0001; }
  else {
    return bend/0.00009; }
}

double meter_bend(double bend, int module) {
  if (module ==0) {
    return bend*0.0001; }
  else {
    return bend*0.00009; }
} 

// calculates the angle theta (down from the z-axis), as a helper function for the general bend calculation
 double get_theta(double x, double y, double z) {
   double arg = sqrt((x*x)+(y*y))/z;
   return atan(arg);
 }

// calculates radius of hit, as helper for general bend calculation
double radius(double x, double y) {
   return sqrt((x*x)+(y*y));
}

// takes a number in degrees and returns it in radians
double get_radian(double deg) {
   return deg * M_PI/180;
}


// returns the tilt (in degrees) of a  module at a given position in the detector for any layer
double get_tilt(double z, int layer) {
  if (layer == 0) {
    if (z < 0.161) {
      return 0.0; }
    else if ( z < 0.2925) {
      return 47.0; }
    else if (z < 0.568) {
      return 60.0; }
    else {
      return 74.0; }
  }
  else if (layer == 1) {
    if (z < 0.246) {
      return 0.0; }
    else if (z < 0.405) {
      return 40.0; }
    else if (z < 0.6855) {
      return 55.0; }
    else {
      return 68.0; }
  }
  else if (layer == 2) {
    if (z < 0.3395) {
      return 0.0; }
    else if (z < 0.696) {
      return 44.0; }
    else {
      return 60.0; }
  }
  else if (layer == 3||layer==4||layer==5) {
    return 0.0; }
  else { //layers==6, 7, 8, 9
    return 90.0; }
}

// gives separation between the two sides of a module in m for any layer
double get_separation(double r, double z, int layer) {
  if (layer == 0) {
    if (z < 0.291) {
      return 0.0026; }
    else {
      return 0.004; }
  }
  else if (layer == 1) {
    if (z < 0.246) {
      return 0.0016; }
    else if (z < 0.6855) {
      return 0.0026; }
    else {
      return 0.004; }
  }
  else if (layer == 2) {
    if (z < 0.3395) {
      return 0.0016; }
    else {
      return 0.0026; }
  }
  else if (layer == 3||layer==4||layer==5) {
    return 0.0018; 
  }
  else if (layer==6||layer==7||layer==8) {
    if (r < 0.662) {
      return 0.004; }
    else {
      return 0.0018; }
  }
  else { //disks
    if (r<0.833) {
      return 0.004; }
    else {
      return 0.0018; }
  }
}

// function calculates the bend based on the pT, not on inner/outer position (used for debugging)
double pt_bend(double x, double y, double z, int layer, double invPt, double charge) {
  double del_phi, del_r, d, theta, beta, r;
  r = radius(x, y);
  beta = get_tilt(z, layer);
  d = get_separation(r, z, layer);
  theta = get_theta(x, y, z);
  del_r = d*sin(theta)/sin(theta+beta);
  del_phi = 0.57*invPt*del_r *charge;
  double bend = del_phi * r;
  if ((layer>=6)&&(z<0)) {
    bend = -1*bend;
  }
  return bend;
}

// recreated bend calculation from Stub.cc file (for debugging)
double stub_bend_calc(double pt, double charge, int module, double r, double z, int layer) {
  double bfield_{3.8112};  //B-field in T
  double c_{0.299792458};  //speed of light m/ns
  double dr = 0.18;
  double rinv = (charge * 0.01 * c_ * bfield_)/pt;
  double pitch;
  if (module == 0) {
    pitch = 0.01; // in cm then...
  }
  else {
    pitch = 0.009;
  }
  double bend = r * dr * 0.5 * rinv / pitch;
  if ((layer>=6)&&(z<0)) {
    bend = -1*bend;
  }
  return bend;
}


// gives bend for any hit within the detector
double new_bend(double inner, double outer, double x, double y, double z, int layer, bool isFlipped) {
  double r = radius(x, y);
  double theta = get_theta(x, y, z);
  double beta = get_radian(get_tilt(z, layer));
  double d = get_separation(r, z, layer);
  double dist = (508 - inner);
  double module = module_type(layer,r);
  double del_r;
  if (layer==0||layer==1||layer==2) {
    del_r = (d*sin(theta))/sin(theta+beta);
  } 
  else if (layer==3||layer==4||layer==5) {
    del_r = d;
  }
  else {
    del_r = d*tan(theta);
  }
  if (not isFlipped) {
    double diff = (inner-outer);
    if (module == 0) {
      return (diff*0.0001) - (((dist*0.0001)/r)*del_r);
    } 
    else {
      return (diff*0.00009) - (((dist*0.00009)/r)*del_r);
    }
  }
  else {
    double diff = (outer - inner);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) - (diff*0.0001);
    } else {
      return (((dist*0.00009)/r)*del_r) - (diff*0.00009);
    }
  }
}

// testing out different variation of the bend calculation (used for debugging)
double bend_isFlipped_variations(double inner, double outer, double x, double y, double z, int layer, bool isFlipped, int variation) {
  double r = radius(x, y);
  double theta = get_theta(x, y, z);
  double beta = get_radian(get_tilt(z, layer));
  double d = get_separation(r, z, layer);
  double module = module_type(layer,r);
  double del_r;
  if (layer==0||layer==1||layer==2) {
    del_r = (d*sin(theta))/sin(theta+beta);
  } 
  else if (layer==3||layer==4||layer==5) {
    del_r = d;
  }
  else {
    del_r = d*tan(theta);
  }
  double dist, diff;
  if (variation == 1) {
    dist = (508 - inner);
    diff = (outer - inner);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) - (diff*0.0001);
    } else {
      return (((dist*0.00009)/r)*del_r) - (diff*0.00009);
    }
  }
  else if (variation == 2) {
    dist = (inner - 508);
    diff = (outer - inner);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) - (diff*0.0001);
    } else {
      return (((dist*0.00009)/r)*del_r) - (diff*0.00009);
    }
  }
  else if (variation == 3) {
    dist = (508 - inner);
    diff = (inner - outer);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) - (diff*0.0001);
    } else {
      return (((dist*0.00009)/r)*del_r) - (diff*0.00009);
    }
  }
  else if (variation == 4){
    dist = (inner - 508);
    diff = (inner - outer);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) - (diff*0.0001);
    } else {
      return (((dist*0.00009)/r)*del_r) - (diff*0.00009);
    }
  }
  else if (variation == 5) {
    dist = (508 - inner);
    diff = (outer - inner);
    if (module == 0) {
      return ((diff*0.0001) - (((dist*0.0001)/r)*del_r));
	} 
    else {
      return ((diff*0.00009) - (((dist*0.00009)/r)*del_r));
	}
  }
  else if (variation == 6) {
    dist = (inner - 508);
    diff = (outer - inner);
    if (module == 0) {
      return ((diff*0.0001) - (((dist*0.0001)/r)*del_r));
	} 
    else {
      return ((diff*0.00009) - (((dist*0.00009)/r)*del_r));
	}
  }
  else if (variation == 7) {
    dist = (508 - inner);
    diff = (inner - outer);
    if (module == 0) {
      return ((diff*0.0001) - (((dist*0.0001)/r)*del_r));
	}
    else {
      return ((diff*0.00009) - (((dist*0.00009)/r)*del_r));
	}
  }
  else if (variation == 8){
    dist = (inner - 508);
    diff = (inner - outer);
    if (module == 0) {
      return ((diff*0.0001) - (((dist*0.0001)/r)*del_r));
	} 
    else {
      return ((diff*0.00009) - (((dist*0.00009)/r)*del_r));
	}
  }
  else if (variation == 9) {
    dist = (508 - inner);
    diff = (outer - inner);
    if (module == 0) {
      return ((diff*0.0001) + (((dist*0.0001)/r)*del_r));
    }
    else {
      return ((diff*0.0001) + (((dist*0.0001)/r)*del_r));
    }
  }
  else if (variation == 10) {
    dist = (inner - 508);
    diff = (outer - inner);
    if (module == 0) {
      return ((diff*0.0001) + (((dist*0.0001)/r)*del_r));
    }
    else {
      return ((diff*0.0001) + (((dist*0.0001)/r)*del_r));
    }     
  }
  else if (variation == 11) {
    dist = (508 - inner);
    diff = (inner - outer);
    if (module == 0) {
      return ((diff*0.0001) + (((dist*0.0001)/r)*del_r));
    }
    else {
      return ((diff*0.0001) + (((dist*0.0001)/r)*del_r));
    }
  }
  else if (variation == 12){
    dist = (inner - 508);
    diff = (inner - outer);
    if (module == 0) {
      return ((diff*0.0001) + (((dist*0.0001)/r)*del_r));
    }
    else {
      return ((diff*0.0001) + (((dist*0.0001)/r)*del_r));
    }
  }
  if (variation == 13) {
    dist = (508 - inner);
    diff = (outer - inner);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) + (diff*0.0001);
    }
    else {
      return (((dist*0.00009)/r)*del_r) + (diff*0.00009);
    }
  }
  else if (variation == 14) {
    dist = (inner - 508);
    diff = (outer - inner);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) + (diff*0.0001);
    } else {
      return (((dist*0.00009)/r)*del_r) + (diff*0.00009);    
    }
  }
  else if (variation == 15) {
    dist = (508 - inner);
    diff = (inner - outer);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) + (diff*0.0001);
    } else {
      return (((dist*0.00009)/r)*del_r) + (diff*0.00009);    
    }
  }
  else {
    dist = (inner - 508);
    diff = (inner - outer);
    if (module == 0) {
      return (((dist*0.0001)/r)*del_r) + (diff*0.0001);
    } else {
      return (((dist*0.00009)/r)*del_r) + (diff*0.00009);    
    }
  }
}

// calculates rotated coordinates for given module tilt angle, used for new isFlipped definition
double r2(double r, double z, double tilt) {
  double angle = get_radian(90.0 - tilt);
  return (-z*sin(angle)) + (r*cos(angle));
}

double z2(double r, double z, double tilt) {
  double angle = get_radian(90.0 - tilt);
  return (z*cos(angle)) + (r*sin(angle));
}

