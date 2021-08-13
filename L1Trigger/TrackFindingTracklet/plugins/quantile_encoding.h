#include <vector>
#include <string>
#include <cmath>

using namespace std;

class Quantile_Encoding {
public:
  //vectors
  vector<string> four_bits{"1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111", "0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111"};
  vector<string> three_bits{"100","101","110","111","000","001","010","011"};

  //quantiles encoding scheme based on all_layer_100.txt, my bend calculation; vector elements represent edges of bins
  // layers 1-3 can be improved now that the error with the isFlipped has been resolved (no more bends of 12 for example)
  vector<double> layer1_0_q{-12.697, -1.706, -1.198, -0.67, 0.0, 0.813, 1.319, 1.86, 12.673};
  vector<double> layer1_47_q{-3.707, -1.409, -0.967, -0.58, 0.0, 0.607, 1.057, 1.47, 3.658};
  vector<double> layer1_60_q{-3.768, -1.558, -1.089, -0.655, 0.0, 0.679, 1.106, 1.567, 3.771};
  vector<double> layer1_74_q{-2/611, -0.917, -0.591, -0.319, 0.0, 0.324, 0.623, 0.932, 2.67};

  vector<double> layer2_0_q{-6.304, -1.728, -1.222, -0.747, 0.0, 0.828, 1.28, 1.714, 6.537};
  vector<double> layer2_40_q{-4.131, -2.008, -1.545, -1.062, 0.0, 1.041, 1.55, 1.96, 4.069};
  vector<double> layer2_55_q{-3.302, -1.555, -1.175, -0.712, 0.0, 0.923, 1.212, 1.58, 3.57};
  vector<double> layer2_68_q{-3.554, -1.622, -1.172, -0.744, 0.0, 0.743, 1.169, 1.597, 3.372};

  vector<double> layer3_0_q{-5.548, -2.242, -1.731, -1.13, 0.0, 1.08, 1.736, 2.242, 5.91};
  vector<double> layer3_44_q{-4.54, -2.615, -2.087, -1.472, 0.0, 1.459, 2.09, 2.581, 4.393};
  vector<double> layer3_60_q{-4.043, -1.863, -1.394, -0.978, 0.0, 1.004, 1.425, 1.839, 3.638};
  
  vector<double> layer4_q{-4.959, -3.85, -3.497, -3.18, -2.854, -2.493, -2.031, -1.271, 0.0, 1.326, 2.078, 2.565, 2.929, 3.244, 3.518, 3.826, 4.993};
  vector<double> layer5_q{-5.961, -4.859, -4.394, -4.029, -3.642, -3.133, -2.564, -1.683, 0.0, 1.755, 2.637, 3.239, 3.7, 4.072, 4.438, 4.825, 5.959};
  vector<double> layer6_q{-7.427, -6.162, -5.64, -5.185, -4.623, -4.056, -3.313, -2.184, 0.0, 2.132, 3.289, 4.095, 4.719, 5.218, 5.71, 6.177, 7.5};

  vector<double> disk7_0_q{-1.678,-0.847, -0.534, -0.272, 0.0, 0.253, 0.496, 0.812, 1.823};
  vector<double> disk7_1_q{-2.926, -1.233, -0.865, -0.512, 0.0, 0.532, 0.89, 1.222, 3.032};
  vector<double> disk7_2_q{-3.849, -1.71, -1.265, -0.821, 0.0, 0.776, 1.228, 1.656, 3.561};
  vector<double> disk7_3_q{-4.815, -2.256, -1.82, -1.35, 0.0, 1.319, 1.866, 2.29, 4.491};
  vector<double> disk7_4_q{-6.557, -3.005, -2.404, -1.841, 0.0, 1.972, 2.426, 2.954, 5.696};
  vector<double> disk7_5_q{-4.85, -2.708, -2.408, -2.107, -1.824, -1.567, -1.334, -0.85, 0.0, 0.819, 1.288, 1.526, 1.717, 2.034, 2.321, 2.6, 4.691};
  vector<double> disk7_6_q{-5.639, -4.148, -3.657, -3.362, -2.929, -2.637, -2.272, -1.555, 0.0, 1.54, 2.346, 2.686, 3.068, 3.419, 3.685, 4.182, 5.664};
  
  vector<double> disk8_0_q{-1.891, -0.743, -0.4, -0.207, 0.0, 0.17, 0.364, 0.735, 1.814};
  vector<double> disk8_1_q{-2.91, -1.046, -0.761, -0.407, 0.0, -.42, 0.815, 1.128, 2.56};
  vector<double> disk8_2_q{-3.103, -1.339, -1.008, -0.679, 0.0, 0.68, 1.018, 1.365, 3.406};
  vector<double> disk8_3_q{-3.87, -1.86, -1.384, -0.937, 0.0, 0.98, 1.407, 1.911, 4.042};
  vector<double> disk8_4_q{-4.273, -2.421, -2.028, -1.57, 0.0, 1.592, 2.001, 2.4931, 3.511};
  vector<double> disk8_5_q{-3.924, -2.338, -1.985, -1.688, -1.499, -1.346, -1.112, -0.729, 0.0, 0.848, 1.131, 1.401, 1.544, 1.769, 2.023, 2.414, 3.735};
  vector<double> disk8_6_q{-5.372, -3.545, -3.21, -2.867, -2.551, -2.277, -1.911, -1.323, 0.0, 1.416, 1.878, 2.273, 2.864, 3.188, 3.533, 5.102};

  vector<double> disk9_0_q{-2.432, -0.961, -0.62, -0.319, 0.0, 0.344, 0.601, 0.924, 2.165};
  vector<double> disk9_1_q{-2.812, -1.341, -0.978, -0.652, 0.0, 0.608, 0.976, 1.314, 2.715};
  vector<double> disk9_2_q{-3.942, -1.858, -1.425, -0.983, 0.0, 0.989, 1.391, 1.839, 3.72};
  vector<double> disk9_3_q{-4.546, -2.703, -2.267, -1.986, -1.709, -1.415, -1.176, -0.812, 0.0, 0.867, 1.177, 1.435, 1.697, 1.932, 2.213, 2.689, 4.511};
  vector<double> disk9_4_q{-4.321, -3.04, -2.601, -2.363, -2.129, -1.889, -1.545, -1.04, 0.0, 1.14, 1.637, 1.931, 2.177, 2.375, 2.636, 3.013, 4.848};
  
  vector <double> disk10_0_q{-1.963 -0.86, -0.537, -0.304, 0.0, 0.312, 0.579, 0.837, 1.647};
  vector <double> disk10_1_q{-2.928, -1.168, -0.765, -0.491, 0.0, 0.485, 0.791, 1.183, 2.878};
  vector <double> disk10_2_q{-3.311, 01.584, -1.241, -0.774, 0.0, 0.814, 1.274, 1.599, 3.383};
  vector <double> disk10_3_q{-3.928, 02.389, -1.976, -1.73, -1.521, -0.931, -0.561, 0.0, 0.556, 0.921, 0.1237, 1.471, 1.701, 1.946, 2.376, 3.868};
  vector <double> disk10_4_q{-4.333, -2.388, -2.19, -1.928, -1.771, -1.598, -1.24, -0.986, 0.0, 1.031, 1.333, 1.661, 1.852, 2.044, 2.234, 2.556, 4.346};

  vector <double> disk11_0_q{-3.009, -1.016, -0.646, -0.394, 0.0, 0.422, 0.682, 1.05, 2.728};
  vector <double> disk11_1_q{-3.009, -1.016, -0.646, -0.394, 0.0, 0.422, 0.682, 1.05, 2.728};
  vector <double> disk11_2_q{-3.411, -1.389, -0.998, -0.602, 0.0, 0.61, 1.0, 1.393, 3.465};
  vector <double> disk11_3_q{-4.111, -2.496, -2.137, -1.77, -1.574, -1.34, -1.091, -0.682, 0.0, 0.671, 1.123, 1.415, 1.622, 1.872, 2.273, 2.553, 4.455};
  vector <double> disk11_4_q{-3.268, -2.082, -1.861, -1.704, -1.325, -1.167, -1.005, -0.753, 0.0, 0.791, 1.041, 1.239, 1.516, 1.757, 1.912, 2.183, 3.773};

  // function utilizes binary search to find closest bin
  int binary_search_q(vector<double> v, double x) {
    double i, j, mid;
    i = 0;
    j = v.size() - 1;
    
    //corner cases
    if (x > v[j]) {
      return v.size() - 2;
    }
    else if (x < v[i]) {
    return 0;
    }
    else {
      while (i <= j) {
	mid = round((i+j)/2);
	if (v[mid] == x){
	  return mid;
	}
	else {
	  if (v[mid] < x) {
	    i = mid + 1;
	  }
	  else {
	    j = mid -1;
	  }
	}
      }
      if (i < j) {
	return i;
      }
      else {
	return j;
      }
    }
  }
  
  // encodes given bend based on layer and tilt angle for layers 1-3
  string encode_layers_123_q(double bend, int layer, double tilt) {
    int index;
    if (layer == 0) {
      if (tilt == 0) {
	index = binary_search_q(layer1_0_q, bend);
      }
     else if (tilt == 47.0) {
       index = binary_search_q(layer1_47_q, bend);
     }
     else if (tilt == 60.0) {
       index = binary_search_q(layer1_60_q, bend);
     }
     else {
       index = binary_search_q(layer1_74_q, bend);
     }
    }
    else if (layer == 1) {
      if (tilt == 0) {
	index = binary_search_q(layer2_0_q, bend);
      }
      else if (tilt == 40.0) {
	index = binary_search_q(layer2_40_q, bend);
      }
      else if (tilt == 55.0 ) {
	index = binary_search_q(layer2_55_q, bend);
      }
      else {
	index = binary_search_q(layer2_68_q, bend);
      }
    }
    else {
      if (tilt == 0) {
	index = binary_search_q(layer3_0_q, bend);
      }
      else if (tilt == 44.0) {
	index = binary_search_q(layer3_44_q, bend);
    }
      else {
	index = binary_search_q(layer3_60_q, bend);
      }
    }
    return three_bits[index];
  }


  // function encodes bend in layers 4-6, based on layer
  string encode_layers_456_q(double bend, int layer) {
    int index;
    if (layer == 3) {
      index = binary_search_q(layer4_q, bend);
    }
    else if (layer == 4) {
      index = binary_search_q(layer5_q, bend);
    }
    else {
      index = binary_search_q(layer6_q, bend);
    }
    return four_bits[index];
  }

  // function encodes bend in disks, based on layer and radius
  string encode_disks_q(double bend, int layer, double r) {
    int index;
    if (layer == 6) {
      if (r < 0.322) {
	index = binary_search_q(disk7_0_q, bend);
	return three_bits[index];
      }
      else if (r < 0.406) {
	index = binary_search_q(disk7_1_q, bend);
	return three_bits[index];
      }
      else if (r < 0.489) {
	index = binary_search_q(disk7_2_q, bend);
	return three_bits[index];
      }
      else if (r < 0.568) {
	index = binary_search_q(disk7_3_q, bend);
	return three_bits[index];
      }
      else if (r < 0.651) {
	index = binary_search_q(disk7_4_q, bend);
	return three_bits[index];
      }
      else if (r < 0.837) {
	index = binary_search_q(disk7_5_q, bend);
	return four_bits[index];
      }
      else {
	index = binary_search_q(disk7_6_q, bend);
	return four_bits[index];
      }
    }
    else if (layer == 7) {
      if (r < 0.322) {
	index = binary_search_q(disk8_0_q, bend);
	return three_bits[index];
      }
      else if (r < 0.406) {
	index = binary_search_q(disk8_1_q, bend);
	return three_bits[index];
      }
      else if (r < 0.489) {
       index = binary_search_q(disk8_2_q, bend);
       return three_bits[index];
      }
      else if (r < 0.568) {
	index = binary_search_q(disk8_3_q, bend);
	return three_bits[index];
      }
      else if (r < 0.651) {
       index = binary_search_q(disk8_4_q, bend);
       return three_bits[index];
      }
      else if (r < 0.837) {
	index = binary_search_q(disk8_5_q, bend);
	return four_bits[index];
      }
      else {
	index = binary_search_q(disk8_6_q, bend);
       return four_bits[index];
      }
    }
    else if (layer == 8) {
      if (r < 0.418) {
	index = binary_search_q(disk9_0_q, bend);
	return three_bits[index];
      }
      else if (r < 0.505) {
	index = binary_search_q(disk9_1_q, bend);
	return three_bits[index];
      }
      else if (r < 0.63) {
	index = binary_search_q(disk9_2_q, bend);
	return three_bits[index];
      }
      else if (r < 0.823) {
	index = binary_search_q(disk9_3_q, bend);
	return four_bits[index];
      }
      else {
	index = binary_search_q(disk9_4_q, bend);
	return four_bits[index];
      }
    }
    else if (layer == 9) {
      if (r < 0.418) {
	index = binary_search_q(disk10_0_q, bend);
	return three_bits[index];
      }
      else if (r < 0.505) {
	index = binary_search_q(disk10_1_q, bend);
	return three_bits[index];
      }
      else if (r < 0.63) {
	index = binary_search_q(disk10_2_q, bend);
	return three_bits[index];
      }
      else if (r < 0.823) {
	index = binary_search_q(disk10_3_q, bend);
	return four_bits[index];
      }
      else {
	index = binary_search_q(disk10_4_q, bend);
	return four_bits[index];
      }
    }
    else {
      if (r < 0.418) {
	index = binary_search_q(disk11_0_q, bend);
      }
      if (r < 0.505) {
	index = binary_search_q(disk11_1_q, bend);
	return three_bits[index];
      }
      else if (r < 0.63) {
	index = binary_search_q(disk11_2_q, bend);
	return three_bits[index];
      }
      else if (r < 0.823) {
	index = binary_search_q(disk11_3_q, bend);
	return four_bits[index];
      }
      else {
	index = binary_search_q(disk11_4_q, bend);
	return four_bits[index];
      }
    }
  }

  // top-level encoding function which calls the other encoding functions based on layer
  string encode_top_q(double bend, int layer, double tilt, double r) {
    if (layer==0||layer==1||layer==2) {
      return encode_layers_123_q(bend, layer, tilt);
    }
    else if (layer==3||layer==4||layer==5) {
      return encode_layers_456_q(bend, layer);
    }
    else {
      return encode_disks_q(bend, layer, r);
    }
  }

} quantile_encoding;

