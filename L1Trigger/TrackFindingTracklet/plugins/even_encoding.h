#include <vector>
#include <string>
#include <cmath>

using namespace std;

class Even_Encoding {
public:
  //vectors
  vector<string> four_bits{"1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111", "0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111"};
  vector<string> three_bits{"100","101","110","111","000","001","010","011"};


  // original encoding sheme (updated below), can toggle to test efficiency
  /*
  vector<double> layer1_0_encoding = bins(-2.8, 2.8, 8);
  vector<double> layer1_47_encoding = bins(-3.6, 3.6, 8);
  vector<double> layer1_60_encoding = bins(-3.8, 3.8, 8);
  vector<double> layer1_74_encoding = bins(-2.6, 2.6, 8);

  vector<double> layer2_0_encoding = bins(-3.1, 3.1, 8);
  vector<double> layer2_40_encoding = bins(-4.0, 4.0, 8);
  vector<double> layer2_55_encoding = bins(-3.6, 3.6, 8);
  vector<double> layer2_68_encoding = bins(-3.62, 3.62, 8);

  vector<double> layer3_0_encoding = bins(-4.1, 4.1, 8);
  vector<double> layer3_44_encoding = bins(-4.5, 4.5, 8);
  vector<double> layer3_60_encoding = bins(-3.85, 3.85, 8);

  vector<double>layer4_encoding_test2 = bins(-5.0, 5.0, 16);
  //vector<double> layer4_encoding_test = bins(-5.0, 5.0, 16);
  vector<double> layer4_encoding = bins(-5.0, 5.0, 16);
  vector<double> layer5_encoding = bins(-6.0, 6.0, 16);
  vector<double> layer6_encoding = bins(-7.5, 7.5, 16);
  
  vector<double> disk7_0_encoding = bins(-2.085, 2.085, 8);
  vector<double> disk7_1_encoding = bins(-3.09, 3.09, 8);
  vector<double> disk7_2_encoding = bins(-4.3, 4.3, 8);
  vector<double> disk7_3_encoding = bins(-6.4, 6.4, 8);
  vector<double> disk7_4_encoding = bins(-6.6, 6.6, 8);
  vector<double> disk7_5_encoding = bins(-6.0, 6.0, 16);
  vector<double> disk7_6_encoding = bins(-6.0, 6.0, 16);
  
  vector<double> disk8_0_encoding = bins(-2.07, 2.07, 8);
  vector<double> disk8_1_encoding = bins(-3.1, 3.1, 8);
  vector<double> disk8_2_encoding = bins(-3.57, 3.57, 8);
  vector<double> disk8_3_encoding = bins(-4.07, 4.07, 8);
  vector<double> disk8_4_encoding = bins(-6.0, 6.0, 8);
  vector<double> disk8_5_encoding = bins(-5.3, 5.3, 16);
  vector<double> disk8_6_encoding = bins(-5.5, 5.5, 16);

  vector<double> disk9_0_encoding = bins(-2.6, 2.6, 8);
  vector<double> disk9_1_encoding = bins(-3.06, 3.06, 8);
  vector<double> disk9_2_encoding = bins(-4.06, 4.06, 8);
  vector<double> disk9_3_encoding = bins(-4.54, 4.54, 16);
  vector<double> disk9_4_encoding = bins(-5.0, 5.0, 16);
  
  vector <double> disk10_0_encoding = bins(-2.5, 2.5, 8);
  vector <double> disk10_1_encoding = bins(-3.05, 3.05, 8);
  vector <double> disk10_2_encoding = bins(-3.55, 3.55, 8);
  vector <double> disk10_3_encoding = bins(-4.0, 4.0, 16);
  vector <double> disk10_4_encoding = bins(-4.41, 4.41, 16);

  vector <double> disk11_0_encoding = bins(-2.5, 2.5, 8);
  vector <double> disk11_1_encoding = bins(-3.04, 3.04, 8);
  vector <double> disk11_2_encoding = bins(-3.54, 3.54, 8);
  vector <double> disk11_3_encoding = bins(-4.5, 4.5, 16);
  vector <double> disk11_4_encoding = bins(-3.84, 3.84, 16);
  */

  // updated encoding scheme based on 100, removed outlier min/max
  vector<double> layer1_0_encoding = bins2(4.980, 8);
  vector<double> layer1_47_encoding = bins2(3.682, 8);
  vector<double> layer1_60_encoding = bins2(3.769, 8);
  vector<double> layer1_74_encoding = bins2(2.383, 8);

  vector<double> layer2_0_encoding = bins2(4.983, 8);
  vector<double> layer2_40_encoding = bins2(4.10, 8);
  vector<double> layer2_55_encoding = bins2(3.435, 8);
  vector<double> layer2_68_encoding = bins2(3.462, 8);

  vector<double> layer3_0_encoding = bins2(5.729, 8);
  vector<double> layer3_44_encoding = bins2(4.467, 8);
  vector<double> layer3_60_encoding = bins2(3.64, 8);

  vector<double>layer4_encoding_test2 = bins2(4.975, 16); 
  //vector<double> layer4_encoding_test = bins(-5.0, 5.0, 16); // this encoding vector occasionally did not work for unknown reasons, hence the test2 vector
  vector<double> layer4_encoding = bins2(4.975, 16);
  vector<double> layer5_encoding = bins2(5.970, 16);
  vector<double> layer6_encoding = bins2(7.463, 16);
  
  vector<double> disk7_0_encoding = bins2(1.750, 8);
  vector<double> disk7_1_encoding = bins2(2.980, 8);
  vector<double> disk7_2_encoding = bins2(3.704, 8);
  vector<double> disk7_3_encoding = bins2(4.652, 8);
  vector<double> disk7_4_encoding = bins2(6.126, 8);
  vector<double> disk7_5_encoding = bins2(4.770, 16);
  vector<double> disk7_6_encoding = bins2(5.651, 16);
  
  vector<double> disk8_0_encoding = bins2(1.477, 8);
  vector<double> disk8_1_encoding = bins2(2.735, 8);
  vector<double> disk8_2_encoding = bins2(3.429, 8);
  vector<double> disk8_3_encoding = bins2(3.955, 8);
  vector<double> disk8_4_encoding = bins2(3.892, 8);
  vector<double> disk8_5_encoding = bins2(3.829, 16);
  vector<double> disk8_6_encoding = bins2(5.236, 16);

  vector<double> disk9_0_encoding = bins2(2.298, 8);
  vector<double> disk9_1_encoding = bins2(2.763, 8);
  vector<double> disk9_2_encoding = bins2(3.830, 8);
  vector<double> disk9_3_encoding = bins2(4.528, 16);
  vector<double> disk9_4_encoding = bins2(4.584, 16);
  
  vector <double> disk10_0_encoding = bins2(1.754, 8);
  vector <double> disk10_1_encoding = bins2(2.903, 8);
  vector <double> disk10_2_encoding = bins2(3.346, 8);
  vector <double> disk10_3_encoding = bins2(3.897, 16);
  vector <double> disk10_4_encoding = bins2(4.339, 16);

  vector <double> disk11_0_encoding = bins2(2.485, 8);
  vector <double> disk11_1_encoding = bins2(2.485, 8);
  vector <double> disk11_2_encoding = bins2(3.438, 8);
  vector <double> disk11_3_encoding = bins2(4.282, 16);
  vector <double> disk11_4_encoding = bins2(3.520, 16);


  //methods    
  // creates a vector where each element is the midpoint of each bin in the encoding scheme
  vector<double> bins(double min, double max, int num_bins) {
    vector<double> bins (num_bins);
    bins[0] = min;
    double increment = (abs(max) + abs(min))/(num_bins);
    double i;
    for (i = 1; i < bins.size()+1; i++) {
      bins[i] = bins[i-1] + increment;
    }
    vector<double> mid (num_bins);
    for (i=0; i < bins.size(); i++) {
      mid[i] = (bins[i] + bins[i+1])/2;
	//  outfile2<<mid[i]<<",";
    }
    return mid;
  }

  vector<double> bins2(double edge, int num_bins) {
    vector<double> bins (num_bins);
    bins[0] = (-1.0*edge);
    double increment = (2*edge)/(num_bins);
    double i;
    for (i = 1; i < bins.size()+1; i++) {
      bins[i] = bins[i-1] + increment;
    }
    vector<double> mid (num_bins);
    for (i=0; i < bins.size(); i++) {
      mid[i] = (bins[i] + bins[i+1])/2;
	//  outfile2<<mid[i]<<",";
    }
    return mid;
  }

// function does a binary search to find the closest midpoint to a given bend value within the bins vector
  int binary_search(vector <double> v, double x) {
    double i, j, mid;
    i = 0;
    j = v.size() - 1;
    
    //corner cases
    if (x > v[j]) {
      return v.size() - 1;
    }
    else if (x < v[i]) {
      return 0;
    }
    else {
      while (i <= j) {
	mid = round((i+j)/2);
	if (v[mid] < x) {
	  i = mid + 1;
	}
	else {
	  j = mid -1;
	}
      }
    }
    double diff_i, diff_j;
    diff_i = abs(x - v[i]);
    diff_j = abs(x - v[j]);
    if (diff_i < diff_j) {
      return i;
    }
    else {
      return j;
    }
  }

// function that encodes bends for layers 1-3, based on the layer and tilt
  string encode_layers_123(double bend, int layer, double tilt) {
    int index;
    if (layer == 0) {
      if (tilt == 0) {
	index = binary_search(layer1_0_encoding, bend);
      }
      else if (tilt == 47.0) {
	index = binary_search(layer1_47_encoding, bend);
    }
      else if (tilt == 60.0) {
	index = binary_search(layer1_60_encoding, bend);
      }
      else {
      index = binary_search(layer1_74_encoding, bend);
      }
    }
    else if (layer == 1) {
      if (tilt == 0) {
      index = binary_search(layer2_0_encoding, bend);
      }
      else if (tilt == 40.0) {
	index = binary_search(layer2_40_encoding, bend);
      }
      else if (tilt == 55.0 ) {
      index = binary_search(layer2_55_encoding, bend);
      }
      else {
	index = binary_search(layer2_68_encoding, bend);
      }
    }
    else {
      if (tilt == 0) {
	index = binary_search(layer3_0_encoding, bend);
      }
      else if (tilt == 44.0) {
	index = binary_search(layer3_44_encoding, bend);
      }
      else {
	index = binary_search(layer3_60_encoding, bend);
      }
    }
    return three_bits[index];
  }
  

  // function encodes bend in layers 4-6, based on layer
  string encode_layers_456(double bend, int layer) {
    //string even_encoding::encode_layers_456(double bend, int layer) {
    int index;
    if (layer == 3) {
      index = binary_search(layer4_encoding, bend);
    }
    else if (layer == 4) {
      index = binary_search(layer5_encoding, bend);
    }
    else {
      index = binary_search(layer6_encoding, bend);
    }
    return four_bits[index];
  }

// function encodes bend in disks, based on layer and radius
  string encode_disks(double bend, int layer, double r) {
    int index;
    if (layer == 6) {
      if (r < 0.322) {
	index = binary_search(disk7_0_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.406) {
	index = binary_search(disk7_1_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.489) {
	index = binary_search(disk7_2_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.568) {
	index = binary_search(disk7_3_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.651) {
	index = binary_search(disk7_4_encoding, bend);
	return three_bits[index];
    }
      else if (r < 0.837) {
	index = binary_search(disk7_5_encoding, bend);
	return four_bits[index];
      }
      else {
	index = binary_search(disk7_6_encoding, bend);
	return four_bits[index];
      }
    }
    else if (layer == 7) {
      if (r < 0.322) {
	index = binary_search(disk8_0_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.406) {
	index = binary_search(disk8_1_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.489) {
	index = binary_search(disk8_2_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.568) {
	index = binary_search(disk8_3_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.651) {
	index = binary_search(disk8_4_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.837) {
	index = binary_search(disk8_5_encoding, bend);
	return four_bits[index];
      }
      else {
	index = binary_search(disk8_6_encoding, bend);
	return four_bits[index];
      }
    }
    else if (layer == 8) {
      if (r < 0.418) {
	index = binary_search(disk9_0_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.505) {
	index = binary_search(disk9_1_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.63) {
	index = binary_search(disk9_2_encoding, bend);
	return three_bits[index];
    }
      else if (r < 0.823) {
	index = binary_search(disk9_3_encoding, bend);
	return four_bits[index];
      }
      else {
	index = binary_search(disk9_4_encoding, bend);
	return four_bits[index];
      }
    }
    else if (layer == 9) {
      if (r < 0.418) {
	index = binary_search(disk10_0_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.505) {
	index = binary_search(disk10_1_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.63) {
	index = binary_search(disk10_2_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.823) {
	index = binary_search(disk10_3_encoding, bend);
	return four_bits[index];
      }
      else {
	index = binary_search(disk10_4_encoding, bend);
	return four_bits[index];
      }
    }
    else {
      if (r < 0.418) {
	index = binary_search(disk11_0_encoding, bend);
      }
      if (r < 0.505) {
	index = binary_search(disk11_1_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.63) {
	index = binary_search(disk11_2_encoding, bend);
	return three_bits[index];
      }
      else if (r < 0.823) {
	index = binary_search(disk11_3_encoding, bend);
	return four_bits[index];
      }
      else {
	index = binary_search(disk11_4_encoding, bend);
	return four_bits[index];
      }
    }
  }

// top-level encoding function which calls the other encoding functions based on layer
  string encode_top(double bend, int layer, double tilt, double r) {
    if (layer==0||layer==1||layer==2) {
      return encode_layers_123(bend, layer, tilt);
    }
    else if (layer==3||layer==4||layer==5) {
      return encode_layers_456(bend, layer);
    }
    else {
      return encode_disks(bend, layer, r);
    }
  }

} even_encoding;


  
