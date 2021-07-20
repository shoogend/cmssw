#ifndef L1Trigger_TrackFindingTracklet_interface_Util_h
#define L1Trigger_TrackFindingTracklet_interface_Util_h

#include <sstream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

namespace trklet {

  //Converts string in binary to hex (used in writing out memory content)
  inline std::string hexFormat(const std::string& binary) {
    std::stringstream ss;

    unsigned int radix = 1, value = 0;
    for (int i = binary.length() - 1; i >= 0; i--) {
      if (binary.at(i) != '0' && binary.at(i) != '1')
        continue;
      value += (binary.at(i) - '0') * radix;
      if (radix == 8) {
        ss << std::hex << value;
        radix = 1;
        value = 0;
      } else
        radix <<= 1;
    }
    if (radix != 1)
      ss << std::hex << value;

    std::string str = ss.str() + "x0";
    std::reverse(str.begin(), str.end());
    return str;
  }

  inline double cosModuleTilt = 0.886454;
  inline double sinModuleTilt = 0.504148;

/*
  inline double tmzinner[3][13] =
    {{172.095, 217.118, 267.987, 315.184, 374.203, 443.696, 526.772, 611.087, 720.246, 850.617, 1004.205, 1182.332, 1200},
     {269.888, 321.982, 378.406, 431.314, 495.745, 567.262, 646.108, 725.236, 820.642, 927.596, 1048.134, 1181.542, 1200},
     {363.047, 416.152, 472.684, 533.009, 597.236, 665.109, 727.825, 804.833, 888.179, 978.285, 1075.169, 1179.356, 1200}};
*/

  inline double tmzinner[3][8] = {{150, 300, 450, 600, 750, 900, 1050, 1200},
                                  {150, 300, 450, 600, 750, 900, 1050, 1200},
                                  {150, 300, 450, 600, 750, 900, 1050, 1200}};


  //inline unsigned int nzbins = 13;
  inline unsigned int nzbins = 8;

  inline unsigned int ztozbin(unsigned int layerdisk, double z){
    //double offset = .25;
    double offset = 0;    

    unsigned int zbin = 0;

    for (unsigned int i = 0; i < nzbins; i++){
      if ( z < tmzinner[layerdisk][i]/10 - offset ){
        zbin = i;
        break;
      }
    }
    assert(zbin<99);
    return zbin;
  }

  inline unsigned int nrbins = 8;

  inline double TErbins[8] = {28.0625, 33.625, 39.1875, 44.75, 50.3125, 55.875, 61.4375, 67.0};

  inline unsigned int rtorbin(unsigned int layerdisk, double r, bool isPS){
    if (!(isPS)){
      return 0;
    }

    assert(layerdisk>5);
    unsigned int rbin = 0;
    for (unsigned int i=0; i < nrbins; i++ ){
      if (r <  TErbins[i]){
        rbin = i;
        break;
      }
    }
    return rbin;
  }


  inline  double bendDisk_TE(double r, double z, int disk, double rinv, double stripPitch, bool isPS, bool useTMCorr){
    double dr = 0.18;
    double CF = 1;

    if (useTMCorr){
      if (((disk==1 || disk==2) && isPS) || ((disk==3 || disk==4) && r<=Settings::diskSpacingCut[1]) || (disk==5 && r<=Settings::diskSpacingCut[2])){
        dr = 0.4;
      }
      CF = r/z;
    }

    double delta = r*dr*0.5*rinv*CF;
    double bend  = delta/stripPitch;

    return bend;
  }

  inline double bendBarrel_TE(double r, double z, int layer, double rinv, double stripPitch, bool useTMCorr){

    double dr = 0.18;
    double CF =1;

    if (useTMCorr){
      if ((layer ==1 && z<=Settings::barrelSpacingCut[3]) || (layer==2 && Settings::barrelSpacingCut[1] <=z && z<=Settings::barrelSpacingCut[4]) || (layer==3  && Settings::barrelSpacingCut[3] <=z && z<=Settings::barrelSpacingCut[5])){
        dr = 0.26;
      }
      else if ((layer==1 && Settings::barrelSpacingCut[2]<=z && z<=Settings::barrelSpacingCut[5]) || (layer==2 && Settings::barrelSpacingCut[4]<=z && z<=Settings::barrelSpacingCut[5])){
        dr = 0.4;
      }
      else if ((layer==2 && z<=Settings::barrelSpacingCut[1]) || (layer==3 && z<=Settings::barrelSpacingCut[3])){
        dr =0.16;
      } 
      if ((layer==1 && Settings::barrelSpacingCut[0]<=z && z<=Settings::barrelSpacingCut[5]) || (layer==2 && Settings::barrelSpacingCut[1] <=z && z<=Settings::barrelSpacingCut[5]) || (layer==3 && Settings::barrelSpacingCut[3]<=z && z<=Settings::barrelSpacingCut[5])){
          CF = cosModuleTilt*(z/r) + sinModuleTilt;
      }
    }
    double delta = r*dr*0.5*rinv;
    double bend  = delta/(stripPitch*CF);

    return bend;
}


  inline double bendstrip(double r, double rinv, double stripPitch) {
    constexpr double dr = 0.18;
    double delta = r * dr * 0.5 * rinv;
    double bend = delta / stripPitch;
    return bend;
  }

  inline double rinv(double phi1, double phi2, double r1, double r2) {
    if (r2 <= r1) {  //FIXME can not form tracklet should not call function with r2<=r1
      return 20.0;
    }

    double dphi = phi2 - phi1;
    double dr = r2 - r1;

    return 2.0 * sin(dphi) / dr / sqrt(1.0 + 2 * r1 * r2 * (1.0 - cos(dphi)) / (dr * dr));
  }

  inline std::string convertHexToBin(const std::string& stubwordhex) {
    std::string stubwordbin = "";

    for (char word : stubwordhex) {
      std::string hexword = "";
      if (word == '0')
        hexword = "0000";
      else if (word == '1')
        hexword = "0001";
      else if (word == '2')
        hexword = "0010";
      else if (word == '3')
        hexword = "0011";
      else if (word == '4')
        hexword = "0100";
      else if (word == '5')
        hexword = "0101";
      else if (word == '6')
        hexword = "0110";
      else if (word == '7')
        hexword = "0111";
      else if (word == '8')
        hexword = "1000";
      else if (word == '9')
        hexword = "1001";
      else if (word == 'A')
        hexword = "1010";
      else if (word == 'B')
        hexword = "1011";
      else if (word == 'C')
        hexword = "1100";
      else if (word == 'D')
        hexword = "1101";
      else if (word == 'E')
        hexword = "1110";
      else if (word == 'F')
        hexword = "1111";
      else {
        throw cms::Exception("Inconsistency")
            << __FILE__ << " " << __LINE__ << " hex string format invalid: " << stubwordhex;
      }
      stubwordbin += hexword;
    }
    return stubwordbin;
  }

  inline int ilog2(double factor) {
    double power = log(factor) / log(2);
    int ipower = round(power);
    assert(std::abs(power - ipower) < 0.1);
    return ipower;
  }

  /******************************************************************************
 * Checks to see if a directory exists. Note: This method only checks the
 * existence of the full path AND if path leaf is a dir.
 *
 * @return   1 if dir exists AND is a dir,
 *           0 if dir does not exist OR exists but not a dir,
 *          -1 if an error occurred (errno is also set)
 *****************************************************************************/
  inline int dirExists(const std::string& path) {
    struct stat info;

    int statRC = stat(path.c_str(), &info);
    if (statRC != 0) {
      if (errno == ENOENT) {
        return 0;
      }  // something along the path does not exist
      if (errno == ENOTDIR) {
        return 0;
      }  // something in path prefix is not a dir
      return -1;
    }

    return (info.st_mode & S_IFDIR) ? 1 : 0;
  }

  //Open file - create directory if not existent.
  inline std::ofstream openfile(const std::string& dir, const std::string& fname, const char* file, int line) {
    if (dirExists(dir) != 1) {
      edm::LogVerbatim("Tracklet") << "Creating directory : " << dir;
      int fail = system((std::string("mkdir -p ") + dir).c_str());
      if (fail) {
        throw cms::Exception("BadDir") << file << " " << line << " could not create directory " << dir;
      }
    }

    std::ofstream out(dir + "/" + fname);

    if (out.fail()) {
      throw cms::Exception("BadFile") << file << " " << line << " could not create file " << fname << " in " << dir;
    }

    return out;
  }

  //Open file - create directory if not existent.
  //If first==true open file in create mode, if first==false open in append mode
  inline void openfile(
      std::ofstream& out, bool first, const std::string& dir, const std::string& fname, const char* file, int line) {
    if (dirExists(dir) != 1) {
      edm::LogVerbatim("Tracklet") << "Creating directory : " << dir;
      int fail = system((std::string("mkdir -p ") + dir).c_str());
      if (fail) {
        throw cms::Exception("BadDir") << file << " " << line << " could not create directory " << dir;
      }
    }

    if (first) {
      out.open(fname);
    } else {
      out.open(fname, std::ofstream::app);
    }

    if (out.fail()) {
      throw cms::Exception("BadFile") << file << " " << line << " could not create file " << fname << " in " << dir;
    }
  }

};  // namespace trklet
#endif
