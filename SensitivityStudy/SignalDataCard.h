#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <string>

#include "ConstantsSnippet.h"

class SignalDataCard
{
 public:
  std::string DC_SignalType;
  int MMass;
  int DMass;
  double DC_all_MC_Signal_avgweight;
  double DC_sb_MC_Signal[NSB] = {0}, DC_sb_MC_Signal_cs[NSB] = {0};
  double DC_sb_MC_Signal_avgweight[NSB] = {0};
  double DC_sb_MC_Signal_statunc[NSB] = {0}, DC_sb_MC_Signal_systunc[NSB] = {0};
  void print_thisSignalDC();
 private:
  void fake_avg_uncs();
};
