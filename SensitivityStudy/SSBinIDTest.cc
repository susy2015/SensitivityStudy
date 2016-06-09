#include <iostream>
//#include "SSBinFunction.h"
#include "SBGeometry.h"

int main()
{
  std::cout << "Testing Search bin ID algorithm..." << std::endl;
  SBGeometry mySBGeometry;
  if(mySBGeometry.SBSelfTest()){ std::cout << "Good SB Def! NSB: " << mySBGeometry.GetNSB() << "; Continue testing..." << std::endl; mySBGeometry.GetNSB(); }
  else return -1;

  for(int i=0;i<127;i++)
  {
    SBBoundaries outBinDef; mySBGeometry.SBIDToBinBoundaries( i, outBinDef );
    std::cout << i << ": Ntop(" << outBinDef.ntop_lo << "," << outBinDef.ntop_hi << "); Nbot(" << outBinDef.nbot_lo << "," << outBinDef.nbot_hi << "); MT2(" << outBinDef.mt2_lo << "," << outBinDef.mt2_hi << "); MET(" << outBinDef.met_lo << "," << outBinDef.met_hi << ");"<< std::endl;

  }
  /*
  for(int i=0;i<5;i++)
  {
    for(int j=0;j<5;j++)
    {
      for(int k=0;k<20;k++)
      {
        for(int l=0;l<30;l++)
        {
          int ntop = i, nbot = j;
          double mt2 = k*20+200, met = l*20+200;
          int sbid = mySBGeometry.GetSBID(ntop,nbot,mt2,met);
          SBBoundaries outBinDef; mySBGeometry.SBIDToBinBoundaries( sbid, outBinDef );

          if(sbid<0) continue;
          std::cout << "(" << ntop << "," << nbot << "," << mt2 << "," << met << "):" << sbid << std::endl; 
          std::cout << sbid << ": Ntop(" << outBinDef.ntop_lo << "," << outBinDef.ntop_hi << "); Nbot(" << outBinDef.nbot_lo << "," << outBinDef.nbot_hi << "); MT2(" << outBinDef.mt2_lo << "," << outBinDef.mt2_hi << "); MET(" << outBinDef.met_lo << "," << outBinDef.met_hi << ");"<< std::endl;

        }
      }
    }
  }
  */
  /*
  for(int i=-1;i<NTOPJETS_BINS+1;i++)
  {
    for(int j=-1;j<NBOTJETS_BINS+1;j++)
    {
      for(int k=-1;k<MT2_BINS+1;k++)
      {
        for(int l=-1;l<MET_BINS+1;l++)
        {
          int sbid = Set_SearchBinID(i,j,k,l);
          if(sbid<0) continue;
          std::cout << "(" << i << "," << j << "," << k << "," << l << "):" << sbid << std::endl; 
        }
      }
    }
  }
  */
  return 1;
}
