#include <iostream>
#include <vector>

//#define NSB 45
//#define NSB 37
//#define NSB 126
#define NSB 59

#define NTOPJETS_BINS 3
#define NBOTJETS_BINS 3

struct SBBoundaries
{
  double ntop_lo, nbot_lo, mt2_lo, met_lo;
  double ntop_hi, nbot_hi, mt2_hi, met_hi;
};

class SBGeometry
{
 public:
  const static int NTOPBINS = 3;
  const static int NBOTBINS = 3;
  const int NMT2BINS[NTOPBINS*NBOTBINS] = {3,3,2,//top 1
                                           3,3,2,//top 2
                                           1,1,1};//top >=3
  const int NMETBINS[NTOPBINS*NBOTBINS] = {4,4,3,//top 1
                                           4,4,3,//top 2
                                           3,3,3};//top >=3
  const double ntopbins_edge[NTOPBINS+1] = {1,2,3,4};
  const double nbotbins_edge[NBOTBINS+1] = {1,2,3,4};
  const std::vector< std::vector<double> > mt2bins_edge = {
                                                           {200.0,350.0,450.0,600},//top 1 bot 1
                                                           {200.0,350.0,450.0,600},//top 1 bot 2
                                                           {200.0,350.0,500.0},//top 1 bot >=3
                                                           {200.0,350.0,450.0,600},//top 2 bot 1
                                                           {200.0,350.0,450.0,600},//top 2 bot 2
                                                           {200.0,350.0,500.0},//top 2 bot >=3
                                                           {200.0,500.0},//top >=3 bot 1
                                                           {200.0,500.0},//top >=3 bot 2
                                                           {200.0,500.0} //top >=3 bot>=3
                                                          };  
  const std::vector< std::vector<double> > metbins_edge = {
                                                           {200.0,350.0,500.0,650.0,1000.0},//top 1 bot 1
                                                           {200.0,350.0,500.0,650.0,1000.0},//top 1 bot 2
                                                           {200.0,350.0,500.0,800.0},//top 1 bot >=3
                                                           {200.0,350.0,500.0,650.0,1000.0},//top 2 bot 1
                                                           {200.0,350.0,500.0,650.0,1000.0},//top 2 bot 2
                                                           {200.0,350.0,500.0,800.0},//top 1 bot >=3
                                                           {200.0,350.0,500.0,800.0},//top >=3 bot 1
                                                           {200.0,350.0,500.0,800.0},//top >=3 bot 2
                                                           {200.0,350.0,500.0,800.0} //top >=3 bot>=3
                                                          };
  int GetTopID(int ntop);
  int GetBotID(int nbot);
  int GetMT2ID(int topbotid, double mt2);
  int GetMETID(int topbotid, double met);

  bool SBSelfTest();
  int GetNSB();
  int GetSBID(int ntop, int nbot, double mt2, double met);//from 0 to last
  void SBIDToBinBoundaries(int inputIdx, SBBoundaries & outBinDef);
};

int SBGeometry::GetTopID(int ntop)
{  
  int topid = -1;

  for(int i=0;i<NTOPBINS;i++)
  { 
    if(i!=NTOPBINS-1)
    { 
      if(ntop >= ntopbins_edge[i] && ntop < ntopbins_edge[i+1]){ topid = i; return topid; }
    }
    else
    { 
      if(ntop >= ntopbins_edge[i]){ topid = i; return topid; }
    }
  }
  return topid;
}

int SBGeometry::GetBotID(int nbot)
{ 
  int botid = -1;
  
  for(int i=0;i<NBOTBINS;i++)
  { 
    if(i!=NBOTBINS-1)
    { 
      if(nbot >= nbotbins_edge[i] && nbot < nbotbins_edge[i+1]){ botid = i; return botid; }
    }
    else
    { 
      if(nbot >= nbotbins_edge[i]){ botid = i; return botid; }
    }
  }
  return botid;
}

int SBGeometry::GetMT2ID(int topbotid, double mt2)
{
  if(topbotid<0 || topbotid>=NTOPBINS*NBOTBINS) return -2;
  int mt2bin_num = -1;
  
  for(int i=0;i<NMT2BINS[topbotid];i++)
  {
    if(i!=NMT2BINS[topbotid]-1)
    {
      if(mt2 >= (mt2bins_edge.at(topbotid)).at(i) && mt2 < (mt2bins_edge.at(topbotid)).at(i+1)){ mt2bin_num = i; return mt2bin_num; }
    }
    else
    {
      if(mt2 >= (mt2bins_edge.at(topbotid)).at(i)){ mt2bin_num = i; return mt2bin_num; }
    }
  }
  return mt2bin_num;
}

int SBGeometry::GetMETID(int topbotid, double met)
{ 
  if(topbotid<0 || topbotid>=NTOPBINS*NBOTBINS) return -2;
  int metbin_num = -1;
  
  for(int i=0;i<NMETBINS[topbotid];i++)
  { 
    if(i!=NMETBINS[topbotid]-1)
    { 
      if(met >= (metbins_edge.at(topbotid)).at(i) && met < (metbins_edge.at(topbotid)).at(i+1)){ metbin_num = i; return metbin_num; }
    }
    else
    { 
      if(met >= (metbins_edge.at(topbotid)).at(i)){ metbin_num = i; return metbin_num; }
    }
  }
  return metbin_num;
}

bool SBGeometry::SBSelfTest()
{
  bool isgoodSB=true;
  for(int i=0;i<NTOPBINS*NBOTBINS;i++)
  {
    if(NMT2BINS[i]!=(mt2bins_edge.at(i)).size()-1) {isgoodSB = false; std::cout << "Bad SB Def on MT2:" << i << std::endl; return isgoodSB;}
    if(NMETBINS[i]!=(metbins_edge.at(i)).size()-1) {isgoodSB = false; std::cout << "Bad SB Def on MET:" << i << std::endl; return isgoodSB;}
  }
  return isgoodSB;
}

int SBGeometry::GetNSB()
{
  int nsb=0;
  for(int i=0;i<NTOPBINS*NBOTBINS;i++){ nsb+=NMETBINS[i]*NMT2BINS[i]; }
  //std::cout << "Total number of search bin: " << nsb << std::endl;

  return nsb;
}

int SBGeometry::GetSBID(int ntop, int nbot, double mt2, double met)
{
  int sbid = -1;
  int topid = -1, botid = -1, metid = -1, mt2id = -1;
  topid=GetTopID(ntop);
  botid=GetBotID(nbot);
  int topbotid = -1;
  (topid>=0 && botid >=0) ? topbotid = topid*NBOTBINS + botid : topbotid = -1; 
  if(topbotid<0 || topid<0 || botid<0) return -1;

  int sbbase = 0;
  if(topbotid!=0)
  {
    for(int i=0;i<topbotid;i++){ sbbase+=NMT2BINS[i]*NMETBINS[i]; }
  }

  mt2id=GetMT2ID(topbotid, mt2);
  metid=GetMETID(topbotid, met);
  if(mt2id<0 || metid<0) return -1;
  sbid = sbbase + mt2id*NMETBINS[topbotid] + metid;

  //std::cout << topid << "," << botid <<"," << mt2id << ","<< metid << ",SBBase:" << sbbase << std::endl;
  return sbid;
}

void SBGeometry::SBIDToBinBoundaries(int inputIdx, SBBoundaries & outBinDef)
{
  if(inputIdx<0){ std::cout << "bad search bin id!" << std::endl; return ; }

  outBinDef.ntop_lo = -1; outBinDef.ntop_hi = -1;
  outBinDef.nbot_lo = -1; outBinDef.nbot_hi = -1;
  outBinDef.mt2_lo = -1; outBinDef.mt2_hi = -1;
  outBinDef.met_lo = -1; outBinDef.met_hi = -1;
 
  double IDntopnbotedge[NTOPBINS*NBOTBINS+1] = {0};
  for(int i=0;i<NTOPBINS*NBOTBINS+1;i++)
  {
    for(int j=0;j<i;j++)
    {
      IDntopnbotedge[i]+=NMT2BINS[j]*NMETBINS[j];
    }
    //if(i!=0) IDntopnbotedge[i]--;

    if(IDntopnbotedge[i]>inputIdx)
    {
      int topid = (int)(i-1)/NBOTBINS;
      int botid = (i-1)%NBOTBINS;
      outBinDef.ntop_lo = ntopbins_edge[topid];
      outBinDef.nbot_lo = nbotbins_edge[botid];
      topid<NTOPBINS-1 ? outBinDef.ntop_hi = ntopbins_edge[topid+1] : outBinDef.ntop_hi = -1;
      botid<NBOTBINS-1 ? outBinDef.nbot_hi = nbotbins_edge[botid+1] : outBinDef.nbot_hi = -1;

      int mt2metid = inputIdx-IDntopnbotedge[i-1];
      int topbotid = topid*NBOTBINS+botid;
      int mt2id = (int)mt2metid/NMETBINS[topbotid];
      int metid = mt2metid%NMETBINS[topbotid];
      outBinDef.mt2_lo = mt2bins_edge.at(topbotid).at(mt2id);
      outBinDef.met_lo = metbins_edge.at(topbotid).at(metid);

      mt2id<NMT2BINS[topbotid]-1 ? outBinDef.mt2_hi = mt2bins_edge.at(topbotid).at(mt2id+1) : outBinDef.mt2_hi = -1;
      metid<NMETBINS[topbotid]-1 ? outBinDef.met_hi = metbins_edge.at(topbotid).at(metid+1) : outBinDef.met_hi = -1;

      break;
    }
  }
  return ;
}

