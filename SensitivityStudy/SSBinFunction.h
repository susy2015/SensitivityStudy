#define NTOPJETS_BINS 3
#define NBOTJETS_BINS 3
#define MET_BINS 6
//#define MT2_BINS 3
#define MT2_BINS 2
//4 NTOP NBOT bins with all MET MT2 bins, 5 with only MET: 4*6*3+5*6, 4*6*2+5*6
//#define NSEARCH_BINS 90
#define NSEARCH_BINS 78

const double ntopbins_edge[NTOPJETS_BINS+1] = {1,2,3,4};
const double nbotbins_edge[NBOTJETS_BINS+1] = {1,2,3,4};
const double metbins_edge[MET_BINS+1] = {200.0,250.0,350.0,450.0,550.0,700.0,1000.0};
const double mt2bins_edge[MT2_BINS+1] = {200.0,300.0,500.0};
//const double mt2bins_edge[MT2_BINS+1] = {200.0,300.0,400.0,500.0};

int Set_ntopjetsbin_number(
		                       int ntopjets
							  					)
{   
	int ntopjetsbin_num = -1;

	for(int i=0;i<NTOPJETS_BINS;i++)
  {
    if(i!=NTOPJETS_BINS-1)
    {
      if(ntopjets >= ntopbins_edge[i] && ntopjets < ntopbins_edge[i+1]){ ntopjetsbin_num = i; return ntopjetsbin_num; }
    }
    else
    {
      if(ntopjets >= ntopbins_edge[i]){ ntopjetsbin_num = i; return ntopjetsbin_num; }
    }
  }
  return ntopjetsbin_num;	
}

int Set_nbotjetsbin_number(
		                       int nbotjets
							  			    )
{
	int nbotjetsbin_num = -1;
 
  for(int i=0;i<NBOTJETS_BINS;i++)
  {
    if(i!=NBOTJETS_BINS-1)
    {
      if(nbotjets >= nbotbins_edge[i] && nbotjets < nbotbins_edge[i+1]){ nbotjetsbin_num = i; return nbotjetsbin_num; }
    }
    else
    {
      if(nbotjets >= nbotbins_edge[i]){ nbotjetsbin_num = i; return nbotjetsbin_num; }
    }
  }
  return nbotjetsbin_num;
}

int Set_metbin_number(
                       double met
                     )
{
  int metbin_num = -1;

  for(int i=0;i<MET_BINS;i++)
  {
    if(i!=MET_BINS-1)
    {
      if(met >= metbins_edge[i] && met < metbins_edge[i+1]){ metbin_num = i; return metbin_num; }
    }
    else
    {
      if(met >= metbins_edge[i]){ metbin_num = i; return metbin_num; }
    }
  }
  return metbin_num;
}

int Set_mt2bin_number(
                       double mt2
                     )
{
  int mt2bin_num = -1;

  for(int i=0;i<MT2_BINS;i++)
  {
    if(i!=MT2_BINS-1)
    {
      if(mt2 >= mt2bins_edge[i] && mt2 < mt2bins_edge[i+1]){ mt2bin_num = i; return mt2bin_num; }
    }
    else
    {
      if(mt2 >= mt2bins_edge[i]){ mt2bin_num = i; return mt2bin_num; }
    }
  }
  return mt2bin_num;
}

int Set_SearchBinID( int topbinid, int botbinid, int mt2binid, int metbinid )
{
  int searchbin_num = -1;
  bool GoodInputID = (topbinid>=0 && topbinid<NTOPJETS_BINS) && (botbinid>=0 && botbinid<NBOTJETS_BINS) && (mt2binid>=0 && mt2binid<MT2_BINS) && (metbinid>=0 && metbinid<MET_BINS);
  if( !GoodInputID ) { std::cout << "Bad input ID! What the fuck is going on!!?? Please check SSBinFunction.h!" << std::endl; return -1;}

  int topbotbinid = topbinid*NBOTJETS_BINS + botbinid;//0,1,2,3,4,5,6,7,8

  if     ( topbotbinid == 0 || topbotbinid == 1 ) searchbin_num = topbotbinid*MT2_BINS*MET_BINS + mt2binid*MET_BINS + metbinid;
  else if( topbotbinid == 2 )                     searchbin_num = 2*MT2_BINS*MET_BINS + metbinid;//no MT2 bin structure
  else if( topbotbinid == 3 || topbotbinid == 4 ) searchbin_num = 2*MT2_BINS*MET_BINS + MET_BINS + (topbotbinid-3)*MT2_BINS*MET_BINS + mt2binid*MET_BINS + metbinid;
  else searchbin_num = 4*MT2_BINS*MET_BINS + MET_BINS + (topbotbinid-5)*MET_BINS + metbinid;//no MT2 bin structure

  return searchbin_num;
}
