#include "HiForestAnalysis/hiForest.h"
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TNtuple.h>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <unordered_set>
#include <sstream> 
// #include <time.h> 

// Example of forest skim

int compare(const void * a, const void * b)
{
  std::pair<int,double> e1 = *(std::pair<int,double>*)a;
  std::pair<int,double> e2 = *(std::pair<int,double>*)b;
  if(e1.second > e2.second) return 1;
  else if (e1.second < e2.second) return -1;
  else return 0;
}

void turnOffBranches(HiForest *c)
{
  c->hasHltTree=0;
  c->hasPFTree=0;
  c->hasPhotonTree=0;
  c->hasTowerTree=0;
  c->hasHbheTree=0;
  c->hasEbTree=0;
  c->hasGenpTree=0;
  c->hasGenParticleTree=0;
  c->hasAk5CaloJetTree=0;
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasAkPu5CaloJetTree=0;
  c->hasAkPu2JetTree=0;
  c->hasAkPu3JetTree=0;
  c->hasAkPu4JetTree=0;
  c->hasAkPu5JetTree=0;
  c->hasAkVs2PFJetTree=0;
  c->hasAkVs3PFJetTree=0;
  c->hasAkVs4PFJetTree=0;
  c->hasAkVs5PFJetTree=0;
}

int findfirstvz(std::pair<int,double> * myevents, double searchvz, double delta)
{
  int arraysize = 986669;
  int index = arraysize/2;
  int granularity = 4;
  for(int i = 0 ; i < 20 ; ++i) // binary search shouldn't take longer than 20 steps
  {
    if( (myevents[index].second < searchvz - delta) ) index = index + arraysize/granularity;
    else if( (myevents[index].second > searchvz - delta) ) index = index - arraysize/granularity;
    granularity *= 2;
  }
  return index;
}

int findlastvz(std::pair<int,double> * myevents, double searchvz, double delta)
{
  int arraysize = 986669;
  int index = arraysize/2;
  int granularity = 4;
  for(int i = 0 ; i < 20 ; ++i) // binary search shouldn't take longer than 20 steps
  {
    // cout<<myevents[index].second<<" vz "<<searchvz + delta<<" at "<<index<<endl;
    if( (myevents[index].second < searchvz + delta) ) index = index + arraysize/granularity;
    else if( (myevents[index].second > searchvz + delta) ) index = index - arraysize/granularity;
    granularity *= 2;
  }
  return index;
}

bool rejectevent( HiForest * c, int centmin, int centmax, double vzrange, double dijetdphicut, double leadingjetpthigh, double leadingjetptlow, double subleadingjetpthigh, double subleadingjetptlow, double jetamax, double jetamin, double ajmin, double ajmax) //all cuts go here
{
  bool willreject = true;
  for (int pseudoindex = 0; pseudoindex < 1; ++pseudoindex)
  {
    
    // ### START Event Selection and quality cuts ###
    if(c->evt.hiBin < centmin || c->evt.hiBin >= centmax) continue;
    if(c->hasGenParticleTree)
    {
      //if is PbPb MC
      if(!c->skim.pcollisionEventSelection) continue;
    }
    else
    {
      //if is PbPb Data
      if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter )) continue;
      // if(!(c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_HIJet80_v1 )) continue;
    }
    if( fabs(c->evt.vz) > vzrange ) continue;
    // ### END Event Selection and quality cuts ###    
    // ### Find leading and subleading jets in |eta| < 2
    int leadindex = -1 , subleadindex = -1;
    bool foundleading = false;
    double oldjetpt = 999999999;
    for(int i = 0 ; i < c->myjet.nref ; ++i)
    {
      double deltajetpt = oldjetpt - c->myjet.jtpt[i];
      if(deltajetpt < 0) { cout<<"FATAL ERROR: unsorted jets"<<endl; exit(1); }
      oldjetpt = c->myjet.jtpt[i];
      if(!foundleading)
      {
        if(fabs(c->myjet.jteta[i])<2.0 )
        {
          leadindex = i;
          foundleading = true;
        }
      }
      else
      {
        if(fabs(c->myjet.jteta[i])<2.0 ) 
        {
          subleadindex = i;
          break;
        }
      }
    }
    if(leadindex < 0 || subleadindex < 0 ) continue; //Didn't find jet pair in |eta|<2
    
    // ### START jet cuts ###
    double dijetdphi = fabs(c->myjet.jtphi[leadindex] - c->myjet.jtphi[subleadindex]);
    if( dijetdphi > TMath::Pi() ) dijetdphi = 2*TMath::Pi() - dijetdphi;
    if( fabs(dijetdphi) < dijetdphicut ) continue;  // dphi cut
    if( c->myjet.jtpt[leadindex] > leadingjetpthigh ) continue;       //lead pt cut
    if( c->myjet.jtpt[leadindex] < leadingjetptlow ) continue;        //lead pt cut
    if( c->myjet.jtpt[subleadindex] > subleadingjetpthigh ) continue; //sublead pt cut
    if( c->myjet.jtpt[subleadindex] < subleadingjetptlow ) continue;  //sublead pt cut
    if( fabs(c->myjet.jteta[leadindex]) > jetamax || fabs(c->myjet.jteta[leadindex]) < jetamin ) continue;        //jet eta cut
    if( fabs(c->myjet.jteta[subleadindex]) > jetamax || fabs(c->myjet.jteta[subleadindex]) < jetamin ) continue;  //jet eta cut
    float aj = (c->myjet.jtpt[leadindex] - c->myjet.jtpt[subleadindex]) / (c->myjet.jtpt[leadindex] + c->myjet.jtpt[subleadindex]);
    if(aj < ajmin || aj > ajmax) continue;  //aj cut
    // ### END jet cuts ###
    willreject = false;
  }
  return willreject;
}

void correlateEventNtuple(char *infname = "0.root", char *mbname = "mb.root", char *outfname = "matched_Jet.root",int parallelization=40,int iteration=0, double deltavz = 0.5)
{
  // srand (time(NULL));
  // ### Some jet cuts ###
  double leadingjetpthigh = 300;
  double leadingjetptlow = 120;
  double subleadingjetpthigh = 300;
  double subleadingjetptlow = 50;
  double jetamax = 1.6;
  double jetamin = 0.0;
  int centmin = 0;
  int centmax = 20;
  int vzrange = 15;
  double dijetdphicut = 5 * TMath::Pi() / 6.0;
  double ajmin = 0.0;
  double ajmax = 1.0;
  // ### End jet cuts  ###
  
  ofstream myfile;
  std::stringstream myfilename; 
  myfilename<<"mixevents_"<<iteration<<".csv";
  myfile.open (myfilename.str());
  // myfile << "Writing this to a file.\n";
  
  cout<<"Loading minbias hiforest... ";
  // Define the input file and HiForest
  HiForest *d = new HiForest(mbname);
  turnOffBranches(d);
  d->SetOutputFile(outfname);
  int nmbentries = d->GetEntries();
  cout<<"done."<<endl;

  cout<<"Loading jet hiforest... ";
  HiForest *c = new HiForest(infname,"forest1",cPPb,0,"akVs3Calo");
  turnOffBranches(c);
  int njetentries = c->GetEntries();
  cout<<"done."<<endl;
  
  std::unordered_set<int> entriestomix;
  std::pair<int,double> * myevents = new std::pair<int,double>[986669];
  std::pair<int,double> * jetevents = new std::pair<int,double>[986669];
  std::vector<int> vjetindices;
  std::vector<int> vmbindices;
  std::vector<int> vmbevt;
  std::vector<int> vmbrun;
  std::vector<int> centralities;
  
  
  int nmb = 0;
  cout<<"Making minbias index,vz pair array... ";
  for (int j = 0 ; j < nmbentries ; ++j)
  {
    d->evtTree->GetEntry(j);
    d->skimTree->GetEntry(j);
    centralities.push_back(d->evt.hiBin);
    vmbevt.push_back(d->evt.evt);
    vmbrun.push_back(d->evt.run);
    if (!(d->skim.pcollisionEventSelection && d->skim.pHBHENoiseFilter )) continue;
    myevents[nmb] = std::make_pair (j,d->evt.vz); //index vz
    vmbindices.push_back(j);
    nmb++;
  }
  cout<<"done."<<endl;
  cout<<"Making minbias index,vz pair array... ";
  int njet = 0;
  for(int i = 0 ; i < njetentries ; ++i)
  {
    c->GetEntry(i);
    if(rejectevent(c,centmin,centmax,vzrange,dijetdphicut,leadingjetpthigh,leadingjetptlow,subleadingjetpthigh,subleadingjetptlow,jetamax,jetamin,ajmin,ajmax)) continue;    
    jetevents[njet] = std::make_pair (c->evt.hiBin,c->evt.vz); //cent vz
    vjetindices.push_back(i);
    njet++;
  }
  cout<<"done."<<endl;

  cout<<"njet = "<<njet<<endl;

  cout<<"Sorting events on vz... ";
  qsort(myevents,986669,sizeof(std::pair<int,double>),compare);
  cout<<"done."<<endl;
  
  TH1D * hnmatched = new TH1D("hnmatched",";n;",1000,0,1000);
  int startindex = (njet/parallelization)*iteration;
  int endindex = (njet/parallelization)*(iteration+1);
  if((njet/parallelization)*(iteration+2) > njet) endindex = njet;
  cout<<"Iterating through jets from "<<startindex<<" to "<<endindex<<" ... ";
  for(int i = startindex ; i < endindex ; ++i)
  {
    if(i%1000==0) 
      cout<<i<<"/"<<njet<<endl;
    int cent = jetevents[i].first;  
    double vz = jetevents[i].second;
    int firstvzinrange = findfirstvz(myevents,vz,deltavz);
    int lastvzinrange = findlastvz(myevents,vz,deltavz);
    int nMatched = 0;
    c->GetEntry(vjetindices[i]);
    std::vector<int> tmp_eventstomix;
    for (int j = firstvzinrange; j <= lastvzinrange; ++j)
    {
      if(cent != centralities[myevents[j].first]) continue;
      // if(c->evt.evt == 17160518) cout<<cent<<" vs "<<centralities[myevents[j].first]<<" vz "<<myevents[j].second<<endl;
      tmp_eventstomix.push_back(myevents[j].first);
    }
    // if(c->evt.evt == 17160518) 
    std::random_shuffle ( tmp_eventstomix.begin(), tmp_eventstomix.end() );
    for (int j = 0 ; j < (int)tmp_eventstomix.size() ; ++j)
    {
      entriestomix.insert(tmp_eventstomix[j]);
      // cout<<c->evt.run<<","<<c->evt.evt<<","<<c->evt.hiBin<<","<<c->evt.vz<<endl;
      
      long long int big = 1000000000;
      long long int thisid = (big*vmbrun[tmp_eventstomix[j] ]) + vmbevt[tmp_eventstomix[j] ];
      myfile<<vjetindices[i]<<","<<thisid<<","<<endl;
      nMatched++;
      if(nMatched>19) break;
    }
    hnmatched->Fill(nMatched);
    // if(i > 2) break;
  }
  cout<<"done."<<endl;
  
  int nfilled = 0;
  cout<<"Iterating minbias events from "<<0<<" to "<<nmbentries<<" ... ";
  for (int j = 0; j < nmbentries; ++j)
  {
    auto search = entriestomix.find(j);
    if(j%10000==0) 
      cout<<j<<"/"<<nmbentries<<" filled: "<<nfilled<<endl;
    if(search == entriestomix.end())
    {
      continue;
    }
    d->GetEntry(j);
    // cout<<d->evt.hiBin<<" "<<d->evt.vz<<endl;
    d->FillOutput();
    nfilled++;
  }
  cout<<"done."<<endl;
  
  delete c;
  delete d;
  
  // HiForest *d = new HiForest(outfname);
  // turnOffBranches(d);
  // d->SetOutputFile("/export/d00/scratch/dav2105/testsort.root");

  myfile.close();
  cout<<"Program completed succesfully."<<endl;
}


int main(int argc, char *argv[])
{
  // if(argc != 7)
  if(argc < 4)
  {
    std::cout << "Usage: ./removeDuplicates.exe <infname> <mbname> <outfname> [parallelization=40] [iteration=0] [deltavz=0.5]" << std::endl;
    return 1;
  }
  correlateEventNtuple(argv[1], argv[2], argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]));
  return 0;
}
