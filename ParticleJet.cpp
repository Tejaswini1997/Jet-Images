#include <stdio.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <math.h>
#define _USE_MATH_DEFINES
 
#include <cmath>
#include <iostream>
 
#include <TFile.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "TCanvas.h"

#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
//#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"

#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
//#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"

#include <sstream>
#include "fastjet/contrib/EnergyCorrelator.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;


#include <fastjet/tools/JHTopTagger.hh>

//.................................................................................
//......... Random number generator .......

#include <chrono>
//main loop begins here
int main()
{
    vector<PseudoJet> inputList;
    vector<PseudoJet> constituents;
    vector<PseudoJet> outputList;
    PseudoJet jettrack, jetphoton, jethadron;
   
    
    
    TFile file("/Users/apple/Delphes-3.5.0/ppjj.root","READ");
    gSystem->Load("libDelphes");
    
    TChain chain("Delphes");
    chain.Add("/Users/apple/Delphes-3.5.0/ppjj.root");
   
    //Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();
    //    TClonesArray *branchphoton = treeReader->UseBranch("photonisolation");//ECAL output
    TClonesArray *branchelectron = treeReader->UseBranch("Electron");//ECAL output
    TClonesArray *branchmuon = treeReader->UseBranch("Muon");//ECAL output
    TClonesArray *branchphoton = treeReader->UseBranch("Photon");
    TClonesArray *branchtower = treeReader->UseBranch("Tower");
    TClonesArray *branchtrack = treeReader->UseBranch("Track");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    
    TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
    TClonesArray *branchGenParticle = treeReader->UseBranch("GenParticle");
    
   
    
    TFile* fileout=new TFile("/Users/apple/Documents/softwares/iyer/ppjjout.root","recreate");
    
    
    TTree *tree = new TTree("TreeB","Background");
    TH1F* h1 = new TH1F("h1", "histogram", 10000, 0.0, 1000);
    
    
    //jet and event variables
    
    
    //Leaves in root browser
    Float_t deltaRljphoton12 =0;
    tree->Branch("deltaRljphoton12", &deltaRljphoton12, "deltaRljphoton12/F");
    Float_t deltaRljphoton22 =0;
    tree->Branch("deltaRljphoton22", &deltaRljphoton22, "deltaRljphoton22/F");
    Float_t deltaRsljphoton12 =0;
    tree->Branch("deltaRsljphoton12", &deltaRsljphoton12, "deltaRsljphoton12/F");
    Float_t deltaRsljphoton22=0;
    tree->Branch("deltaRsljphoton22", &deltaRsljphoton22, "deltaRsljphoton22/F");
    Float_t ljtracksize2=0;
    tree->Branch("ljtracksize2", &ljtracksize2, "ljtracksize2/F");
    Float_t sljtracksize2=0;
    tree->Branch("sljtracksize2", &sljtracksize2, "sljtracksize2/F");
    Float_t deltaR2 =0;
    tree->Branch("deltaR2", &deltaR2, "deltaR2/F");
    Float_t deltaRphotons2 =0;
    tree->Branch("deltaRphotons2", &deltaRphotons2, "deltaRphotons2/F");
    Float_t dRleadgenlj2=0;
    tree->Branch("dRleadgenlj2", &dRleadgenlj2, "dRleadgenlj2/F");
    Float_t dRleadgenslj2=0;
    tree->Branch("dRleadgenslj2", &dRleadgenslj2, "dRleadgenslj2/F");
    Float_t ptw2=0;
    tree->Branch("ptw2", &ptw2, "ptw2/F");
    
    
    
    ofstream variables7;
    variables7.open("/Users/apple/Documents/softwares/iyer/text_files/ppjjout.txt",  std::ios::app);
    
    
   // ofstream variables8;
    //variables8.open("/Users/apple/Documents/softwares/iyer/text_files/zpwwjjWbosonpt.txt",  std::ios::app);
    cout<<numberOfEntries<<"\n";
    
    
    Int_t numjets =0 ;
   
    Float_t electroneta =0;
    Float_t electronphi =0;
    Float_t positroneta =0;
    Float_t positronphi =0;
    Float_t electronpt =0;
    Float_t positronpt =0;
    Float_t leadgenparticleeta =0;
    Float_t leadgenparticlephi =0;
   
    Float_t sumpt =0;
    
    
    fstream jet7_out;
    jet7_out.open("/Users/apple/Documents/softwares/iyer/csv_files/ppjj_100events/jet7_data.csv", ios::out | ios::app); ///
    jet7_out << "deltaetalj" << ","
    << "deltaphilj"<< ","
    << "ljtransverse_energy"<< "\n";

  
    
    
    //event loop begins here
    for(Int_t entry=0; entry < numberOfEntries; ++entry)
    //for(Int_t entry=0; entry < 100; ++entry)
    { treeReader->ReadEntry(entry);
        Float_t photon1eta =0;
        Float_t photon1phi =0;
        Float_t photon2eta =0;
        Float_t photon2phi =0;
       
        Float_t ljconeta = 0;
        Float_t ljconphi = 0;
        Float_t sljconeta = 0;
        Float_t sljconphi = 0;
            Float_t deltaeta =0;
           Float_t  deltaphi =0;
            Float_t deltaR =0;
       Float_t deltaetaelepos =0 ;
        Float_t deltaphielepos=0 ;
        Float_t  deltaRelepos=0;
        Float_t particleindex =0;
        Float_t photonindex =0;
        Float_t deltaRdaughters=0;
        
        Float_t sljtransverse_energy =0;
        Float_t ljtransverse_energy =0;
        
        cout << "\n";
        cout << "================================================================ " <<  "\n";
        cout << entry << "\n";
        
        
        //branchparticle loop begins here
        for (unsigned int jj=0; jj<branchParticle -> GetEntries(); ++jj)
        {
            
            
            //cout << entry << "\t" << "print " << jj <<   "\n";
            GenParticle* particle1 = static_cast<GenParticle*>(branchParticle->At(jj));
            //GenParticle *particle1;
            
           // cout << "first loop"<< "\n";
           
            // get the pid of the particle
            Int_t particlepid = particle1->PID;
            //cout << "particle pid " << particlepid << "\n";
           
            // Get the index of the mother particle
            Int_t motherIndex = particle1->M1;
           
           if (motherIndex >= 0 && motherIndex < branchParticle->GetEntries())
            {
                  GenParticle* motherParticle = static_cast<GenParticle*>(branchParticle->At(motherIndex));

                  // Access the PID of the mother particle
                  Int_t motherPID = motherParticle->PID;

                  if((particlepid == 22) && (motherPID == 9000005))
                  {
                   particleindex = jj;
                     
                      cout << "it is photon" << "\n";
                    //  cout << "second loop"<< "\n";
                      cout << "particle index is " << particleindex << "\n";
                      electroneta = particle1->Eta;
                      electronphi = particle1->Phi;
                      electronpt = particle1->PT;
                      cout << "electron eta and phi " << "\t" << electroneta << "\t" << electronphi << "\n";
                    
                  }
                if((particlepid == -22  )  && (motherPID == 9000005))
                {
                  cout << "it is weird" << "\n";
                    positroneta = particle1->Eta;
                    positronphi = particle1->Phi;
                    positronpt = particle1->PT;
                   // cout << "positron eta and phi " << "\t" << positroneta << "\t" << positronphi << "\n";
                }
               
                deltaetaelepos = electroneta - positroneta;
                deltaetaelepos = electronphi - positronphi;
                deltaRelepos = sqrt(deltaetaelepos*deltaetaelepos + deltaphielepos*deltaphielepos);
                deltaRelepos2 = deltaRelepos;
                
                if(electronpt > positronpt)
                {
                    leadgenparticleeta = electroneta;
                    leadgenparticlephi = electronphi;
                }
                else if (positronpt > electronpt)
                {
                    leadgenparticleeta = positroneta;
                    leadgenparticlephi = positronphi;
                    
                }
                
                
            }
            
          
        } //branchparticle loop ends here  
        
        Int_t photonCount = 0;
        Int_t photonIndex1 = -1;
        Int_t photonIndex2 = -1;
        Double_t photonEta1 = 0.0;
        Double_t photonPhi1 = 0.0;
        Double_t photonEta2 = 0.0;
        Double_t photonPhi2 = 0.0;
        Float_t photonPt1 =0;
        Float_t photonPt2 =0;
        
      
        
    /*    gen_particles = event.UseBranch("Particle")

// Find the W boson
    for particle in gen_particles:
        if particle.PID == 24:  # Assuming W+ boson
            w_boson = particle
            break
    
    // Find the daughter particles of the W boson
    daughters = []
    for particle in gen_particles:
        if particle.Mother1 == w_boson.Index or particle.Mother2 == w_boson.Index:
            daughters.append(particle)
    
    // Calculate delta R between daughter particles
    if len(daughters) == 2:
        delta_eta = daughters[0].Eta - daughters[1].Eta
        delta_phi = ROOT.TVector2.Phi_mpi_pi(daughters[0].Phi - daughters[1].Phi)
        delta_r = ROOT.TMath.Sqrt(delta_eta**2 + delta_phi**2)
        print(f"Delta R between daughter particles: {delta_r}")
        
        */
        

        for (Int_t i = 0; i < branchParticle->GetEntries(); ++i) {
            GenParticle* particle = static_cast<GenParticle*>(branchParticle->At(i));
            
            
            if (particle->PID ==24)
            {   Float_t ptw = particle->PT;
                int daughter1index = particle->D1;
                int daughter2index = particle->D2;
                cout << "checking if the particle is W" << "\n";
                
                
                // Check if both daughters exist
                if (daughter1index >= 0 && daughter2index >= 0) {
                    GenParticle* daughter1 = static_cast<GenParticle*>(branchParticle->At(daughter1index));
                    GenParticle* daughter2 = static_cast<GenParticle*>(branchParticle->At(daughter2index));
                    
                    // Now you have the daughter particles, and you can access their properties
                    double eta1 = daughter1->Eta;
                    double phi1 = daughter1->Phi;
                    
                    double eta2 = daughter2->Eta;
                    double phi2 = daughter2->Phi;
                    
                    cout << "Daughter 1: Eta=" << eta1 << ", Phi=" << phi1 << endl;
                    cout << "Daughter 2: Eta=" << eta2 << ", Phi=" << phi2 << endl;
                    double   deltaetadaughters = eta1 -eta2;
                    double     deltaphidaughters = phi1 - phi2;
                    if(abs(deltaphidaughters) >= pi)
                    {
                        deltaphidaughters = 2*pi - abs(deltaphidaughters);
                    }
                    else if(abs(deltaphidaughters) < pi)
                    {
                        deltaphidaughters = deltaphidaughters;
                        
                    }
                    
                       deltaRdaughters = sqrt(deltaetadaughters*deltaetadaughters + deltaphidaughters*deltaphidaughters);
                    
                    cout << "delta R of daughters " << deltaRdaughters << "\n";
                    
                }
                deltaR2 = deltaRdaughters;
                ptw2=ptw;
               // variables8<<setprecision(3) << "pt of W : " << ptw2 << "\n";
                
            }
        }
        
/*
            for (Int_t i = 0; i < branchParticle->GetEntries(); ++i) {
                GenParticle* particle = static_cast<GenParticle*>(branchParticle->At(i));
            
            // Check if the particle is a photon with the desired mother ID
            if (particle->PID == 22 && particle->M1 >= 0) {
                GenParticle* motherParticle = static_cast<GenParticle*>(branchParticle->At(particle->M1));

                // Check if the mother particle ID matches the desired ID
                if (motherParticle->PID == 9000005) {
                    if (photonCount == 0) {
                        photonIndex1 = i;
                        photonEta1 = particle->Eta;
                        photonPhi1 = particle->Phi;
                        photonPt1 = particle->PT;
                        
                        
                        photonCount++;
                        cout << "photon1 eta and phi " << "\t" << photonEta1 << "\t" << photonPhi1 << "\t"<< photonPt1 << "\n";
                        
                    } else if (photonCount == 1) {
                        photonIndex2 = i;
                        photonEta2 = particle->Eta;
                        photonPhi2 = particle->Phi;
                        photonPt2 = particle->PT;
                        photonCount++;
                        cout << "photon2 eta and phi " << "\t" << photonEta2 << "\t" << photonPhi2 <<"\t" << photonPt2 <<  "\n";
                        break; // Break the loop once we have found both photons
                    }
                }
            }
            
           
        }

        if (photonCount == 2) {
            // Calculate deltaR between the two photons
            Double_t deltaEta = photonEta1 - photonEta2;
            Double_t deltaPhi = std::abs(photonPhi1 - photonPhi2);
            if (deltaPhi > M_PI) deltaPhi = 2 * M_PI - deltaPhi;
            Double_t deltaRphotons = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
            deltaRphotons2 =deltaRphotons;
            // Perform further analysis or print the deltaR value
            std::cout << "DeltaR between the two photons: " << deltaRphotons << std::endl;
        }
        else {
            std::cout << "Did not find two photons with the desired mother ID." << std::endl;
        }
        if(photonPt1>photonPt2)
        {
            leadgenparticleeta = photonEta1;
            leadgenparticlephi = photonPhi1;
        }
        
        else if (photonPt2 > photonPt1)
        {
            leadgenparticleeta = photonEta2;
            leadgenparticlephi = photonPhi2;
            
        }
        cout << "lead gen particle eta and phi" << "\t" << leadgenparticleeta <<"\t"<< leadgenparticlephi << "\n";

        
        */
        
        
        //branchtrack loop begins here
        for(int ii=0; ii < branchtrack -> GetEntriesFast(); ++ii)
        {
            Track *track1=(Track*) branchtrack->At(ii);
            
            Float_t pt1 = track1->PT;
            Float_t eta1 = track1->Eta;
            Float_t phi1 = track1->Phi;
            
            
            //defining 4 momentum variables
            Float_t eT1 = sqrt(pow(pt1,2));
            Float_t pe1 = pt1*cosh(eta1);
            Float_t pz1 = pt1*sinh(eta1);
            Float_t px1 = pt1*cos(phi1);
            Float_t py1 = pt1*sin(phi1);
            
            //         re-scale the track four momentum
            //  if (pt1 >=.1)
            if (pt1 >=2)
            {  //tracks = tracks+1;
                
                Float_t eps = pow(10,5);
                
                jettrack = PseudoJet(px1/eps, py1/eps, pz1/eps, pe1/eps);
                jettrack.set_user_index(0);
                inputList.push_back(jettrack);
            }
        }//branchtrack loop ends here
        
       //branchtower loop begins here
        for(int i=0; i < branchtower -> GetEntriesFast(); ++i)
        {
            Tower *tower2= (Tower*) branchtower->At(i);
            
            Float_t eta1 = tower2->Eta;
            Float_t phi1 = tower2->Phi;
            Float_t  Ehad = tower2->Ehad;
            Float_t Eem = tower2->Eem;
            
            
            //defining 4 momentum variables for ecal
            Float_t eTecal = Eem/cosh(eta1);
            Float_t peecal = Eem;
            Float_t pzecal = eTecal*sinh(eta1);
            Float_t pxecal =  eTecal*cos(phi1);
            Float_t pyecal =  eTecal*sin(phi1);
            
            
            if(eTecal>0.1)
            {
                jetphoton = PseudoJet(pxecal, pyecal, pzecal, peecal);
                jetphoton.set_user_index(1);
                inputList.push_back(jetphoton);
                
                
            }
            //defining 4 momentum variables for hcal
            Float_t eThcal = Ehad/cosh(eta1);
            Float_t pehcal = Ehad;
            Float_t pzhcal = eThcal*sinh(eta1);
            Float_t pxhcal =  eThcal*cos(phi1);
            Float_t pyhcal =  eThcal*sin(phi1);
            Float_t caltowerE = tower2->E;
            
            
            if(eThcal>0.5)
            {
                jethadron = PseudoJet(pxhcal, pyhcal, pzhcal, pehcal);
                jethadron.set_user_index(2);
                inputList.push_back(jethadron);
                
            }
            
           // cout << "cal 4 mom" << pxecal << "\t" << pyecal<< "\t" <<pzecal << "\t" <<
            
            //   Etot = caltowerE ;
            //   thetaJ = log(Ehad/Etot);
            
            //    cout << "thetaj value " << Ehad << "\n";
            
            
    }// branchtower loop ends here
        
     
        //vector pseudojet outputlist;
             double R = 0.4;
                     JetDefinition jet_def1(antikt_algorithm, R);
                    //non area clustering
                     ClusterSequence sequence1(inputList,jet_def1);
             
             outputList.clear();
             outputList = sorted_by_pt(sequence1.inclusive_jets(2.0));
        
        inputList.clear();
        
        
        numjets = outputList.size();
       
        cout << entry << "    number of jets" << "\t" << numjets << "\n";
        
        
        //looping over all the jets in each event
       
        
       
        Int_t ljtracksize =0;
        Int_t sljtracksize =0;
        
        
       
       if(outputList.size()>=2)
       {
           
        //   cout << "outputlist check" << "\n";
          
           Float_t ljeta =0;
           Float_t ljphi =0;
           Float_t sljeta =0;
           Float_t sljphi =0;
           Float_t deltaetaljphoton1 =0;
           Float_t deltaphiljphoton1 =0;
           Float_t deltaetaljphoton2 =0;
           Float_t deltaphiljphoton2 =0;
           Float_t deltaetasljphoton1 =0;
           Float_t deltaphisljphoton1 =0;
           Float_t deltaetasljphoton2 =0;
           Float_t deltaphisljphoton2 =0;
           Float_t deltaRljphoton1 =0;
           Float_t deltaRsljphoton1 =0;
           Float_t deltaRljphoton2 =0;
           Float_t deltaRsljphoton2 =0;
           
           Int_t ljtracksize =0;
           Int_t sljtracksize =0;
           
           Float_t ljconeta =0;
           Float_t ljconphi =0;
           Float_t ljconenergy =0;
           Float_t ljcontransenergy =0;
           Float_t ljdeltaeta =0;
           Float_t ljdeltaphi=0;
           
           Float_t sljconeta =0;
           Float_t sljconphi =0;
           Float_t sljconenergy =0;
           Float_t sljcontransenergy =0;
           Float_t sljdeltaeta =0;
           Float_t sljdeltaphi=0;
           
           Float_t  detaleadgenlj =0;
           Float_t dphileadgenlj =0;
           Float_t detaleadgenslj =0;
           Float_t dphileadgenslj =0;
           Float_t dRleadgenlj =0;
           Float_t dRleadgenslj =0;
          
          
         
           vector<PseudoJet> ljconstituents = outputList[0].constituents();
           vector<PseudoJet> sljconstituents = outputList[1].constituents();
           fstream jet_constituents_out;
      
          
           
           
           
           
           for(unsigned int s =0; s < ljconstituents.size(); s++)
           {
               ljconeta = ljconstituents[s].eta();
               ljconphi = ljconstituents[s].phi();
               ljconenergy = ljconstituents[s].E();
               ljcontransenergy = ljconenergy/cosh(ljconeta);
               
               Float_t ljdeltaeta =0;
               Float_t ljdeltaphi =0;
               Float_t ljdeltaphinew =0;
              
               ljdeltaeta = ljeta - ljconeta;
               ljdeltaphi = ljphi - ljconphi;
               
               
              
               //cout << entry <<"\t" << "slj consttutient size is " << "\t" << sljconstituents.size() << "\n";
               // cout << " slj constituent loop is working" <<"\n";
               if(ljconstituents[s].user_index()==0)
               {
                   ljtracksize = ljtracksize +1;
                  
               }
               ljtracksize2 =ljtracksize;
           }
               // cout << " lj tracksize  is " << ljtracksize << "\n";
           
           for(unsigned int s =0; s < sljconstituents.size(); s++)
           {
               sljconeta = sljconstituents[s].eta();
               sljconphi = sljconstituents[s].phi();
               sljconenergy = sljconstituents[s].E();
               sljcontransenergy = sljconenergy/cosh(sljconeta);
               
               Float_t sljdeltaeta =0;
               Float_t sljdeltaphi =0;
               Float_t sljdeltaphinew =0;
              
               sljdeltaeta = sljeta - sljconeta;
               sljdeltaphi = sljphi - sljconphi;
              
               
               
               //cout << entry <<"\t" << "slj consttutient size is " << "\t" << sljconstituents.size() << "\n";
               // cout << " slj constituent loop is working" <<"\n";
               if(sljconstituents[s].user_index()==0)
               {
                   sljtracksize = sljtracksize +1;
                   
               }
               sljtracksize2 =sljtracksize;
           }
               // cout << " slj tracksize  is " << sljtracksize << "\n";
           
           
           
           //loop for all the jets begins here
           for( unsigned int k =0 ; k < outputList.size(); k++)
           {
              
               // cout << "numjets " << outputList.size() << "\n";
               
               ljeta = outputList[0].eta();
               ljphi = outputList[0].phi();
             
               sljeta = outputList[1].eta();
               sljphi = outputList[1].phi();
               
               deltaetaljphoton1 = ljeta - photonEta1;
               deltaphiljphoton1 = ljphi - photonPhi1;
               if(abs(deltaphiljphoton1) >= pi)
               {
                   deltaphiljphoton1 = 2*pi - abs(deltaphiljphoton1);
               }
               else if(abs(deltaphiljphoton1) < pi)
               {
                   deltaphiljphoton1 = deltaphiljphoton1;
                   
               }
               
               deltaRljphoton1 =sqrt(pow(deltaetaljphoton1,2) + pow(deltaphiljphoton1,2));
               deltaRljphoton12 =deltaRljphoton1;
               
               deltaetaljphoton2 = ljeta - photonEta2;
               deltaphiljphoton2 = ljphi - photonEta2;
               if(abs(deltaphiljphoton2) >= pi)
               {
                   deltaphiljphoton2 = 2*pi - abs(deltaphiljphoton2);
               }
               else if(abs(deltaphiljphoton2) < pi)
               {
                   deltaphiljphoton2 = deltaphiljphoton2;
                   
               }
               deltaRljphoton2 = sqrt(pow(deltaetaljphoton2,2) + pow(deltaphiljphoton2,2));
               deltaRljphoton22 = deltaRljphoton2;
               
               deltaetasljphoton1 = sljeta - photonEta1;
               deltaphisljphoton1 = sljphi - photonPhi1;
               if(abs(deltaphisljphoton1) >= pi)
               {
                   deltaphisljphoton1 = 2*pi - abs(deltaphisljphoton1);
               }
               else if(abs(deltaphisljphoton1) < pi)
               {
                   deltaphisljphoton1 = deltaphisljphoton1;
                   
               }
               deltaRsljphoton1 = sqrt(pow(deltaetasljphoton1,2) + pow(deltaphisljphoton1,2));
               deltaRsljphoton12 = deltaRsljphoton1;
               
               deltaetasljphoton2 = sljeta - photonEta2;
               deltaphisljphoton2 = sljphi - photonPhi2;
               if(abs(deltaphisljphoton2) >= pi)
               {
                   deltaphisljphoton2 = 2*pi - abs(deltaphisljphoton2);
               }
               else if(abs(deltaphisljphoton2) < pi)
               {
                   deltaphisljphoton2 = deltaphisljphoton2;
                   
               }
               
               deltaRsljphoton2 = sqrt(pow(deltaetasljphoton2,2) + pow(deltaphisljphoton2,2));
               deltaRsljphoton22 = deltaRsljphoton2;
               
            //   cout << "delta R of slj and electron" << "\t" << deltaRsljele << "\n";
               
           }  //end of looping over all the jets in each event
           
        /*   cout << "leading jet eta and phi" << "\t"<< ljeta << "\t" << ljphi << "\n";
           cout << "delta eta of lj and electron" << "\t"<< deltaetaljele << "\n";
           cout << "delta phi of lj and electron" << "\t"<< deltaphiljele << "\n";
           cout << "delta eta of lj and positron" << "\t"<< deltaetaljpos << "\n";
           cout << "delta phi of lj and positron" << "\t"<< deltaphiljpos << "\n";
           cout << "delta R of lj and electron" << "\t" << deltaRljele << "\n";
           cout << "delta R of lj and positron" << "\t" << deltaRljpos << "\n";
           cout << "-----------------------------------------------------------------"<<"\n";
           cout << "sub-leading jet eta and phi" << "\t"<< sljeta << "\t" << sljphi << "\n";
           cout << "delta eta of slj and electron" << "\t"<< deltaetasljele << "\n";
           cout << "delta phi of slj and electron" << "\t"<< deltaphisljele << "\n";
           cout << "delta R of slj and electron" << "\t" << deltaRsljele << "\n";
           cout << "delta eta of slj and positron" << "\t"<< deltaetasljpos << "\n";
           cout << "delta phi of slj and positron" << "\t"<< deltaphisljpos << "\n";
           cout << "delta R of slj and positron" << "\t" << deltaRsljpos << "\n";
           cout << "-----------------------------------------------------------------"<<"\n";
          
           */
          
           
         //  if(ljtracksize ==2 && deltaRphotons2 <0.4)
           
          // if(ljtracksize ==2)
         //  cout << "deltaRdaughters values ------------------------------------ " << deltaRdaughters << "\n";
          // if( deltaRdaughters > 0.2)
           //if(deltaRdaughters> 0)
           {   cout << "deltaR check in lj  " << deltaRdaughters<< "\n";
               deltaeta = ljeta - electroneta;
               deltaphi = ljphi -electronphi;
               if(abs(deltaphi) >= pi)
               {
                   deltaphi = 2*pi - abs(deltaphi);
               }
               else if(abs(deltaphi) < pi)
               {
                   deltaphi = deltaphi;
                   
               }
             //  cout << " delta eta and phi  after checking tracks from lj " << "\t"<< deltaeta << "\t" << deltaphi << "\n";
               deltaR = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
              // cout << " delta R after checking tracks from lj " << "\t"  << deltaR << "\n";
               //variables7<<setprecision(3) << deltaR << "\n";
               //if(deltaR < 0.5)
               //{// deltaR2 = deltaR;}
               
               detaleadgenlj = ljeta - leadgenparticleeta;
               dphileadgenlj = ljphi - leadgenparticlephi;
            
               
              if (abs(dphileadgenlj) >= pi)
              {   dphileadgenlj = 2*pi - abs(dphileadgenlj);
                  dRleadgenlj = sqrt(detaleadgenlj*detaleadgenlj + dphileadgenlj * dphileadgenlj);
               dRleadgenlj2 = dRleadgenlj;
              }
               else if(abs(dphileadgenlj) < pi)
               {
                   dRleadgenlj = sqrt(detaleadgenlj*detaleadgenlj + dphileadgenlj * dphileadgenlj);
                dRleadgenlj2 = dRleadgenlj;
               }

              
                  
                  
              
               jet_constituents_out.open("/Users/apple/Documents/softwares/iyer/ppjj_allevents/jet_constituents_data_"+ to_string(entry) +".csv", ios::out | ios::app);
               jet_constituents_out << "deltaetalj" << "," << "deltaphilj" << "," <<"ljtransverse_energy" << "\n" ;
               
               for (int l = 0; l < ljconstituents.size(); l++)
                   
               {
                   //    cout << " constituent " << l << "’s eta: "<< constituents[l].perp() << endl;
                   // cout << "checking 1 " << tracksize << "\n";
                   ljconeta = ljconstituents[l].eta();
                   ljconphi = ljconstituents[l].phi();
                   ljconenergy = ljconstituents[l].E();
                   ljtransverse_energy = ljconenergy/cosh(ljconeta);
                   //  ljconinvmass = ljconstituents[l].m2();
                   //  cout << "coneta and conphi " << coneta << "\t" << conphi <<  "\n" ;
                   
                   
                   Float_t deltaetalj =0;
                   Float_t deltaphilj = 0;
                   Float_t newdeltaphilj =0;
                   deltaetalj = ljeta - ljconeta;
                   deltaphilj = ljphi - ljconphi;
                   Float_t deltaRlj = sqrt(deltaetalj*deltaetalj + deltaphilj*deltaphilj);
                  // cout << "+++++++++++++++++ deltaRlj << " << "\t" << deltaRlj<< "\n";
                   cout << "deltaetalj and deltaphilj " << deltaetalj << "\t" << deltaphilj <<  "\n" ;
                   
                  
                 //  if(deltaRlj < 0.1 && deltaRlj < 0.02)
                 //  if( deltaRlj < 0.1)
                    sumpt = sumpt + ljtransverse_energy;
                       
                       cout << "checking if loop is working" << "\n";
                       if  (abs(deltaphilj) >= pi)
                       {   deltaphilj = 2*pi - abs(deltaphilj);
                           // deltaRsljcon = sqrt(pow(deltaetalj,2) + pow(deltaphilj,2));
                           
                           cout << "delta eta phi and pt of lj" <<  "\t" << deltaetalj << "\t" << deltaphilj << "\t" << ljtransverse_energy<< "\n";
                           variables7<<setprecision(3) << deltaetalj << "\t" << deltaphilj << "\t" << ljtransverse_energy << "\n";
                         //  variables7<<setprecision(3) << ljtransverse_energy << "\n";
                          // jet1_constituents_out << "deltaetalj" << "," << "deltaphilj" << "," <<"ljtransverse_energy" << "\n" ;
                          // jet6_constituents_out.open("jet6_constituents_data_"+ to_string(entry) +".csv", ios::out | ios::app);
                          
                           jet_constituents_out<< deltaetalj << "," << deltaphilj << "," <<ljtransverse_energy << "\n" ;
                            ///
                           }
                           
                       
                       else if(abs(deltaphilj < pi))
                       {  //deltaRljcon = sqrt(pow(deltaetalj,2) + pow(deltaphilj,2));
                           deltaphilj = deltaphilj;
                           cout << "delta eta phi and pt of lj" << "\t"  <<  deltaetalj << "\t" << deltaphilj << "\t" << ljtransverse_energy <<  "\n";
                           variables7<<setprecision(3) << deltaetalj << "\t" << deltaphilj << "\t" << ljtransverse_energy << "\n";
                          // variables7<<setprecision(3) << ljtransverse_energy << "\n";
                          // jet1_constituents_out << "deltaetalj" << "," << "deltaphilj" << "," <<"ljtransverse_energy" << "\n" ;
                          // jet6_constituents_out.open("jet6_constituents_data_"+ to_string(entry) +".csv", ios::out | ios::app);
                           jet_constituents_out << deltaetalj << "," << deltaphilj << "," <<ljtransverse_energy << "\n" ;
                           ///
                           
                         //  jet1_constituents_out << deltaetalj << "," << deltaphilj << "," <<ljtransverse_energy << "\n" ;
                       }
                        
                   }
               }
               
               
           }
         
           //else if(sljtracksize ==2 && deltaRphotons2<0.4)
          // else if(sljtracksize ==2)
        /*   else if(sljtracksize >=0 && deltaRdaughters> 0.1 && deltaRdaughters < 0.4)
           //if(deltaRdaughters> 0)
           {    cout << "deltaR check in slj" << deltaRphotons2 << "\n";
               deltaeta = sljeta -photonEta1;
               deltaphi = sljphi -photonPhi1;
               if(abs(deltaphi) >= pi)
               {
                   deltaphi = 2*pi - abs(deltaphi);
               }
               else if(abs(deltaphi) < pi)
               {
                   deltaphi = deltaphi;
                   
               }
           //    cout << " delta eta and phi after checking tracks from sublj " << "\t" << deltaeta << "\t" << deltaphi << "\n";
               deltaR = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
              // cout << " delta R after checking tracks from slj " << "\t"  << deltaR << "\n";
               //variables7<<setprecision(3) << deltaR << "\n";
              // if(deltaR < 0.5)
              // {// deltaR2 = deltaR;}
               
               detaleadgenslj = sljeta - leadgenparticleeta;
               dphileadgenslj = sljphi - leadgenparticlephi;
               if (abs(dphileadgenslj) >= pi)
               {   dphileadgenslj = 2*pi - abs(dphileadgenslj);
                   dRleadgenslj = sqrt(detaleadgenslj*detaleadgenslj + dphileadgenslj * dphileadgenslj);
                dRleadgenslj2 = dRleadgenslj;
               }
                else if(abs(dphileadgenslj) < pi)
                {
                    dRleadgenslj = sqrt(detaleadgenslj*detaleadgenslj + dphileadgenslj * dphileadgenslj);
                 dRleadgenslj2 = dRleadgenslj;
                }
               
               
               
               for (int l = 0; l < sljconstituents.size(); l++)
                   
               {
               //    cout << " constituent " << l << "’s eta: "<< constituents[l].perp() << endl;
                  // cout << "checking 1 " << tracksize << "\n";
                   sljconeta = sljconstituents[l].eta();
                   sljconphi = sljconstituents[l].phi();
                   sljconenergy = sljconstituents[l].E();
                   sljtransverse_energy = sljconenergy/cosh(sljconeta);
                   
                 //  cout << "coneta and conphi " << coneta << "\t" << conphi <<  "\n" ;
                  
                   
                   Float_t deltaetaslj =0;
                   Float_t deltaphislj = 0;
                   Float_t newdeltaphislj =0;
                   deltaetaslj = sljeta - sljconeta;
                   deltaphislj = sljphi - sljconphi;
                  // cout << "deltaetaslj and deltaphislj " << deltaetaslj << "\t" << deltaphislj <<  "\n" ;
                   
                   if  (abs(deltaphislj) >= pi)
                    {   newdeltaphislj = 2*pi - abs(deltaphislj);
                       // deltaRsljcon = sqrt(pow(deltaetalj,2) + pow(deltaphilj,2));
                        
                        cout << "delta eta phi and pt of slj" <<  "\t" << deltaetaslj << "\t" << newdeltaphislj <<"\t" << sljtransverse_energy<< "\n";
                        variables7<<setprecision(3) << deltaetaslj << "\t" << newdeltaphislj << "\t" << sljtransverse_energy << "\n";
                        
                    }
                  else if(abs(deltaphislj < pi))
                   {  //deltaRljcon = sqrt(pow(deltaetalj,2) + pow(deltaphilj,2));
                       deltaphislj = deltaphislj;
                       cout << "delta eta phi and pt of slj" << "\t"  <<  deltaetaslj << "\t" << deltaphislj << "\t" << sljtransverse_energy <<  "\n";
                   variables7<<setprecision(3) << deltaetaslj << "\t" << deltaphislj << "\t" << sljtransverse_energy << "\n";
                       
                   }
                 
               }
           } */
           
       } //if(outputlistsize >1) ends here
        
        
        tree ->Fill();
        
    } //event loop ends here
    
   // cout << "sum of transverse energy +++++++++++++++++++++++++++++++++++++++++" << sumpt << "\n";
 //  tree->Write();
    
//} //main loop ends here

