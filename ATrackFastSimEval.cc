//the script has been slightly modified from https://github.com/eic/EicToyModel/blob/master/fun4all_with_eicroot/sandbox/TrackFastSimEval.cc



#include <iostream>

#include <TVector3.h>
#include <TH1D.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>  

#include <phool/getClass.h>

#include "TrackFastSimEval.h"

using namespace std;

// ---------------------------------------------------------------------------------------

TrackFastSimEval::TrackFastSimEval(const string &name, const string &filename, 
				   const string &trackmapname)
  : SubsysReco(name)
  , _outfile_name(filename)
  , _trackmapname(trackmapname)
  , _event(0)
  , _h1d_Delta_mom(nullptr)
  , h2(nullptr)
  ,h3(nullptr)

{
} // TrackFastSimEval::TrackFastSimEval()

// ---------------------------------------------------------------------------------------

int TrackFastSimEval::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Opening file " << _outfile_name << endl;
  PHTFileServer::get().open(_outfile_name, "RECREATE");

  _h1d_Delta_mom = new TH1D("_h1d_Delta_mom", "#frac{#Delta p}{truth p}", 100, -0.2, 0.2);
 h2 = new TH1D ("eta_dist","#eta distribution",100,0,4);
  h3 = new TH1D ("eta_gen"," #eta disrtibution",100, 0,4);


  return Fun4AllReturnCodes::EVENT_OK;
} // TrackFastSimEval::Init()

// ---------------------------------------------------------------------------------------

int TrackFastSimEval::process_event(PHCompositeNode *topNode)
{
  if (++_event % 100 == 0) cout << PHWHERE << "Events processed: " << _event << endl;

  auto _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  auto _trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname);
  if (!_truth_container || !_trackmap) {
    cout << PHWHERE << " Either PHG4TruthInfoContainer or SvtxTrackMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  } //if

  {
    auto range = _truth_container->GetPrimaryParticleRange();

    // Loop through all the truth particles;
    for (auto truth_itr = range.first; truth_itr != range.second; ++truth_itr) {
      auto g4particle = truth_itr->second;
      if (!g4particle) continue;
	    
	      TVector3 truth_mom(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz());
          TVector3 reco_mom (     track->get_px(),      track->get_py(),      track->get_pz());
          float theta = acos(g4particle->get_pz()/truth_mom.Mag());
          float eta;
          if(theta>0){
            eta = -log(tan(theta/2));}
          else{ eta = -(-log(tan((M_PI + theta)/2)));
          }
          std::cout<<" theta "<<theta<<" eta "<<eta<<std::endl;

          h3->Fill(theta);


      // Loop through all the reconstructed particles and find a match; how about std::map?;
      for (auto track_itr = _trackmap->begin(); track_itr != _trackmap->end(); track_itr++) {
	auto track = dynamic_cast<SvtxTrack_FastSim *>(track_itr->second);
	if (!track) {
	  std::cout << "ERROR CASTING PARTICLE!" << std::endl;
	  continue;
	} //if
	      
	 std::cout<<" eta and id "<<track->get_truth_track_id()<<"   "<<track->get_eta()<<std::endl;
        float track_theta = 2*atan(exp(-track->get_eta()));
        h2->Fill(track_theta);


	// Matching reconstructed particle found, use it and break;
	if ((track->get_truth_track_id() - g4particle->get_track_id()) == 0) {
	  TVector3 truth_mom(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz());
	  TVector3 reco_mom (     track->get_px(),      track->get_py(),      track->get_pz());

	  _h1d_Delta_mom->Fill((reco_mom.Mag() - truth_mom.Mag()) / truth_mom.Mag());

	  break;
	} //if
      } //for 
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
} // TrackFastSimEval::process_event()

// ---------------------------------------------------------------------------------------

int TrackFastSimEval::End(PHCompositeNode *topNode)
{
  PHTFileServer::get().cd(_outfile_name);

  _h1d_Delta_mom->Write();
  h2->Write();
  h3->Write();


  return Fun4AllReturnCodes::EVENT_OK;
} // TrackFastSimEval::End()

// ---------------------------------------------------------------------------------------
