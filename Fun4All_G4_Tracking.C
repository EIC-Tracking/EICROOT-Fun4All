//new version of the Fun4All_G4_Tracking.C

//                                                                                                                                                    
//  EicRoot VST & MM barrel geometry creation and tracking;                                                                                           
//                                                                                                                                                    

#include <phgenfit/Track.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4histos/G4HitNtuple.h>
//#include <gdmlimporter/GdmlImportDetectorSubsystem.h>                                                                                               
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <phool/recoConsts.h>
#include <g4detectors/PHG4GDMLSubsystem.h>

#include <g4tpc/PHG4TpcEndCapSubsystem.h>
#include <g4main/HepMCNodeReader.h>
#include <phhepmc/Fun4AllHepMCInputManager.h>

//#include <gdmlimporter/GdmlImportDetectorSubsystem.h>                                                                                               
//#include <gdmlimporter/SimpleNtuple.h>                                                                                                              
//#include <gdmlimporter/TrackFastSimEval.h>                                                                                                          



#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4Detector.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <phool/recoConsts.h>

#include <EicToyModelSubsystem.h>
#include <EicRootVstSubsystem.h>
#include <EicRootGemSubsystem.h>
#include <EicRootMuMegasSubsystem.h>
#include <EtmOrphans.h>

#include <TrackFastSimEval.h>

R__LOAD_LIBRARY(libeicdetectors.so)
// FIXME: add to CMakeLists.txt;                                                                                                                      
R__LOAD_LIBRARY(libg4trackfastsim.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)

	
R__LOAD_LIBRARY(libg4tpc.so)


// This scrip is simple, sorry: either Qt display or tracking; uncomment if want to see the geometry;                                                 
#define _QT_DISPLAY_
#define _TPC_ //comment out this line for no tpc                                                                                                      
#define _TPC_ENDCAPS_ //comment out this line for no TPC endcaps                                                                                      


namespace G4TPC
{
  int n_tpc_layer_inner = 16;
  double tpc_layer_thick_inner = 1.25;  // EIC- recover default inner radius of TPC vol.                                                              
  int tpc_layer_rphi_count_inner = 1152;

  int n_tpc_layer_mid = 16;
  double tpc_layer_thick_mid = 1.25;
  int tpc_layer_rphi_count_mid = 1536;

  int n_tpc_layer_outer = 16;
  double tpc_layer_thick_outer = 1.125;  // outer later reach from 60-78 cm (instead of 80 cm), that leads to radial thickness of 1.125 cm            
  int tpc_layer_rphi_count_outer = 2304;

  double outer_radius = 78.;
  double inner_cage_radius = 20.;
  double cage_length = 211.0;  // From TPC group, gives eta = 1.1 at 78 cm                                                                            
  double n_rad_length_cage = 1.13e-02;
  double cage_thickness = 28.6 * n_rad_length_cage;  // Kapton X_0 = 28.6 cm  // mocks up Kapton + carbon fiber structure                             

  string tpcgas = "sPHENIX_TPC_Gas";  //  Ne(90%) CF4(10%) - defined in g4main/PHG4Reco.cc                                                            

}  // namespace G4TPC                                                                                                                                 



void Fun4All_G4_Tracking(int nEvents = 1000)
{
  // The fun4all Server;                                                                                                                              
  Fun4AllServer *se = Fun4AllServer::instance();

  // May want to fix the random seed for reproducibility;                                                                                             
  recoConsts::instance()->set_IntFlag("RANDOMSEED", 12345);

  // Particle Generator Setup;                                                                                                                        
  {
    auto gen = new PHG4ParticleGenerator();

    gen->set_name("pi+");                 // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ...     
    gen->set_vtx(0,0,0);                  // Vertex generation range                                                                                  
    gen->set_mom_range(10., 10.);         // Momentum generation range in GeV/c                                                                       
    gen->set_eta_range(-1.5,1.5);                 // Detector coverage around theta ~ 90 degrees in this example                                      
    gen->set_phi_range(0.,2.*TMath::Pi());
    se->registerSubsystem(gen);
  }

  // fun4all Geant4 wrapper;                                                                                                                          
  PHG4Reco* g4Reco = new PHG4Reco();

  // BeAST magnetic field;                                                                                                                            
  g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/mfield.4col.dat"), PHFieldConfig::kFieldBeast);

  // EicRoot media import; neither bound to EicToyModel nor to a particular EicRoot detector;                                                         
  EicGeoParData::ImportMediaFile("../../examples/eicroot/media.geo");

  // EicRoot vertex tracker; be aware: "VST" will also become a SuperDetector name;                                                                   
  auto vst = new EicRootVstSubsystem("VST");
  {
    vst->SetGeometryType(EicGeoParData::NoStructure);
    vst->SetActive(true);

    // Barrel layers; hits belonging to these layers will be labeled internally                                                                       
    // according to the sequence of these calls;                                                                                                      
    {
      auto ibcell = new MapsMimosaAssembly();
      // See other MapsMimosaAssembly class POD entries in MapsMimosaAssembly.h;                                                                      
      ibcell->SetDoubleVariable("mAssemblyBaseWidth", 17.5 * etm::mm);

      // Compose barrel layers; parameters are:                                                                                                       
      //  - cell assembly type;                                
      //  - number of staves in this layer;                                                                                                           
      //  - number of chips in a stave;                                                                                                               
      //  - chip center installation radius;                                                                                                          
      //  - additional stave slope around beam line direction; [degree];                                                                              
      //  - layer rotation around beam axis "as a whole"; [degree];                                                                                   
      vst->AddBarrelLayer(ibcell, 1*3*12,  1*9, 1*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 2*3*12,  1*9, 2*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 3*3*12,  2*9, 3*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 4*3*12,  2*9, 4*3*23.4 * etm::mm, 12.0, 0.0);
    }

    g4Reco->registerSubsystem(vst);
  }
  /*                                                                                                                                                  
  // EicRoot micromegas central tracker barrels;                                                                                                      
  auto mmt = new EicRootMuMegasSubsystem("MMT");                                                                                                      
  {                                                                                                                                                   
    mmt->SetActive(true);                                                                                                                             
                                                                                                                                                      
    {                                                                                                                                                 
      auto layer = new MuMegasLayer();                                                                                                                
      // See other MuMegasLayer class POD entries in MuMegasGeoParData.h;                                                                             
      layer->SetDoubleVariable("mOuterFrameWidth", 20 * etm::mm);                                                                                     
                                                                                                                                                      
      // Compose barrel layers; parameters are:                                                                                                       
      //   - layer description (obviously can mix different geometries);                                                                              
      //   - length along Z;                                                                                                                          
      //   - segmentation in Z;                                                                                                                       
      //   - radius;                                                                                                                                  
      //   - segmentation in phi;                                                                                                                     
      //   - Z offset from 0.0 (default);  
       //   - azimuthal rotation from 0.0 (default);                                                                                                   
      mmt->AddBarrel(layer, 600 * etm::mm, 2, 300 * etm::mm, 3, 0.0, 0.0);                                                                            
      mmt->AddBarrel(layer, 600 * etm::mm, 3, 400 * etm::mm, 4, 0.0, 0.0);                                                                            
      //      mmt->AddBarrel(layer, 2400 * etm::mm, 4, 600 * etm::mm, 4, 0.0, 0.0);                                                                   
      mmt->SetTransparency(50);                                                                                                                       
    }                                                                                                                                                 
                                                                                                                                                      
    g4Reco->registerSubsystem(mmt);                                                                                                                   
  }                                                                                                                                                   
  */

#ifdef _TPC_
  // time projection chamber layers --------------------------------------------                                                                      
  //This is just cylinder subsystems basically.                                                                                                       
  //What layer numbers to use? This will be the detector (group) number, right? So.. We can start at 60.                                              
  PHG4CylinderSubsystem* cylTPC;
  cylTPC->CheckOverlap();

  double radius = G4TPC::inner_cage_radius;

  // inner field cage                                                                                                                                 
  cylTPC = new PHG4CylinderSubsystem("FIELDCAGE", 60);
  cylTPC->set_double_param("radius", radius);
  cylTPC->set_double_param("length", G4TPC::cage_length);
  cylTPC->set_string_param("material", "G4_KAPTON");
  cylTPC->set_double_param("thickness", G4TPC::cage_thickness);
  cylTPC->SuperDetector("FIELDCAGE");
  cylTPC->Verbosity(0);
  g4Reco->registerSubsystem(cylTPC);


  int n_tpc_layers[3] = {16, 16, 16};
  int tpc_layer_rphi_count[3] = {1152, 1536, 2304};

  double tpc_region_thickness[3] = {20., 20., 18.};
  // Active layers of the TPC (inner layers)                                                                                                          
  int nlayer = 60;
  for (int irange = 0; irange < 3; irange++)
    {
      double tpc_layer_thickness = tpc_region_thickness[irange] / n_tpc_layers[irange];  // thickness per layer                                       
      for (int ilayer = 0; ilayer < n_tpc_layers[irange]; ilayer++)
        {
          cylTPC = new PHG4CylinderSubsystem("TPC", nlayer);
          cylTPC->set_double_param("radius", radius);
          cylTPC->set_double_param("length", G4TPC::cage_length);
          cylTPC->set_string_param("material", G4TPC::tpcgas);
          cylTPC->set_double_param("thickness", tpc_layer_thickness - 0.01);
          cylTPC->SetActive();
          cylTPC->SuperDetector("TPC");
          g4Reco->registerSubsystem(cylTPC);

          radius += tpc_layer_thickness;
          nlayer++;
        }
    }

  // outer field cage                                                                                                                                 
  cylTPC = new PHG4CylinderSubsystem("FIELDCAGE", nlayer);
  cylTPC->set_double_param("radius", radius);
  cylTPC->set_int_param("lengthviarapidity", 0);
  cylTPC->set_double_param("length", G4TPC::cage_length);
  cylTPC->set_string_param("material", "G4_KAPTON");
  cylTPC->set_double_param("thickness", G4TPC::cage_thickness);  // Kapton X_0 = 28.6 cm                                                              
  cylTPC->SuperDetector("FIELDCAGE");
    g4Reco->registerSubsystem(cylTPC);
#endif //End of TPC definition                                                                                                                        

#ifdef _TPC_ENDCAPS_
  bool AbsorberActive = false;

  PHG4TpcEndCapSubsystem* tpc_endcap = new PHG4TpcEndCapSubsystem("TPC_ENDCAP");
  tpc_endcap->SuperDetector("TPC_ENDCAP");

  if (AbsorberActive) {
    tpc_endcap->SetActive();
  }
  tpc_endcap->CheckOverlap();

  //  tpc_endcap->set_int_param("construction_verbosity", 2);                                                                                         

  g4Reco->registerSubsystem(tpc_endcap);
  #endif





  // Forward GEM tracker module(s);                                                                                                                   
  auto fgt = new EicRootGemSubsystem("FGT");
  {
    fgt->SetActive(true);
    fgt->CheckOverlap();
    //fgt->SetTGeoGeometryCheckPrecision(0.000001 * etm::um);                                                                                         

    {
      auto sbs = new GemModule();
  // See other GemModule class data in GemGeoParData.h;                                                                                           
      sbs->SetDoubleVariable("mFrameBottomEdgeWidth", 30 * etm::mm);

      // Compose sectors; parameters are:                                                                                                             
      //   - layer description (obviously can mix different geometries);                                                                              
      //   - azimuthal segmentation;                                                                                                                  
      //   - gas volume center radius;                                                                                                                
      //   - Z offset from 0.0 (default);                                                                                                             
      //   - azimuthal rotation from 0.0 (default);                                                                                                   
      fgt->AddWheel(sbs, 12, 420 * etm::mm, 1200 * etm::mm, 0);
      fgt->AddWheel(sbs, 12, 420 * etm::mm, 1300 * etm::mm, 0);
    }

    g4Reco->registerSubsystem(fgt);
  }

  // Truth information;                                                                                                                               
  g4Reco->registerSubsystem(new PHG4TruthSubsystem());
  se->registerSubsystem(g4Reco);

#ifdef _QT_DISPLAY_
  g4Reco->InitRun(se->topNode());
  g4Reco->ApplyDisplayAction();
  g4Reco->StartGui();
#else
  // Ideal track finder and Kalman filter;                                                                                                            
  {
 auto kalman = new PHG4TrackFastSim("PHG4TrackFastSim");

    kalman->set_use_vertex_in_fitting(false);

    // Silicon tracker hits;                                                                                                                          
    kalman->add_phg4hits(vst->GetG4HitName(),           // const std::string& phg4hitsNames                                                           
                         PHG4TrackFastSim::Cylinder,    // const DETECTOR_TYPE phg4dettype                                                            
                         999.,                          // radial-resolution [cm] (this number is not used in cylindrical geometry)                   
                         // 20e-4/sqrt(12) cm = 5.8e-4 cm, to emulate 20x20 um pixels;                                                                
                         5.8e-4,                        // azimuthal (arc-length) resolution [cm]                                                     
                         5.8e-4,                        // longitudinal (z) resolution [cm]                                                           
                         1,                             // efficiency (fraction)                                                                      
                         0);                            // hit noise                                                                                  
    /*                                                                                                                                                
    // MM tracker hits;                                                                                                                               
    kalman->add_phg4hits(mmt->GetG4HitName(),           // const std::string& phg4hitsNames                                                           
                         PHG4TrackFastSim::Cylinder,    // const DETECTOR_TYPE phg4dettype                                                            
                         999.,                          // radial-resolution [cm] (this number is not used in cylindrical geometry)                   
                         // Say 50um resolution?; [cm];                                                                                               
                         50e-4,                         // azimuthal (arc-length) resolution [cm]                                                     
                         50e-4,                         // longitudinal (z) resolution [cm]                                                           
                         1,                             // efficiency (fraction)                                                                      
                         0);                            // hit noise                                                                                  
    */

    //TPC simpler, since it's fundamental cylinders                                                                                                   
    //No need to here keep track of the layer numbers                                                                                                 
    #ifdef _TPC_
    kalman->add_phg4hits(
                         "G4HIT_TPC",                //      const std::string& phg4hitsNames,                                                        
                         PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,      
	                 1,                           //      const float radres,                                                                     
                         200e-4,                      //      const float phires,                                                                     
                         500e-4,                      //      const float lonres,                                                                     
                         1,                           //      const float eff,                                                                        
                         0                            //      const float noise                                                                       
                         );
    #endif





    // GEM tracker hits; should work;                                                                                                                 
    kalman->add_phg4hits(fgt->GetG4HitName(),           // const std::string& phg4hitsNames                                                           
                         PHG4TrackFastSim::Vertical_Plane,      // const DETECTOR_TYPE phg4dettype                                                    
                         999.,                          // radial-resolution [cm] (this number is not used in cylindrical geometry)                   
                         // Say 70um resolution?; [cm];                                                                                               
                         70e-4,                         // azimuthal (arc-length) resolution [cm]                                                     
                         70e-4,                         // longitudinal (z) resolution [cm]                                                           
                         1,                             // efficiency (fraction)                                                                      
                         0);                            // hit noise                                                                                  

    se->registerSubsystem(kalman);
  }



  // User analysis code: just a single dp/p histogram;                                                                                                
  se->registerSubsystem(new TrackFastSimEval());
  /* SimpleNtuple * hits = new SimpleNtuple("Hits");                                                                                                  
#ifdef _TPC_                                                                                                                                          
           hits->AddNode("TPC", 60);                                                                                                                  
  #endif                                                                                                                                              
                                                                                                                                                      
           se->registerSubsystem(hits);                                                                                                               
  */
  // Run it all, eventually;                                                                                                                          
  se->run(nEvents);
  se->End();
#endif
  delete se;
  gSystem->Exit(0);
} // Fun4All_G4_Sandbox()                                                                                                                             



