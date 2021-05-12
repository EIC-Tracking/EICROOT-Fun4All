#pragma once
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
#include <gdmlimporter/GdmlImportDetectorSubsystem.h>
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

#include <gdmlimporter/GdmlImportDetectorSubsystem.h>
#include <gdmlimporter/SimpleNtuple.h>
#include <gdmlimporter/TrackFastSimEval.h>

//to include gems                                                                                                                                     
#include <EicToyModelSubsystem.h>
#include <EicRootGemSubsystem.h>
#include <EtmOrphans.h>



//For MPGD from QH                                                                                                                                    
#include <g4exampledetector/PHG4CylinderStripSubsystem.h>
#include <g4exampledetector/CreateCZHitContainer.h>
#include "detector_setup.h"

R__LOAD_LIBRARY(libeicdetectors.so)

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libgdmlimportdetector.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)
R__LOAD_LIBRARY(libg4tpc.so)

R__LOAD_LIBRARY(libg4example01detector.so) //For MPGD from QH                                                                                         

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




void Fun4All_G4_HybridBaseline(
                        int nEvents = 1000,
                        const char *outputFile = "out_hybridBaseline",
                        const char *genpar = "e-")
{
 // ======================================================================================================                                     
        // Input from the user                                                                                                                        
        const int particle_gen = 1;     // 1 = particle generator, 2 = particle gun, 3 = nothing (i.e. HepMC input, for example)                      
        const int magnetic_field = 1;   // 1 = uniform 1.5T, 2 = uniform 3.0T, 3 = sPHENIX 1.4T map, 4 = Beast 3.0T map                               
        // ======================================================================================================                                     
        // Parameters for black holes                                                                                                                 
        string radialBH   = "DIRC";             // Cylindrical surface object name                                                                    
        double radialBHradius = 80.0;           // [cm]                                                                                               
        double radialBHlength     = 400.0;              // [cm]                                                                                       
        // ---                                                                                                                                        
        double thickness    = 0.1;              // black hole thickness, needs to be taken into account for the z positions                           
        // ---                                                                                                                                        
        string forwardBH   = "FOR";             // Forward plane object name                                                                          
        double forwardBHz   = 130.0+thickness/2.0;      // [cm]                                                                                       
        double forwardBHradius = 80.0;          // [cm]                                                                                               
        // ---                                                                                                                                        
        string backwardBH   = "BACK";           // Backward plane object name                                                                         
        double backwardBHz   = -(130.0+thickness/2.0);// [cm]                                        
        double backwardBHradius = 80.0;         // [cm]                                                                                               
        // ======================================================================================================                                     
        // Make the Server                                                                                                                            
        Fun4AllServer *se = Fun4AllServer::instance();
        se->Verbosity(1);
        // If you want to fix the random seed for reproducibility                                                                                     
        recoConsts *rc = recoConsts::instance();
        rc->set_IntFlag("RANDOMSEED", 12345);
        // ======================================================================================================                                     
        // Particle Generation                                                                                                                        
        cout << "Particle that will be generated: " << std::string(genpar) << endl;
        // --------------------------------------------------------------------------------------                                                     
        // Particle Generator Setup                                                                                                                   
        PHG4ParticleGenerator * gen = new PHG4ParticleGenerator();
        gen->set_name(std::string(genpar));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)                  
        gen->set_vtx(0,0,0);                    // Vertex generation range                                                                            
        gen->set_mom_range(0.0, 30.0);          // Momentum generation range in GeV/c                                                                 
        gen->set_z_range(0.0, 0.0); //THIS NEEDS TO BE SET, TO NOT DEFAULT TO -10 TO 10                                                               
        gen->set_eta_range(-4.0, 4.0);          // Detector coverage corresponds to |η|< 4                                                            
        gen->set_phi_range(0.0, 2.0*TMath::Pi());

        // --------------------------------------------------------------------------------------                                                     
        // Particle Gun Setup                                                                                                                         
        PHG4ParticleGun *gun = new PHG4ParticleGun();
        gun->set_name(std::string(genpar));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ...                                                 
        gun->set_vtx(0,0,0);
        gun->set_mom(0,1,0);
        // --------------------------------------------------------------------------------------                                                     
        if     (particle_gen==1){
                se->registerSubsystem(gen); cout << "Using particle generator" << endl;
        } else if(particle_gen==2){
                     se->registerSubsystem(gun); cout << "Using particle gun"       << endl;
        }


        // ======================================================================================================                                     
        PHG4Reco *g4Reco = new PHG4Reco();

        EicGeoParData::ImportMediaFile("../../examples/eicroot/media.geo");

        // ======================================================================================================                                     
        // Magnetic field setting                                                                                                                     
        TString B_label;
        if(magnetic_field==1){          // uniform 1.5T                                                                                               
                B_label = "_B_1.5T";
                g4Reco->set_field(1.5);
        }
        else if(magnetic_field==2){     // uniform 3.0T                                                                                               
                B_label = "_B_3.0T";
                g4Reco->set_field(3.0);
        }
        else if(magnetic_field==3){     // sPHENIX 1.4T map                                                                                           
                B_label = "_sPHENIX";
                g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"), PHFieldConfig::kField2D);
                g4Reco->set_field_rescale(-1.4/1.5);
        }
        else if(magnetic_field==4){     // Beast 3.0T map                                                                                             
                B_label = "_Beast";
                g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/mfield.4col.dat"), PHFieldConfig::kFieldBeast);
        }
        else{                           // The user did not provide a valid B field setting                                                           
                cout << "User did not provide a valid magnetic field setting. Set 'magnetic_field'. Bailing out!" << endl;
        }
 #ifdef _BEAMPIPE_
        // Central pipe dimension                                                                                                                     
        // Extracted via mechanical model: Detector chamber 3-20-20                                                                                   
        // directly implimenting the central Be section in G4 cylinder for max speed simulation in the detector region.                               
        // The jointer lip structure of the pipe R = 3.2cm x L=5mm is ignored here                                                                    
        const double be_pipe_radius = 3.1000;
        const double be_pipe_thickness = 3.1762 - be_pipe_radius;  // 760 um for sPHENIX                                                              
        const double be_pipe_length_plus = 66.8;                   // +z beam pipe extend.                                                            
        const double be_pipe_length_neg = -79.8;                   // -z beam pipe extend.                                                            
        bool overlapcheck = false;
        bool absorberactive = false;

        // process pipe extentions?                                                                                                                   
        bool use_forward_pipes = false;
        #ifdef _USE_FORWARD_PIPES_
        use_forward_pipes = true;
        #endif
        const bool do_pipe_hadron_forward_extension = use_forward_pipes && true;
        const bool do_pipe_electron_forward_extension = use_forward_pipes && true;

        const double be_pipe_length = be_pipe_length_plus - be_pipe_length_neg;  // pipe length                                                       
        const double be_pipe_center = 0.5 * (be_pipe_length_plus + be_pipe_length_neg);

        // mid-rapidity beryillium pipe                                                                                                               
        PHG4CylinderSubsystem* cyl = new PHG4CylinderSubsystem("VAC_BE_PIPE", 0);
        cyl->set_double_param("radius", 0.0);
        cyl->set_int_param("lengthviarapidity", 0);
        cyl->set_double_param("length", be_pipe_length);
        cyl->set_double_param("place_z", be_pipe_center);
        cyl->set_string_param("material", "G4_Galactic");
 cyl->set_double_param("thickness", be_pipe_radius);
        cyl->SuperDetector("PIPE");
        cyl->OverlapCheck(overlapcheck);
        //cyl->set_color(0, 0, 0, 0); //To make vacuum invisible                                                                                      
        if (absorberactive) cyl->SetActive();
        g4Reco->registerSubsystem(cyl);

        cyl = new PHG4CylinderSubsystem("BE_PIPE", 1);
        cyl->set_double_param("radius", be_pipe_radius);
        cyl->set_int_param("lengthviarapidity", 0);
        cyl->set_double_param("length", be_pipe_length);
        cyl->set_double_param("place_z", be_pipe_center);
        cyl->set_string_param("material", "G4_Be");
        cyl->set_double_param("thickness", be_pipe_thickness);
        cyl->SuperDetector("PIPE");
        cyl->OverlapCheck(overlapcheck);
        //cyl->set_color(192.0/255.0, 192.0/255.0, 192.0/255.0, 1); //To give beryllium greyish colour                                                
        if (absorberactive) cyl->SetActive();
        g4Reco->registerSubsystem(cyl);

        if (do_pipe_electron_forward_extension)
        {
        PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("ElectronForwardEnvelope");
        gdml->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/Beam/Detector chamber 3-20-20.G4Import.gdml");
        gdml->set_string_param("TopVolName", "ElectronForwardEnvelope");
        gdml->set_int_param("skip_DST_geometry_export", 1);  // do not export extended beam pipe as it is not supported by TGeo and outside Kalman filter acceptance                                                                                                                                       
        gdml->OverlapCheck(1);
        g4Reco->registerSubsystem(gdml);
        }
if (do_pipe_hadron_forward_extension)
        {
        PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("HadronForwardEnvelope");
        gdml->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/Beam/Detector chamber 3-20-20.G4Import.gdml");
        gdml->set_string_param("TopVolName", "HadronForwardEnvelope");
        gdml->set_int_param("skip_DST_geometry_export", 1);  // do not export extended beam pipe as it is not supported by TGeo and outside Kalman fi\
lter acceptance                                                                                                                                       
        gdml->OverlapCheck(1);
        g4Reco->registerSubsystem(gdml);
        }
        #endif


        #ifdef _BARREL_
	// Loading barrel from gdml file                                                                                                              
        GdmlImportDetectorSubsystem * svtPart = new GdmlImportDetectorSubsystem();
        svtPart->set_string_param("GDMPath","Vst_GDML_.gdml");
        //svtPart->OverlapCheck(); //Doesn't do anything, it seems.                                                                                   
        svtPart->AddAssemblyVolume("VST");      // Barrel. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********                        
        svtPart->SuperDetector("SVT");
        svtPart->SetActive();          // this saves hits in the MimosaCore volumes                                                                   
        svtPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)                                                       
        g4Reco->registerSubsystem(svtPart);
        #endif

        #ifdef _FORWARD_SILICON_DISKS_
        //Forward disks                                                                                                                               
        GdmlImportDetectorSubsystem * fstPart = new GdmlImportDetectorSubsystem();
        fstPart->set_string_param("GDMPath","Fst_GDML_.gdml");
        //fstPart->OverlapCheck(); //Doesn't do anything, it seems.        
   fstPart->AddAssemblyVolume("FST");      // Disks. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********                         
        fstPart->SuperDetector("FstDisks");
        fstPart->SetActive();          // this saves hits in the MimosaCore volumes                                                                   
        fstPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)                                                       
        g4Reco->registerSubsystem(fstPart);
        #endif

        #ifdef _BACKWARD_SILICON_DISKS_
        //Backward disks                                                                                                                              
        GdmlImportDetectorSubsystem * bstPart = new GdmlImportDetectorSubsystem();
        bstPart->set_string_param("GDMPath","Bst_GDML_.gdml");
        //bstPart->OverlapCheck(); //Doesn't do anything, it seems.                                                                                   
        bstPart->AddAssemblyVolume("BST");      // Disks. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********                         
        bstPart->SuperDetector("BstDisks");
        bstPart->SetActive();          // this saves hits in the MimosaCore volumes                                                                   
        bstPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)                                                       
        g4Reco->registerSubsystem(bstPart);
        #endif

        #ifdef _TPC_
        // time projection chamber layers --------------------------------------------                                                                
        //This is just cylinder subsystems basically.                                                                                                 
        //What layer numbers to use? This will be the detector (group) number, right? So.. We can start at 60.                                        
        PHG4CylinderSubsystem* cylTPC;

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
        //tpc_endcap->OverlapCheck(OverlapCheck);                                                                                                     

        //  tpc_endcap->set_int_param("construction_verbosity", 2);                                                                                   

        g4Reco->registerSubsystem(tpc_endcap);
        #endif

        #ifdef _MPGD_
        double gap_betweenCZ = 1.5, Gap_betweenlayer = 1.5;
        //double thickness = 0.355199;                                                                                                                
        double thicknessMPGD = 0.36499;
 int nCZlayer = 2;
        bool use_2Dreadout = true;
        if (use_2Dreadout) {
                gap_betweenCZ = 0;
                nCZlayer = 1;
        }
  //double BMT_r[6] = {20., 20.+nCZlayer*thicknessMPGD+gap_betweenCZ+Gap_betweenlayer, 50-nCZlayer*thicknessMPGD-gap_betweenCZ-Gap_betweenlayer/2, 50\
+Gap_betweenlayer/2, 80-(nCZlayer*thicknessMPGD+gap_betweenCZ)*2-Gap_betweenlayer, 80-nCZlayer*thicknessMPGD-gap_betweenCZ};                          
  double BMT_r[6] = {
    47.7153,
    49.5718,
    71.8958,
    73.7523,
    75.6088,
    77.4653
  };

        PHG4CylinderStripSubsystem *example01;
        const double prapidity =1;
        double bmt_length = (1-exp(-2*prapidity))/exp(-prapidity)*80;
        for (int ilayer = 0; ilayer< 6; ilayer++){
                example01 = new PHG4CylinderStripSubsystem("BMT",ilayer);
                example01->set_double_param("radius", BMT_r[ilayer]);
                example01->set_string_param("gas", "myMMGas");
                //example01->set_double_param("steplimits", 300e-4);                                                                                  
                example01->set_double_param("phi0", 15*ilayer);
                example01->set_double_param("gap", gap_betweenCZ);
                example01->SetActive();
                example01->SuperDetector("BMT");
                example01->set_int_param("lengthviarapidity",0);
                example01->set_double_param("length", bmt_length);
                example01->set_double_param("deadzone", 0.2);
                example01->set_int_param("nhit", 2);
                example01->OverlapCheck(true);
                example01->set_int_param("use_2Dreadout",use_2Dreadout);
                g4Reco->registerSubsystem(example01);
                //example01->Print();                                                                                                                 
        }
        #endif


  g4Reco->SetWorldMaterial("G4_Galactic"); //For material scans as in EICROOT                                                                         

        #ifdef _BLACKHOLES_
        // ======================================================================================================                                     
        //Setting up black holes in the forward, backward, and radial directions                                                                      
        PHG4CylinderSubsystem * bhCyl;
        bhCyl = new PHG4CylinderSubsystem(radialBH,0);
        bhCyl->set_double_param("length", radialBHlength);
        bhCyl->set_double_param("radius", radialBHradius); // dirc radius                                                                             
        bhCyl->set_double_param("thickness", 0.1); // needs some thickness                                                                            
        bhCyl->set_string_param("material", "G4_AIR");
        bhCyl->SetActive(1);
        bhCyl->SuperDetector(radialBH);
        bhCyl->BlackHole();
        bhCyl->set_color(0,1,1,0.3); //blue                                                                                                           
        g4Reco->registerSubsystem(bhCyl);

        bhCyl = new PHG4CylinderSubsystem(forwardBH,0);
        bhCyl->set_double_param("length", thickness);
        bhCyl->set_double_param("radius", 0); //                                                                                                      
        bhCyl->set_double_param("thickness", forwardBHradius); //                                                                                     
        bhCyl->set_string_param("material", "G4_AIR");
 bhCyl->set_double_param("place_z", forwardBHz);
        bhCyl->SetActive(1);
        bhCyl->SuperDetector(forwardBH);
        bhCyl->BlackHole();
        bhCyl->set_color(0,1,1,0.3); //blue                                                                                                           
        g4Reco->registerSubsystem(bhCyl);

        bhCyl = new PHG4CylinderSubsystem(backwardBH,0);
        bhCyl->set_double_param("length", thickness);
        bhCyl->set_double_param("radius", 0); //                                                                                                      
        bhCyl->set_double_param("thickness", backwardBHradius); //                                                                                    
        bhCyl->set_string_param("material", "G4_AIR");
        bhCyl->set_double_param("place_z", backwardBHz);
        bhCyl->SetActive(1);
        bhCyl->SuperDetector(backwardBH);
        bhCyl->BlackHole();
        bhCyl->set_color(0,1,1,0.3); //blue                                                                                                           
        g4Reco->registerSubsystem(bhCyl);
        #endif

#ifdef _GEMS_

        // Forward GEM tracker module(s);                                                                                                            \
                                                                                                                                                      
        auto fgt = new EicRootGemSubsystem("FGT");
        {
          fgt->SetActive(true);
          //fgt->CheckOverlap();                                                                                                                     \
                                                                                                                                                      
          //fgt->SetTGeoGeometryCheckPrecision(0.000001 * etm::um);                                                                                  \
                                                                             
                    {
            auto sbs = new GemModule();
            // See other GemModule class data in GemGeoParData.h;                                                                                    \
                                                                                                                                                      
            sbs->SetDoubleVariable("mFrameBottomEdgeWidth", 30 * etm::mm);

            // Compose sectors; parameters are:                                                                                                      \
                                                                                                                                                      
            //   - layer description (obviously can mix different geometries);                                                                       \
                                                                                                                                                      
            //   - azimuthal segmentation;                                                                                                            
            //   - gas volume center radius;                                                                                                         \
                                                                                                                                                      
            //   - Z offset from 0.0 (default);                                                                                                      \
                                                                                                                                                      
            //   - azimuthal rotation from 0.0 (default);                                                                                            \
                                                                                                                                                      
            fgt->AddWheel(sbs, 12, 420 * etm::mm, 1200 * etm::mm, 0);
            fgt->AddWheel(sbs, 12, 420 * etm::mm, 1300 * etm::mm, 0);
            fgt->AddWheel(sbs, 24, 840 * etm::mm, 1500 * etm::mm, 0);
            fgt->AddWheel(sbs, 12, 420 * etm::mm, -1200 * etm::mm, 0);
            fgt->AddWheel(sbs, 12, 420 * etm::mm, -1300 * etm::mm, 0);
          }

          g4Reco->registerSubsystem(fgt);
        }
#endif
  // ======================================================================================================                                     

        //For HepMC                                                                                                                                   
        HepMCNodeReader * hr = new HepMCNodeReader();
        se->registerSubsystem(hr);

        //Register the truth hits                                                                                                                     
        PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
        g4Reco->registerSubsystem(truth);


        se->registerSubsystem(g4Reco);

  // ======================================================================================================                                     
        //Reconstruction                                                                                                                              
        // Fast pattern recognition and full Kalman filter                                                                                            
        // output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer                                                          
        PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
        kalman->set_use_vertex_in_fitting(false);
        kalman->set_sub_top_node_name("SVTX");
        kalman->set_trackmap_out_name("SvtxTrackMap");

        //10: SVT                                                                                                                                     
        //20: Time-stamping layer(s)                                                                                                                  
        //30: Silicon TPC                                                                                                                             
        //40: Forward disks                                                                                                                           
        //50: Backward disks                                                                                                                          
	double pixelResolution = 10e-4/sqrt(12);
        std::cout << pixelResolution << std::endl;
        #ifdef _BARREL_
        //SVT                                                                                                                                         
        //This adds the hit nodes for EACH LAYER OF THE DETECTOR.                                                                                     
        for (int i = 10; i < 10+_NO_OF_BARREL_LAYERS_; i++) {
                std::string nodeName = "G4HIT_SVT_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.\
cc ******                                                                                                                                             
                kalman->add_phg4hits(
                                nodeName,                               // const std::string& phg4hitsNames                                           
                                PHG4TrackFastSim::Cylinder,             // const DETECTOR_TYPE phg4dettype                                            
                                999.,                                   // radial-resolution [cm] (this number is not used in cylindrical geometry)   
                                pixelResolution,                                        // azimuthal (arc-length) resolution [cm]                     
                                pixelResolution,                                        // longitudinal (z) resolution [cm]                           
                                1,                                      // efficiency (fraction)                                                      
                                0                                       // hit noise                                                                  
                );
  }
        #endif

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

    #ifdef _MPGD_
        if (use_2Dreadout) {
                kalman->add_phg4hits(
                        "G4HIT_BMT",                //      const std::string& phg4hitsNames,                                                         
                        PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,                                                         
                        2.5/2/sqrt(12),                      //       radial-resolution [cm], only used for Vertical Plane Detector Type              
                        150e-4,                       //        azimuthal-resolution [cm]                                                             
                        150e-4,                           //      z-resolution [cm]                                                                   
                        1,                           //      efficiency,                                                                              
                        0                            //      noise hits                                                                               
                );
          }
          else {
                //Actually create the hits                                
             CreateCZHitContainer * cz = new CreateCZHitContainer("BMT");
                se->registerSubsystem(cz);
                //Then do the kalman as usual                                                                                                         
                kalman->add_phg4hits(
                        "G4HIT_CZBMT",                //      const std::string& phg4hitsNames,                                                       
                        PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,                                                         
                        2.5/2/sqrt(12),                      //       radial-resolution [cm], only used for Vertical Plane Detector Type              
                        150e-4,                       //        azimuthal-resolution [cm]                                                             
                        150e-4,                           //      z-resolution [cm]                                                                   
                        1,                           //      efficiency,                                                                              
                        0                            //      noise hits                                                                               
                );
          }
    #endif


        #ifdef _FORWARD_SILICON_DISKS_
        //Forward silicon disks                                                                                                                       
        for (int i = 40; i < 40+_NO_OF_FORWARD_DISKS_ ; i++) {
                std::string nodeName = "G4HIT_FstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDete\
ctor.cc ******                                                                                                                                        
                                                                                                                                                     \
     //Comes from SuperDetector name when adding it to Geant4.                                                                                        
                kalman->add_phg4hits(
                                nodeName,                               // const std::string& phg4hitsNames                                           
                                PHG4TrackFastSim::Vertical_Plane,               // const DETECTOR_TYPE phg4dettype                                    
                                pixelResolution,                                        // radial-resolution [cm]                                     
                                pixelResolution,                                        // azimuthal (arc-length) resolution [cm]                     
                                999.0,                                  // longitudinal (z) resolution [cm](not used for vertical plane)  
                          1,                                      // efficiency (fraction)                                                      
                                0                                       // hit noise                                                                  
                );
        }
        #endif
        #ifdef _BACKWARD_SILICON_DISKS_
        //Backward silicon disks                                                                                                                      
        for (int i = 50; i < 50+_NO_OF_BACKWARD_DISKS_ ; i++) {
                std::string nodeName = "G4HIT_BstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******                                                                                                                                        
                                                                                                                                                     //Comes from SuperDetector name when adding it to Geant4.                                                                                        
                kalman->add_phg4hits(
                                nodeName,                               // const std::string& phg4hitsNames                                           
                                PHG4TrackFastSim::Vertical_Plane,               // const DETECTOR_TYPE phg4dettype                                    
                                pixelResolution,                                        // radial-resolution [cm]                                     
                                pixelResolution,                                        // azimuthal (arc-length) resolution [cm]                     
                                999.0,                                  // longitudinal (z) resolution [cm] (not used for vertical plane)             
                                1,                                      // efficiency (fraction)                                                      
                                0                                       // hit noise                                                                  
                );
        }
        #endif
#ifdef _GEMS_
        // GEM tracker hits; should work;                                                                                                            \
                                                                                                                                                      
        kalman->add_phg4hits(fgt->GetG4HitName(),// const std::string& phg4hitsNames                                                                 \
                                                                                                                                                      
                             PHG4TrackFastSim::Vertical_Plane,// const DETECTOR_TYPE phg4dettype                                                     \
                                                                                                         
                                                                                                                                                       
                             999.,// radial-resolution [cm] (this number is not used in cylindrical geometry)                                        \
                                                                                                                                                      
                             // Say 70um resolution?; [cm];                                                                                          \
                                                                                                                                                      
                             70e-4,        // azimuthal (arc-length) resolution [cm]                                                                 \
                                                                                                                                                      
                             70e-4,        // longitudinal (z) resolution [cm]                                                                       \
                                                                                                                                                      
                             1,// efficiency (fraction)                                                                                              \
                                                                                                                                                      
                             0);// hit noise                                                                                                         \
                                                                                                                                                      

#endif


        se->registerSubsystem(kalman);

        // ======================================================================================================                                     
        TrackFastSimEval *fast_sim_eval = new TrackFastSimEval("FastTrackingEval");
        fast_sim_eval->set_filename(TString(outputFile)+B_label+"_FastTrackingEval.root");
        //fast_sim_eval->AddProjection(radialBH);                                                                                                     
        //fast_sim_eval->AddProjection(forwardBH);                                                                                                    
        //fast_sim_eval->AddProjection(backwardBH);                                                                                                   
        se->registerSubsystem(fast_sim_eval);


  SimpleNtuple * hits = new SimpleNtuple("Hits");
        #ifdef _BARREL_
        //SVT starts at 10, to allow for absorbers before                                                                                             
        //hits->AddNode("ABSORBER_SVT",0); // hits in the passive volumes                                                                             
        for (int i = 10; i < 10+_NO_OF_BARREL_LAYERS_; i++) { // hits in the  MimosaCore volumes                                                      
                std::string nodeName = "SVT_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of t\
he G4HIT, because currently SimpleNtuple adds that.                                                                                                   
                                                                                                                   // See SimpleNtuple.cc ******      
                hits->AddNode(nodeName, i);
        }
        #endif
        #ifdef _TPC_
                hits->AddNode("TPC", 60);
        #endif
        #ifdef _MPGD_
                if (use_2Dreadout) {
                        hits->AddNode("BMT", 60);
                } else {
                        hits->AddNode("CZBMT", 60);
                }
        #endif
 #ifdef _FORWARD_SILICON_DISKS_
        //Forward silicon disks                                                                                                                       
        //hits->AddNode("ABSORBER_FstDisks",3); // hits in the passive volumes                                                                        
        for (int i = 40; i < 40+_NO_OF_FORWARD_DISKS_; i++) { // hits in the  MimosaCore volumes                                                      
                std::string nodeName = "FstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped\
 of the G4HIT, because currently SimpleNtuple adds that.                                                                                              
                                                                                                                   // See SimpleNtuple.cc ******      
                hits->AddNode(nodeName, i);
        }
        #endif
 #ifdef _BACKWARD_SILICON_DISKS_
        //Backward silicon disks                                                                                                                      
        //hits->AddNode("ABSORBER_BstDisks",4); // hits in the passive volumes                                                                        
        for (int i = 50; i < 50+_NO_OF_FORWARD_DISKS_; i++) { // hits in the  MimosaCore volumes                                                      
                std::string nodeName = "BstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped\
 of the G4HIT, because currently SimpleNtuple adds that.                                                                                              
                                                                                                                   // See SimpleNtuple.cc ******      
                hits->AddNode(nodeName, i);
        }
        #endif
        /*#ifdef _GEMS_                                                                                                                               
        hits->AddNode("GEMS",2);                                                                                                                      
        #endif  */
        se->registerSubsystem(hits);

  ///////////////////////////////////////////                                                                                                   
        // IOManagers...                                                                                                                              
        ///////////////////////////////////////////                                                                                                   
        const std::string dst_name = std::string(outputFile)+"_gdmlimporter.root";
        //Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",TString(outputFile)+"_gdmlimporter.root");                              
        Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",dst_name);
        out->Verbosity(10);
        //Let's not output full DST tree right now                                                                                                    
        //se->registerOutputManager(out);                                                                                                             

        if (particle_gen != 1 && particle_gen != 2) {

                Fun4AllInputManager *in2 = new Fun4AllHepMCInputManager("DSTIN");
                se->registerInputManager(in2);
                se->fileopen(in2->Name().c_str(), "../../eventFiles/ggm2ccbar_28_1_100000.dat"); //HepMC event file                                   
        } else {
                Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
                se->registerInputManager(in);
        }

  if (nEvents <= 0)
        {
                return;
        }
        se->run(nEvents);
        se->End();
        delete se;
        gSystem->Exit(0);
}






