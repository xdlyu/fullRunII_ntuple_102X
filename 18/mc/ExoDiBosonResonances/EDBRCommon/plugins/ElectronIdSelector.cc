
/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   Chayanit Asawatangtrakuldee chayanit@cern.ch 
 *
 * Description:
 *   - Selects "loose" and "tight" electrons needed for V-boson analysis.
 *   - Saves collection of the reference vectors of electrons passing the 
 *     required electron ID.
 * History:
 *   
 *
 *****************************************************************************/
////////////////////////////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
class ElectronIdSelector : public edm::EDProducer
{
public:
  // construction/destruction
  ElectronIdSelector(const edm::ParameterSet& iConfig);
  virtual ~ElectronIdSelector();
  
  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();

private:  
  // member data
  //edm::InputTag  src_;
  std::string    moduleLabel_;
  std::string    idLabel_;  
  bool           useDetectorIsolation_;
  bool           applyTightID_;
  bool           applyMediumID_;
  bool           applyLooseID_;
  bool           applyVetoID_;
  unsigned int nTot_;
  unsigned int nPassed_;
  edm::EDGetTokenT<edm::View<pat::Electron>> ElectronToken_;
  edm::EDGetTokenT<reco::VertexCollection> VertexToken_;
  edm::EDGetTokenT<double> RhoToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > vidToken_; //VID is versioned ID, is the standard E/gamma ID producer which we have configured for HEEP
  EffectiveAreas effectiveAreas_;
};



////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
ElectronIdSelector::ElectronIdSelector(const edm::ParameterSet& iConfig)
//  : src_    (iConfig.getParameter<edm::InputTag>     ("src"))
  : moduleLabel_(iConfig.getParameter<std::string>   ("@module_label"))
  , idLabel_(iConfig.existsAs<std::string>("idLabel") ? iConfig.getParameter<std::string>("idLabel") : "loose")
  , useDetectorIsolation_(iConfig.existsAs<bool>("useDetectorIsolation") ? iConfig.getParameter<bool>("useDetectorIsolation") : false)
  , nTot_(0)
  , nPassed_(0)
  , ElectronToken_ (consumes<edm::View<pat::Electron>> (iConfig.getParameter<edm::InputTag>( "src" ) ) )
  , VertexToken_ (consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>( "vertex" ) ) )
  , RhoToken_ (consumes<double> (iConfig.getParameter<edm::InputTag>( "rho") ) )
  , vidToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("vid")))
  , effectiveAreas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )
{
  produces<std::vector<pat::Electron> >();

  /// ------- Decode the ID criteria --------
  applyTightID_ = false;
  applyMediumID_ = false;
  applyLooseID_ = false;
  applyVetoID_ = false;

  if( (idLabel_.compare("tight")==0) || 
      (idLabel_.compare("Tight")==0) || 
      (idLabel_.compare("TIGHT")==0) ||
      (idLabel_.compare("WP70")==0) ||
      (idLabel_.compare("wp70")==0) )  
    applyTightID_ = true;
  else if( (idLabel_.compare("medium")==0) ||
      (idLabel_.compare("Medium")==0) ||
      (idLabel_.compare("MEDIUM")==0) ||
      (idLabel_.compare("WP80")==0) ||
      (idLabel_.compare("wp80")==0) )  applyMediumID_ = true;
  else if( (idLabel_.compare("loose")==0) || 
      (idLabel_.compare("Loose")==0) || 
      (idLabel_.compare("LOOSE")==0) ||
      (idLabel_.compare("WP90")==0) ||
      (idLabel_.compare("wp90")==0) )  applyLooseID_ = true;
  else if( (idLabel_.compare("veto")==0) || 
      (idLabel_.compare("Veto")==0) || 
      (idLabel_.compare("VETO")==0) ||
      (idLabel_.compare("VETOid")==0) ||
      (idLabel_.compare("VetoId")==0) )  applyVetoID_ = true;
}

 
//______________________________________________________________________________
ElectronIdSelector::~ElectronIdSelector(){}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////
 
//______________________________________________________________________________
void ElectronIdSelector::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{

   edm::Handle<reco::VertexCollection> vtxs;
   iEvent.getByToken(VertexToken_, vtxs);

   std::unique_ptr<std::vector<pat::Electron> > passingElectrons(new std::vector<pat::Electron >);
  
   edm::Handle<edm::View<pat::Electron> > electrons;
   iEvent.getByToken(ElectronToken_, electrons);  

   edm::Handle<edm::ValueMap<bool> > vid;
   iEvent.getByToken(vidToken_,vid);

   bool* isPassing = new bool[electrons->size()];


   for(unsigned int iElec=0; iElec<electrons->size(); iElec++) { 

    isPassing[iElec]=false;

    const pat::Electron& ele = electrons->at(iElec);
    const auto elo = electrons->ptrAt(iElec);

    bool isTight  = false;  
    bool isMedium = false;  
    bool isLoose  = false;  
    bool isVeto    = false; 

    float pt  = ele.pt();

    isTight = (pt > 55) && ((*vid)[elo]);
    isLoose = (pt > 35) && ((*vid)[elo]);

    /// ------- Finally apply selection --------
    if(applyTightID_ && isTight)   isPassing[iElec]= true;
    if(applyMediumID_ && isMedium) isPassing[iElec]= true;
    if(applyLooseID_ && isLoose)   isPassing[iElec]= true;
    if(applyVetoID_ && isVeto) isPassing[iElec]= true;
    
 }
  

 for (unsigned int iElectron = 0; iElectron < electrons -> size(); iElectron ++)
   {     if(isPassing[iElectron]) passingElectrons->push_back( electrons -> at(iElectron) );       
  }
 
  nTot_  +=electrons->size();
  nPassed_+=passingElectrons->size();

  delete [] isPassing;  
  iEvent.put(std::move(passingElectrons));
}

 
//______________________________________________________________________________
void ElectronIdSelector::endJob()
{
  std::stringstream ss;
  ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_
    <<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   <<"\n"<<moduleLabel_<<"(ElectronIdSelector) SUMMARY:\n"<<ss.str()
	   <<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////
typedef ElectronIdSelector   			    PATElectronIdSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATElectronIdSelector);
