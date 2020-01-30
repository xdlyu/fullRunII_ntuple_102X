// system include files

#include <iostream>
#include <memory>
#include "TMath.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include<algorithm>
#define Pi 3.141593
#include "Math/VectorUtil.h"
#include "TMath.h"
#include <TFormula.h>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
struct ZsortPt
{
	bool operator()(TLorentzVector* s1, TLorentzVector* s2) const
	{
		return s1->Pt() >= s2->Pt();
	}
} ZmysortPt;
//
// class declaration
//

class ZPKUTreeMaker : public edm::EDAnalyzer {
	public:
		explicit ZPKUTreeMaker(const edm::ParameterSet&);
		~ZPKUTreeMaker();
		//static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

		enum PhotonMatchType {UNMATCHED = 0,
			MATCHED_FROM_GUDSCB,
			MATCHED_FROM_PI0,
			MATCHED_FROM_OTHER_SOURCES};

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
		virtual void endRun(const edm::Run&, const edm::EventSetup&) override;
		virtual void addTypeICorr( edm::Event const & event );
//		virtual void addTypeICorr_user( edm::Event const & event );//---for MET, Meng
		virtual double getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
		virtual double getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
		math::XYZTLorentzVector getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType);
		bool hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<edm::View<pat::Electron> > &eleCol,const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,float lxyMin=2.0, float probMin=1e-6, unsigned int nHitsBeforeVtxMax=0);
		int matchToTruth(const reco::Photon &pho,const edm::Handle<edm::View<reco::GenParticle>>  &genParticles, bool &ISRPho, double &dR, int &isprompt);

		void findFirstNonPhotonMother(const reco::Candidate *particle,int &ancestorPID, int &ancestorStatus);
		// muon station2 retrieve, L1 issue, Meng 2017/3/26
		std::pair<double,double> EtaPhiAtME2X(const pat::Muon *iM, const PropagateToMuon *propagatetomuon);
		std::pair<double,double> lep1_etaphi_;
		std::pair<double,double> lep2_etaphi_;
		double lep1_eta_station2;
		double lep1_phi_station2;
		int lep1_sign;
		double lep2_eta_station2;
		double lep2_phi_station2;
		int lep2_sign;
		edm::EDGetTokenT<edm::View<pat::Muon> > goodmuonToken_;
		PropagateToMuon *muPropagator2nd_;
		// Lu

		float EAch(float x); 
		float EAnh(float x);
		float EApho(float x);
		std::vector<std::string> offsetCorrLabel_;
		FactorizedJetCorrector* jecOffset_;
		std::vector<std::string> jetCorrLabel_;
		edm::Handle< double >  rho_;
		edm::EDGetTokenT<double> rhoToken_;
		edm::EDGetTokenT<pat::METCollection>  metInputToken_;
		std::vector<edm::EDGetTokenT<pat::METCollection>> mettokens;
		edm::EDGetTokenT<pat::METCollection> metToken_;
		edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_ ;
		edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;
		edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
		edm::EDGetTokenT<std::vector<reco::Conversion> > conversionsToken_;
		edm::EDGetTokenT<edm::View<pat::Electron> > looseelectronToken_ ; 
		edm::EDGetTokenT<edm::View<pat::Muon> > loosemuonToken_; 

		// Filter
		edm::EDGetTokenT<edm::TriggerResults> 		     noiseFilterToken_;
		edm::Handle< edm::TriggerResults> 			     noiseFilterBits_;
		std::string HBHENoiseFilter_Selector_;
		std::string HBHENoiseIsoFilter_Selector_;
		std::string ECALDeadCellNoiseFilter_Selector_;
		std::string GoodVtxNoiseFilter_Selector_;
		std::string EEBadScNoiseFilter_Selector_;
		std::string globalTightHaloFilter_Selector_;
		edm::EDGetTokenT<bool>  badMuon_Selector_;
		edm::EDGetTokenT<bool>  badChargedHadron_Selector_;

		// Filter, Meng, bad muon, duplicate muon
		std::string badMuonFilter_Selector_;
		std::string duplicateMuonFilter_Selector_;


		edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMapToken_;
		edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_;
		edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_;
		edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_;
		EffectiveAreas effAreaChHadrons_;
		EffectiveAreas effAreaNeuHadrons_;
		EffectiveAreas effAreaPhotons_;

		//L1 prefiring
		edm::EDGetTokenT< double > prefweight_token;
		edm::EDGetTokenT< double > prefweightup_token;
		edm::EDGetTokenT< double > prefweightdown_token;

		// ----------member data ---------------------------
		TH1D*hout;
		TTree* outTree_;

		//L1 prefiring
		double prefiringweight, prefiringweightup, prefiringweightdown;

		int nevent, run, ls;
		int nVtx;
		double triggerWeight, lumiWeight, pileupWeight;
		double theWeight;
		double  nump=0.;
		double  numm=0.;
//		double pweight[703];
		double  npT, npIT;
		int     nBX;
		double ptVlep, yVlep, phiVlep, massVlep;
		double Mla, Mla2, Mva;
		double Mla_f, Mla2_f, Mva_f;
		double ptlep1, etalep1, philep1;
		double ptlep2, etalep2, philep2;
		// for muon rochester correction
		int muon1_trackerLayers;
		double matchedgenMu1_pt;
		int muon2_trackerLayers;
		double matchedgenMu2_pt;
		double getDR(double eta1, double phi1, double eta2, double phi2);
		// for muon rochester correction
		int  lep, nlooseeles,nloosemus, ngoodmus;
		double met, metPhi, j1metPhi, j2metPhi;
		double j1metPhi_f, j2metPhi_f;
		//Met JEC
		double METraw_et, METraw_phi, METraw_sumEt;
		double genMET, MET_et, MET_phi, MET_sumEt, MET_corrPx, MET_corrPy;
		double useless;
		// AK4 Jets
		double ak4jet_pt[6],ak4jet_eta[6],ak4jet_phi[6],ak4jet_e[6];
		double ak4jet_pt_jer[6];
		double ak4jet_csv[6],ak4jet_icsv[6];
		double drjetlep[6], drjetphoton[6];
		double genphoton_pt[6],genphoton_eta[6],genphoton_phi[6];
		double genjet_pt[6],genjet_eta[6],genjet_phi[6],genjet_e[6];
		double genmuon_pt[6],genmuon_eta[6],genmuon_phi[6];
		int genmuon_pid[6];
		double genelectron_pt[6],genelectron_eta[6],genelectron_phi[6];
		//Photon
		double photon_pt[6],photon_eta[6],photon_phi[6],photon_e[6], photon_sc_eta[6];
		bool   photon_pev[6],photon_pevnew[6],photon_ppsv[6],photon_iseb[6],photon_isee[6];
		double photon_hoe[6],photon_sieie[6],photon_sieie2[6],photon_chiso[6],photon_nhiso[6],photon_phoiso[6],photon_drla[6],photon_drla2[6],photon_mla[6],photon_mla2[6],photon_mva[6];
		int      photon_istrue[6], photon_isprompt[6];
		double photonet, photoneta, photonphi, photone, photoneta_sc;
		double photonet_f, photoneta_f, photonphi_f, photone_f, photoneta_sc_f;
		double photonsieie, photonphoiso, photonchiso, photonnhiso;
		double photonsieie_f, photonphoiso_f, photonchiso_f, photonnhiso_f;
		int iphoton;
		int iphoton_f;
		double drla,drla2;
		double drla_f,drla2_f;
		bool passEleVeto, passEleVetonew, passPixelSeedVeto;
		//Photon gen match
		int   isTrue_;
		bool ISRPho;
		int isprompt_;
		double dR_;
		//Jets
		double jet1pt, jet1eta, jet1phi, jet1e, jet1csv, jet1icsv;
		double jet1pt_f, jet1eta_f, jet1phi_f, jet1e_f, jet1csv_f, jet1icsv_f;
		double jet2pt, jet2eta, jet2phi, jet2e, jet2csv, jet2icsv;
		double jet2pt_f, jet2eta_f, jet2phi_f, jet2e_f, jet2csv_f, jet2icsv_f;
		double drj1a, drj2a, drj1l, drj2l, drj1l2,drj2l2;
		double drj1a_f, drj2a_f, drj1l_f, drj2l_f, drj1l2_f, drj2l2_f;
		double Mjj, deltaetajj, zepp;
		double Mjj_f, deltaetajj_f, zepp_f; 

		void setDummyValues();

		/// Parameters to steer the treeDumper
		int originalNEvents_;
		double crossSectionPb_;
		double targetLumiInvPb_;
		std::string PKUChannel_;
		bool isGen_ , RunOnMC_;
		std::vector<std::string> jecAK4Labels_;
		std::vector<std::string> jecAK4chsLabels_;
		//correction jet
		FactorizedJetCorrector* jecAK4_;
		std::string gravitonSrc_;
		std::map<std::string,double>  TypeICorrMap_;
		edm::InputTag mets_;
		//High Level Trigger
		HLTConfigProvider hltConfig;
		edm::EDGetTokenT<edm::TriggerResults> hltToken_;
		std::vector<std::string> elPaths1_, elPaths2_, elPaths3_, elPaths4_, elPaths5_, elPaths6_, elPaths7_;
		std::vector<std::string> muPaths1_, muPaths2_, muPaths3_, muPaths4_, muPaths5_, muPaths6_, muPaths7_, muPaths8_, muPaths9_, muPaths10_, muPaths11_;;
		std::vector<std::string> elPaths1, elPaths2, elPaths3, elPaths4, elPaths5, elPaths6, elPaths7;
		std::vector<std::string> muPaths1, muPaths2, muPaths3, muPaths4,muPaths5, muPaths6, muPaths7, muPaths8,muPaths9, muPaths10, muPaths11;
		int  HLT_Ele1;
		int  HLT_Ele2;
		int  HLT_Ele3;
		int  HLT_Ele4;
		int  HLT_Ele5;
		int  HLT_Ele6;
		int  HLT_Ele7;
		int  HLT_Mu1;
		int  HLT_Mu2;
		int  HLT_Mu3;
		int  HLT_Mu4;
		int  HLT_Mu5;
		int  HLT_Mu6;
		int  HLT_Mu7;
		int  HLT_Mu8;
		int  HLT_Mu9;
		int  HLT_Mu10;
		int  HLT_Mu11;
		// filter
		bool passFilter_HBHE_                   ;
		bool passFilter_HBHEIso_                ;
		bool passFilter_globalTightHalo_ ;
		bool passFilter_ECALDeadCell_           ;
		bool passFilter_GoodVtx_                ;
		bool passFilter_EEBadSc_                ;
		bool passFilter_badMuon_                ;
		bool passFilter_badChargedHadron_       ;
		// Meng
		bool passFilter_MetbadMuon_		  ;
		bool passFilter_duplicateMuon_	  ; 

		edm::EDGetTokenT<GenEventInfoProduct> GenToken_;
		edm::EDGetTokenT<LHEEventProduct> LheToken_;
		edm::EDGetTokenT<reco::GenJetCollection> genJet_;
		edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUToken_;
		edm::EDGetTokenT<edm::View<reco::Candidate>> leptonicVSrc_;
		edm::EDGetTokenT<edm::View<pat::Jet>> ak4jetsSrc_;
		edm::EDGetTokenT<edm::View<pat::Photon>> photonSrc_;
		edm::EDGetTokenT<edm::View<reco::GenParticle>> genSrc_;
		edm::EDGetTokenT<edm::View<reco::Candidate>> metSrc_;
		edm::EDGetTokenT<reco::VertexCollection> VertexToken_;
		edm::EDGetTokenT<pat::JetCollection> t1jetSrc_;
		edm::EDGetTokenT<edm::View<pat::Muon>> t1muSrc_;

};


float ZPKUTreeMaker::EAch( float x){
	float EA = 0.0360;
	if(x>1.0)   EA = 0.0377;
	if(x>1.479) EA = 0.0306;
	if(x>2.0)   EA = 0.0283;
	if(x>2.2)   EA = 0.0254;
	if(x>2.3)   EA = 0.0217;
	if(x>2.4)   EA = 0.0167;
	return EA;
}
float ZPKUTreeMaker::EAnh( float x){
	float EA = 0.0597;
	if(x>1.0)   EA = 0.0807;
	if(x>1.479) EA = 0.0629;
	if(x>2.0)   EA = 0.0197;
	if(x>2.2)   EA = 0.0184;
	if(x>2.3)   EA = 0.0284;
	if(x>2.4)   EA = 0.0591;
	return EA;
}
float ZPKUTreeMaker::EApho( float x){
	float EA = 0.1210;
	if(x>1.0)   EA = 0.1107;
	if(x>1.479) EA = 0.0699;
	if(x>2.0)   EA = 0.1056;
	if(x>2.2)   EA = 0.1457;
	if(x>2.3)   EA = 0.1719;
	if(x>2.4)   EA = 0.1998;
	return EA;
}

//muon station2 retrieve, L1 issue, Meng 2017/3/26
std::pair<double,double> ZPKUTreeMaker::EtaPhiAtME2X(const pat::Muon *iM, const PropagateToMuon *propagatetomuon){
	std::pair<double,double>  etaphi(0.,0.);
	TrajectoryStateOnSurface stateAtMuSt2 = propagatetomuon->extrapolate(*iM);
	if (stateAtMuSt2.isValid()) {
		etaphi.first = stateAtMuSt2.globalPosition().eta();
		etaphi.second = stateAtMuSt2.globalPosition().phi();
	}
	return etaphi;
}
// Lu

//
// constructors and destructor
//
ZPKUTreeMaker::ZPKUTreeMaker(const edm::ParameterSet& iConfig)//:
	:effAreaChHadrons_((iConfig.getParameter<edm::FileInPath>("effAreaChHadFile")).fullPath() )
	 ,effAreaNeuHadrons_((iConfig.getParameter<edm::FileInPath>("effAreaNeuHadFile")).fullPath() )
	 ,effAreaPhotons_((iConfig.getParameter<edm::FileInPath>("effAreaPhoFile")).fullPath() )
{
	hltToken_=consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltToken"));
	elPaths1_=iConfig.getParameter<std::vector<std::string>>("elPaths1");
	elPaths2_=iConfig.getParameter<std::vector<std::string>>("elPaths2");
	elPaths3_=iConfig.getParameter<std::vector<std::string>>("elPaths3");
	elPaths4_=iConfig.getParameter<std::vector<std::string>>("elPaths4");
	elPaths5_=iConfig.getParameter<std::vector<std::string>>("elPaths5");
	elPaths6_=iConfig.getParameter<std::vector<std::string>>("elPaths6");
	elPaths7_=iConfig.getParameter<std::vector<std::string>>("elPaths7");
	muPaths1_=iConfig.getParameter<std::vector<std::string>>("muPaths1");
	muPaths2_=iConfig.getParameter<std::vector<std::string>>("muPaths2");
	muPaths3_=iConfig.getParameter<std::vector<std::string>>("muPaths3");
	muPaths4_=iConfig.getParameter<std::vector<std::string>>("muPaths4");
	muPaths5_=iConfig.getParameter<std::vector<std::string>>("muPaths5");
	muPaths6_=iConfig.getParameter<std::vector<std::string>>("muPaths6");
	muPaths7_=iConfig.getParameter<std::vector<std::string>>("muPaths7");
	muPaths8_=iConfig.getParameter<std::vector<std::string>>("muPaths8");
	muPaths9_=iConfig.getParameter<std::vector<std::string>>("muPaths9");
	muPaths10_=iConfig.getParameter<std::vector<std::string>>("muPaths10");
	muPaths11_=iConfig.getParameter<std::vector<std::string>>("muPaths11");
	GenToken_=consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "generator") ) ;
	LheToken_=consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>( "lhe") ) ;
	genJet_=consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJet"));
	PUToken_=consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup") ) ;
	leptonicVSrc_=consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>( "leptonicVSrc") ) ;
	ak4jetsSrc_      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>( "ak4jetsSrc") ) ;
	photonSrc_      = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>( "photonSrc") ) ;
	genSrc_      = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genSrc") ) ;
	metSrc_      = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>( "metSrc") ) ;
	VertexToken_ =consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>( "vertex" ) ) ;
	t1jetSrc_      = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>( "t1jetSrc") ) ;
	t1muSrc_      = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>( "t1muSrc") ) ;
	originalNEvents_ = iConfig.getParameter<int>("originalNEvents");
	crossSectionPb_  = iConfig.getParameter<double>("crossSectionPb");
	targetLumiInvPb_ = iConfig.getParameter<double>("targetLumiInvPb");
	PKUChannel_     = iConfig.getParameter<std::string>("PKUChannel");
	isGen_           = iConfig.getParameter<bool>("isGen");
	RunOnMC_           = iConfig.getParameter<bool>("RunOnMC");
	rhoToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
	jecAK4chsLabels_   =  iConfig.getParameter<std::vector<std::string>>("jecAK4chsPayloadNames");
	jecAK4Labels_   =  iConfig.getParameter<std::vector<std::string>>("jecAK4PayloadNames");
	metToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metSrc"));
	mettokens.push_back( metToken_ );
	metInputToken_ = mettokens[0];
	electronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("electrons")))            ;
	looseelectronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("looseelectronSrc"))) ;
	loosemuonToken_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("loosemuonSrc")))              ;
	goodmuonToken_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("goodmuonSrc")))              ;
	beamSpotToken_    = (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))) ;
	conversionsToken_ = (consumes<std::vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))) ;


	jetCorrLabel_ = jecAK4chsLabels_;
	offsetCorrLabel_.push_back(jetCorrLabel_[0]);

	// filter
	noiseFilterToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"));
	HBHENoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseFilter");
	HBHENoiseIsoFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseIsoFilter");  
	globalTightHaloFilter_Selector_ = iConfig.getParameter<std::string> ("noiseFilterSelection_globalTightHaloFilter"); 
	ECALDeadCellNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter");
	GoodVtxNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_goodVertices");
	EEBadScNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_eeBadScFilter");
	badMuon_Selector_ =  consumes<bool>(iConfig.getParameter<edm::InputTag> ("noiseFilterSelection_badMuon"));
	badChargedHadron_Selector_ =  consumes<bool>(iConfig.getParameter<edm::InputTag> ("noiseFilterSelection_badChargedHadron"));
	full5x5SigmaIEtaIEtaMapToken_=(consumes <edm::ValueMap<float> >
			(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap")));
	phoChargedIsolationToken_=(consumes <edm::ValueMap<float> >
			(iConfig.getParameter<edm::InputTag>("phoChargedIsolation")));
	phoNeutralHadronIsolationToken_=(consumes <edm::ValueMap<float> >
			(iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation")));
	phoPhotonIsolationToken_=(consumes <edm::ValueMap<float> >
			(iConfig.getParameter<edm::InputTag>("phoPhotonIsolation")));
	// filter Meng
	badMuonFilter_Selector_ = iConfig.getParameter<std::string> ("badMuonFilterSelection"); 
	duplicateMuonFilter_Selector_ = iConfig.getParameter<std::string> ("duplicateMuonFilterSelection");

	// L1 prefiring
	prefweight_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProb"));
	prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProbUp"));
	prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProbDown"));

	//now do what ever initialization is needed
	edm::Service<TFileService> fs;
	hout = fs->make<TH1D>("hist","hist",2,-1,1);
	outTree_ = fs->make<TTree>("ZPKUCandidates","ZPKU Candidates");

	/// Basic event quantities
	outTree_->Branch("event"           ,&nevent         ,"event/I"          );
	outTree_->Branch("run"           ,&run         ,"run/I"          );
	outTree_->Branch("ls"           ,&ls         ,"ls/I"          );
	outTree_->Branch("nVtx"            ,&nVtx           ,"nVtx/I"           );
	outTree_->Branch("theWeight"           ,&theWeight         ,"theWeight/D"          );
	outTree_->Branch("nump"           ,&nump         ,"nump/D"          );
	outTree_->Branch("numm"           ,&numm         ,"numm/D"          );
//	outTree_->Branch("pweight"           ,pweight         ,"pweight[703]/D"          );
	outTree_->Branch("npT"           ,&npT         ,"npT/D"          );
	outTree_->Branch("lep"             ,&lep            ,"lep/I"            );
	outTree_->Branch("ptVlep"          ,&ptVlep         ,"ptVlep/D"         );
	outTree_->Branch("yVlep"           ,&yVlep          ,"yVlep/D"          );
	outTree_->Branch("phiVlep"         ,&phiVlep        ,"phiVlep/D"        );
	outTree_->Branch("massVlep"        ,&massVlep       ,"massVlep/D"       );
	outTree_->Branch("Mla"          ,&Mla         ,"Mla/D"         );
	outTree_->Branch("Mla_f"          ,&Mla_f         ,"Mla_f/D"         );
	outTree_->Branch("Mla2"          ,&Mla2         ,"Mla2/D"         );
	outTree_->Branch("Mla2_f"          ,&Mla2_f         ,"Mla2_f/D"         );
	outTree_->Branch("Mva"          ,&Mva         ,"Mva/D"         );
	outTree_->Branch("Mva_f"          ,&Mva_f         ,"Mva_f/D"         );
	outTree_->Branch("nlooseeles"          ,&nlooseeles         ,"nlooseeles/I"         );
	outTree_->Branch("nloosemus"          ,&nloosemus         ,"nloosemus/I"         );
	outTree_->Branch("ngoodmus"          ,&ngoodmus         ,"ngoodmus/I"         );
	outTree_->Branch("genphoton_pt"        , genphoton_pt       ,"genphoton_pt[6]/D"       );
	outTree_->Branch("genphoton_eta"        , genphoton_eta       ,"genphoton_eta[6]/D"       );
	outTree_->Branch("genphoton_phi"        , genphoton_phi       ,"genphoton_phi[6]/D"       );
	outTree_->Branch("genjet_pt"            ,genjet_pt             ,"genjet_pt[6]/D"           );
	outTree_->Branch("genjet_eta"            ,genjet_eta             ,"genjet_eta[6]/D"           );
	outTree_->Branch("genjet_phi"            ,genjet_phi             ,"genjet_phi[6]/D"           );
	outTree_->Branch("genjet_e"           , genjet_e             ,"genjet_e[6]/D"           );
	outTree_->Branch("genmuon_pt"        , genmuon_pt       ,"genmuon_pt[6]/D"       );
	outTree_->Branch("genmuon_eta"        , genmuon_eta       ,"genmuon_eta[6]/D"       );
	outTree_->Branch("genmuon_phi"        , genmuon_phi       ,"genmuon_phi[6]/D"       );
	outTree_->Branch("genmuon_pid"        , genmuon_pid       ,"genmuon_pid[6]/I"       );
	outTree_->Branch("genelectron_pt"        , genelectron_pt       ,"genelectron_pt[6]/D"       );
	outTree_->Branch("genelectron_eta"        , genelectron_eta       ,"genelectron_eta[6]/D"       );
	outTree_->Branch("genelectron_phi"        , genelectron_phi       ,"genelectron_phi[6]/D"       );
	/// Photon
	outTree_->Branch("photon_pt"        , photon_pt       ,"photon_pt[6]/D"       );
	outTree_->Branch("photon_eta"        , photon_eta       ,"photon_eta[6]/D"       );
	outTree_->Branch("photon_sc_eta"        , photon_sc_eta       ,"photon_sc_eta[6]/D"       );
	outTree_->Branch("photon_phi"        , photon_phi       ,"photon_phi[6]/D"       );
	outTree_->Branch("photon_e"        , photon_e       ,"photon_e[6]/D"       );
	outTree_->Branch("photon_pev"        , photon_pev       ,"photon_pev[6]/O"       );
	outTree_->Branch("photon_pevnew"        , photon_pevnew       ,"photon_pevnew[6]/O"       );
	outTree_->Branch("photon_ppsv"        , photon_ppsv       ,"photon_ppsv[6]/O"       );
	outTree_->Branch("photon_iseb"        , photon_iseb       ,"photon_iseb[6]/O"       );
	outTree_->Branch("photon_isee"        , photon_isee       ,"photon_isee[6]/O"       );
	outTree_->Branch("photon_hoe"        , photon_hoe       ,"photon_hoe[6]/D"       );
	outTree_->Branch("photon_sieie"        , photon_sieie       ,"photon_sieie[6]/D"       );
	outTree_->Branch("photon_sieie2"        , photon_sieie2       ,"photon_sieie2[6]/D"       );
	outTree_->Branch("photon_chiso"        , photon_chiso       ,"photon_chiso[6]/D"       );
	outTree_->Branch("photon_nhiso"        , photon_nhiso       ,"photon_nhiso[6]/D"       );
	outTree_->Branch("photon_phoiso"        , photon_phoiso       ,"photon_phoiso[6]/D"       );
	outTree_->Branch("photon_istrue"        , photon_istrue       ,"photon_istrue[6]/I"       );
	outTree_->Branch("photon_isprompt"        , photon_isprompt       ,"photon_isprompt[6]/I"       );
	outTree_->Branch("photon_drla"        , photon_drla       ,"photon_drla[6]/D"       );
	outTree_->Branch("photon_drla2"        , photon_drla2       ,"photon_drla2[6]/D"       );
	outTree_->Branch("photon_mla"        , photon_mla       ,"photon_mla[6]/D"       );
	outTree_->Branch("photon_mla2"        , photon_mla2       ,"photon_mla2[6]/D"       );
	outTree_->Branch("photon_mva"        , photon_mva       ,"photon_mva[6]/D"       );
	outTree_->Branch("passEleVeto"        , &passEleVeto       ,"passEleVeto/O"       );
	outTree_->Branch("passEleVetonew"        , &passEleVetonew       ,"passEleVetonew/O"       );
	outTree_->Branch("passPixelSeedVeto"        , &passPixelSeedVeto       ,"passPixelSeedVeto/O"       );
	outTree_->Branch("photonet"          ,&photonet         ,"photonet/D"         );
	outTree_->Branch("photonet_f"          ,&photonet_f         ,"photonet_f/D"         );
	outTree_->Branch("photoneta"          ,&photoneta         ,"photoneta/D"         );
	outTree_->Branch("photoneta_f"          ,&photoneta_f         ,"photoneta_f/D"         );
	outTree_->Branch("photoneta_sc"          ,&photoneta_sc         ,"photoneta_sc/D"         );
	outTree_->Branch("photoneta_sc_f"          ,&photoneta_sc_f         ,"photoneta_sc_f/D"         );
	outTree_->Branch("photonphi"          ,&photonphi         ,"photonphi/D"         );
	outTree_->Branch("photonphi_f"          ,&photonphi_f         ,"photonphi_f/D"         );
	outTree_->Branch("photone"          ,&photone         ,"photone/D"         );
	outTree_->Branch("photone_f"          ,&photone_f         ,"photone_f/D"         );
	outTree_->Branch("photonsieie"          ,&photonsieie         ,"photonsieie/D"         );
	outTree_->Branch("photonsieie_f"          ,&photonsieie_f         ,"photonsieie_f/D"         );
	outTree_->Branch("photonphoiso"          ,&photonphoiso         ,"photonphoiso/D"         );
	outTree_->Branch("photonphoiso_f"          ,&photonphoiso_f         ,"photonphoiso_f/D"         );
	outTree_->Branch("photonchiso"          ,&photonchiso         ,"photonchiso/D"         );
	outTree_->Branch("photonchiso_f"          ,&photonchiso_f         ,"photonchiso_f/D"         );
	outTree_->Branch("photonnhiso"          ,&photonnhiso         ,"photonnhiso/D"         );
	outTree_->Branch("photonnhiso_f"          ,&photonnhiso_f         ,"photonnhiso_f/D"         );
	outTree_->Branch("iphoton"             ,&iphoton            ,"iphoton/I"            );
	outTree_->Branch("iphoton_f"             ,&iphoton_f            ,"iphoton_f/I"            );
	outTree_->Branch("drla"          ,&drla         ,"drla/D"         );
	outTree_->Branch("drla_f"          ,&drla_f         ,"drla_f/D"         );
	outTree_->Branch("drla2"          ,&drla2         ,"drla2/D"         );
	outTree_->Branch("drla2_f"          ,&drla2_f         ,"drla2_f/D"         );
	outTree_->Branch("isTrue", &isTrue_, "isTrue/I");
	outTree_->Branch("isprompt"    , &isprompt_, "isprompt/I");
	//jets
	outTree_->Branch("ak4jet_pt"        , ak4jet_pt       ,"ak4jet_pt[6]/D"       );
	outTree_->Branch("ak4jet_eta"        , ak4jet_eta       ,"ak4jet_eta[6]/D"       );
	outTree_->Branch("ak4jet_phi"        , ak4jet_phi       ,"ak4jet_phi[6]/D"       );
	outTree_->Branch("ak4jet_e"        , ak4jet_e       ,"ak4jet_e[6]/D"       );
	outTree_->Branch("ak4jet_csv"        , ak4jet_csv       ,"ak4jet_csv[6]/D"       );
	outTree_->Branch("ak4jet_icsv"        , ak4jet_icsv       ,"ak4jet_icsv[6]/D"       );
	outTree_->Branch("jet1pt"          ,&jet1pt         ,"jet1pt/D"         );
	outTree_->Branch("jet1pt_f"          ,&jet1pt_f         ,"jet1pt_f/D"         );
	outTree_->Branch("jet1eta"          ,&jet1eta         ,"jet1eta/D"         );
	outTree_->Branch("jet1eta_f"          ,&jet1eta_f         ,"jet1eta_f/D"         );
	outTree_->Branch("jet1phi"          ,&jet1phi         ,"jet1phi/D"         );
	outTree_->Branch("jet1phi_f"          ,&jet1phi_f         ,"jet1phi_f/D"         );
	outTree_->Branch("jet1e"          ,&jet1e         ,"jet1e/D"         );
	outTree_->Branch("jet1e_f"          ,&jet1e_f         ,"jet1e_f/D"         );
	outTree_->Branch("jet1csv"          ,&jet1csv         ,"jet1csv/D"         );
	outTree_->Branch("jet1csv_f"          ,&jet1csv_f         ,"jet1csv_f/D"         );   
	outTree_->Branch("jet1icsv"          ,&jet1icsv         ,"jet1icsv/D"         );
	outTree_->Branch("jet1icsv_f"          ,&jet1icsv_f         ,"jet1icsv_f/D"         ); 
	outTree_->Branch("jet2pt"          ,&jet2pt         ,"jet2pt/D"         );
	outTree_->Branch("jet2pt_f"          ,&jet2pt_f         ,"jet2pt_f/D"         );
	outTree_->Branch("jet2eta"          ,&jet2eta         ,"jet2eta/D"         );
	outTree_->Branch("jet2eta_f"          ,&jet2eta_f         ,"jet2eta_f/D"         );
	outTree_->Branch("jet2phi"          ,&jet2phi         ,"jet2phi/D"         );
	outTree_->Branch("jet2phi_f"          ,&jet2phi_f         ,"jet2phi_f/D"         );
	outTree_->Branch("jet2e"          ,&jet2e         ,"jet2e/D"         );
	outTree_->Branch("jet2e_f"          ,&jet2e_f         ,"jet2e_f/D"         );
	outTree_->Branch("jet2csv"          ,&jet2csv         ,"jet2csv/D"         );
	outTree_->Branch("jet2csv_f"          ,&jet2csv_f         ,"jet2csv_f/D"         );
	outTree_->Branch("jet2icsv"          ,&jet2icsv         ,"jet2icsv/D"         );
	outTree_->Branch("jet2icsv_f"          ,&jet2icsv_f         ,"jet2icsv_f/D"         );
	outTree_->Branch("drj1a"          ,&drj1a         ,"drj1a/D"         );
	outTree_->Branch("drj1a_f"          ,&drj1a_f         ,"drj1a_f/D"         );
	outTree_->Branch("drj2a"          ,&drj2a         ,"drj2a/D"         );
	outTree_->Branch("drj2a_f"          ,&drj2a_f         ,"drj2a_f/D"         );
	outTree_->Branch("drj1l"          ,&drj1l         ,"drj1l/D"         );
	outTree_->Branch("drj1l_f"          ,&drj1l_f         ,"drj1l_f/D"         );
	outTree_->Branch("drj2l"          ,&drj2l         ,"drj2l/D"         );
	outTree_->Branch("drj2l_f"          ,&drj2l_f         ,"drj2l_f/D"         );
	outTree_->Branch("drj1l2"          ,&drj1l2         ,"drj1l2/D"         );
	outTree_->Branch("drj1l2_f"          ,&drj1l2_f         ,"drj1l2_f/D"         );
	outTree_->Branch("drj2l2"          ,&drj2l2         ,"drj2l2/D"         );
	outTree_->Branch("drj2l2_f"          ,&drj2l2_f         ,"drj2l2_f/D"         );
	outTree_->Branch("Mjj"          ,&Mjj         ,"Mjj/D"         );
	outTree_->Branch("Mjj_f"          ,&Mjj_f         ,"Mjj_f/D"         );
	outTree_->Branch("deltaetajj"          ,&deltaetajj         ,"deltaetajj/D"         );
	outTree_->Branch("deltaetajj_f"          ,&deltaetajj_f         ,"deltaetajj_f/D"         );
	outTree_->Branch("zepp"          ,&zepp         ,"zepp/D"         );
	outTree_->Branch("zepp_f"          ,&zepp_f         ,"zepp_f/D"         );
	outTree_->Branch("ptlep1"          ,&ptlep1         ,"ptlep1/D"         );
	outTree_->Branch("etalep1"         ,&etalep1        ,"etalep1/D"        );
	outTree_->Branch("philep1"         ,&philep1        ,"philep1/D"        );
	outTree_->Branch("ptlep2"          ,&ptlep2         ,"ptlep2/D"         );
	outTree_->Branch("etalep2"         ,&etalep2        ,"etalep2/D"        );
	outTree_->Branch("philep2"         ,&philep2        ,"philep2/D"        );
	//for muon rochester correction
	outTree_->Branch("muon1_trackerLayers"          ,&muon1_trackerLayers         ,"muon1_trackerLayers/I"         );
	outTree_->Branch("matchedgenMu1_pt"          ,&matchedgenMu1_pt         ,"matchedgenMu1_pt/D"         );
	outTree_->Branch("muon2_trackerLayers"          ,&muon2_trackerLayers         ,"muon2_trackerLayers/I"         );
	outTree_->Branch("matchedgenMu2_pt"          ,&matchedgenMu2_pt         ,"matchedgenMu2_pt/D"         );
	//for muon rochester correction
	outTree_->Branch("j1metPhi"          ,&j1metPhi         ,"j1metPhi/D"         );
	outTree_->Branch("j1metPhi_f"          ,&j1metPhi_f         ,"j1metPhi_f/D"         );
	outTree_->Branch("j2metPhi"          ,&j2metPhi         ,"j2metPhi/D"         );
	outTree_->Branch("j2metPhi_f"          ,&j2metPhi_f         ,"j2metPhi_f/D"         );
	// MET
	outTree_->Branch("MET_et",&MET_et,"MET_et/D");
	outTree_->Branch("MET_phi",&MET_phi,"MET_phi/D");
	//HLT bits
	outTree_->Branch("HLT_Ele1"  ,&HLT_Ele1 ,"HLT_Ele1/I" );
	outTree_->Branch("HLT_Ele2"  ,&HLT_Ele2 ,"HLT_Ele2/I" );
	outTree_->Branch("HLT_Ele3"  ,&HLT_Ele3 ,"HLT_Ele3/I" );
	outTree_->Branch("HLT_Ele4"  ,&HLT_Ele4 ,"HLT_Ele4/I" );
	outTree_->Branch("HLT_Ele5"  ,&HLT_Ele5 ,"HLT_Ele5/I" );
	outTree_->Branch("HLT_Ele6"  ,&HLT_Ele6 ,"HLT_Ele6/I" );
	outTree_->Branch("HLT_Ele7"  ,&HLT_Ele7 ,"HLT_Ele7/I" );
	outTree_->Branch("HLT_Mu1"   ,&HLT_Mu1  ,"HLT_Mu1/I"  );
	outTree_->Branch("HLT_Mu2"   ,&HLT_Mu2  ,"HLT_Mu2/I"  );
	outTree_->Branch("HLT_Mu3"   ,&HLT_Mu3  ,"HLT_Mu3/I"  );
	outTree_->Branch("HLT_Mu4"   ,&HLT_Mu4  ,"HLT_Mu4/I"  );
	outTree_->Branch("HLT_Mu5"   ,&HLT_Mu5  ,"HLT_Mu5/I"  );
	outTree_->Branch("HLT_Mu6"   ,&HLT_Mu6  ,"HLT_Mu6/I"  );
	outTree_->Branch("HLT_Mu7"   ,&HLT_Mu7  ,"HLT_Mu7/I"  );
	outTree_->Branch("HLT_Mu8"   ,&HLT_Mu8  ,"HLT_Mu8/I"  );
	outTree_->Branch("HLT_Mu9"   ,&HLT_Mu9  ,"HLT_Mu9/I"  );
        outTree_->Branch("HLT_Mu10"   ,&HLT_Mu10  ,"HLT_Mu10/I"  );
        outTree_->Branch("HLT_Mu11"   ,&HLT_Mu11  ,"HLT_Mu11/I"  );
	// filter
	outTree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
	outTree_->Branch("passFilter_HBHEIso"                 ,&passFilter_HBHEIso_                ,"passFilter_HBHEIso_/O");
	outTree_->Branch("passFilter_globalTightHalo"	,&passFilter_globalTightHalo_ ,"passFilter_globalTightHalo_/O");
	outTree_->Branch("passFilter_ECALDeadCell"         ,&passFilter_ECALDeadCell_        ,"passFilter_ECALDeadCell_/O");
	outTree_->Branch("passFilter_GoodVtx"              ,&passFilter_GoodVtx_             ,"passFilter_GoodVtx_/O");
	outTree_->Branch("passFilter_EEBadSc"              ,&passFilter_EEBadSc_             ,"passFilter_EEBadSc_/O");
	outTree_->Branch("passFilter_badMuon"                 ,&passFilter_badMuon_                ,"passFilter_badMuon_/O");
	outTree_->Branch("passFilter_badChargedHadron"                 ,&passFilter_badChargedHadron_                ,"passFilter_badChargedHadron_/O");

	// Meng, badmuon, duplicate muon
	outTree_->Branch("passFilter_MetbadMuon"		,&passFilter_MetbadMuon_		,"passFilter_MetbadMuon_/O");
	outTree_->Branch("passFilter_duplicateMuon"		,&passFilter_duplicateMuon_	,"passFilter_duplicateMuon_/O");
	outTree_->Branch("lumiWeight"      ,&lumiWeight     ,"lumiWeight/D"     );
	outTree_->Branch("pileupWeight"    ,&pileupWeight   ,"pileupWeight/D"   );

	// muon station2 retrieve, L1 issue, Meng 2017/3/26
	outTree_->Branch("lep1_eta_station2"      ,&lep1_eta_station2     ,"lep1_eta_station2/D"     ); 
	outTree_->Branch("lep1_phi_station2"      ,&lep1_phi_station2     ,"lep1_phi_station2/D"     );
	outTree_->Branch("lep1_sign"      ,&lep1_sign     ,"lep1_sign/I"     );
	outTree_->Branch("lep2_eta_station2"      ,&lep2_eta_station2     ,"lep2_eta_station2/D"     );
	outTree_->Branch("lep2_phi_station2"      ,&lep2_phi_station2     ,"lep2_phi_station2/D"     );
	outTree_->Branch("lep2_sign"      ,&lep2_sign     ,"lep2_sign/I"     );
	//L1 prefiring
	outTree_->Branch("prefiringweight", &prefiringweight, "prefiringweight/D");
	outTree_->Branch("prefiringweightup", &prefiringweightup, "prefiringweightup/D");
	outTree_->Branch("prefiringweightdown", &prefiringweightdown, "prefiringweightdown/D");
	//Lu
}

//------------------------------------
double ZPKUTreeMaker::getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
	double jetCorrFactor = 1.;
	if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
		jecAK4_->setJetEta( rawJetP4.eta() );
		jecAK4_->setJetPt ( rawJetP4.pt() );
		jecAK4_->setJetE  ( rawJetP4.energy() );
		jecAK4_->setJetPhi( rawJetP4.phi()    );
		jecAK4_->setJetA  ( jet.jetArea() );
		jecAK4_->setRho   ( *(rho_.product()) );
		jecAK4_->setNPV   ( nVtx );
		jetCorrFactor = jecAK4_->getCorrection();
	}
	reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
	corrJetP4 *= jetCorrFactor;
	return jetCorrFactor;
}
//------------------------------------
double ZPKUTreeMaker::getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
	double jetCorrFactor = 1.;
	if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
		jecOffset_->setJetEta( rawJetP4.eta()     );
		jecOffset_->setJetPt ( rawJetP4.pt()      );
		jecOffset_->setJetE  ( rawJetP4.energy()  );
		jecOffset_->setJetPhi( rawJetP4.phi()     );
		jecOffset_->setJetA  ( jet.jetArea()      );
		jecOffset_->setRho   ( *(rho_.product())  );
		jecOffset_->setNPV   ( nVtx  );
		jetCorrFactor = jecOffset_->getCorrection();
	}
	reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
	corrJetP4 *= jetCorrFactor;
	return jetCorrFactor;
}
//------------------------------------
void ZPKUTreeMaker::addTypeICorr( edm::Event const & event ){
	TypeICorrMap_.clear();
	edm::Handle<pat::JetCollection> jets_;
	event.getByToken(t1jetSrc_, jets_);
	event.getByToken(rhoToken_      , rho_     );
	edm::Handle<reco::VertexCollection> vertices_;
	event.getByToken(VertexToken_, vertices_);
	edm::Handle<edm::View<pat::Muon>> muons_;
	event.getByToken(t1muSrc_,muons_);
	bool skipEM_                    = true;
	double skipEMfractionThreshold_ = 0.9;
	bool skipMuons_                 = true;
	std::string skipMuonSelection_string = "isGlobalMuon | isStandAloneMuon";
	StringCutObjectSelector<reco::Candidate>* skipMuonSelection_ = new StringCutObjectSelector<reco::Candidate>(skipMuonSelection_string,true);
	double jetCorrEtaMax_           = 9.9;
	double type1JetPtThreshold_     = 15.0;
	double corrEx    = 0;
	double corrEy    = 0;
	double corrSumEt = 0;
	std::vector<JetCorrectorParameters> vPar;
	for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4chsLabels_.begin(), payloadEnd = jecAK4chsLabels_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters pars(*ipayload);
		vPar.push_back(pars);
	}
	jecAK4_ = new FactorizedJetCorrector(vPar);
	vPar.clear();
	for ( std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters pars(*ipayload);
		vPar.push_back(pars);
	}
	jecOffset_ = new FactorizedJetCorrector(vPar);
	vPar.clear();
	for (const pat::Jet &jet : *jets_) {
		double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
		if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;
		reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);
		double corr = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);
		
		   if ( skipMuons_ && jet.muonMultiplicity() != 0 ) {
		   for (const pat::Muon &muon : *muons_) {
		   if( !muon.isGlobalMuon() && !muon.isStandAloneMuon() ) continue;
		   TLorentzVector muonV; muonV.SetPtEtaPhiE(muon.p4().pt(),muon.p4().eta(),muon.p4().phi(),muon.p4().e());
		   TLorentzVector jetV; jetV.SetPtEtaPhiE(jet.p4().pt(),jet.p4().eta(),jet.p4().phi(),jet.p4().e());
		   if( muonV.DeltaR(jetV) < 0.5 ){
		   reco::Candidate::LorentzVector muonP4 = muon.p4();
		   rawJetP4 -= muonP4;
		   }
		   }
		   }
		   
		if ( skipMuons_ ) {
			const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
			for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
					cand != cands.end(); ++cand ) {
				const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
				const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
				if ( mu != 0 && (*skipMuonSelection_)(*mu) ) {
					reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
					rawJetP4 -= muonP4;
				}
			}
		}

		reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;
		if ( corrJetP4.pt() > type1JetPtThreshold_ ) {
			reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
			corr = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
			reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;
			corrEx    -= (corrJetP4.px() - rawJetP4offsetCorr.px());
			corrEy    -= (corrJetP4.py() - rawJetP4offsetCorr.py());
			corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
		}
	}
	TypeICorrMap_["corrEx"]    = corrEx;
	TypeICorrMap_["corrEy"]    = corrEy;
	TypeICorrMap_["corrSumEt"] = corrSumEt;
	delete jecAK4_;
	jecAK4_=0;
	delete jecOffset_;
	jecOffset_=0;
	delete skipMuonSelection_;
	skipMuonSelection_=0;
}
//------------------------------------
bool ZPKUTreeMaker::hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<edm::View<pat::Electron> > &eleCol, const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,  float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax) {
	//check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
	//and not matching any conversion in the collection passing the quality cuts
	if (sc.isNull()) return false;
	for (edm::View<pat::Electron>::const_iterator it = eleCol->begin(); it!=eleCol->end(); ++it) {
		//match electron to supercluster
		if (it->superCluster()!=sc) continue;
		//check expected inner hits
		if (it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
		//check if electron is matching to a conversion
		if (ConversionTools::hasMatchedConversion(*it,convCol,beamspot)) continue;
		return true;
	}
	return false;
}
//------------------------------------
int ZPKUTreeMaker::matchToTruth(const reco::Photon &pho,
		const edm::Handle<edm::View<reco::GenParticle>>
		&genParticles, bool &ISRPho, double &dR, int &isprompt)
{
	//
	// Explicit loop and geometric matching method
	// Find the closest status 1 gen photon to the reco photon
	dR = 999;
	const reco::Candidate *closestPhoton = 0;
	int im=0;
	for(size_t i=0; i<genParticles->size();i++){
		const reco::Candidate *particle = &(*genParticles)[i];
		if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
			continue;
		double dRtmp = deltaR(pho.eta(),pho.phi(),particle->eta(),particle->phi());  
		if( dRtmp < dR ){
			dR = dRtmp;
			im=i;    
			closestPhoton = particle;
		}
	}
	// See if the closest photon (if it exists) is close enough.
	// If not, no match found.
	if( !(closestPhoton != 0 && dR < 0.3) ) {
		return UNMATCHED;
		// ISRPho = false;
	}
	isprompt=(*genParticles)[im].isPromptFinalState();
	// Find ID of the parent of the found generator level photon match
	int ancestorPID = -999;
	int ancestorStatus = -999;
	findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);
	// Allowed parens: quarks pdgId 1-5, or a gluon 21
	std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, 21,-11,11,-13,13,-15,15,23,-24,24 };
	if( !(std::find(allowedParents.begin(),
					allowedParents.end(), ancestorPID)
				!= allowedParents.end()) ){
		// So it is not from g, u, d, s, c, b. Check if it is from pi0 or not.
		if( abs(ancestorPID) == 111 )
			return MATCHED_FROM_PI0;
		// ISRPho =true;
		else
			//std::cout<<"Mother = "<<abs(ancestorPID)<<std::endl;
			return MATCHED_FROM_OTHER_SOURCES;
		//  ISRPho =true;
	}
	return MATCHED_FROM_GUDSCB;
	//   ISRPho =true;
}
//------------------------------------
void ZPKUTreeMaker::findFirstNonPhotonMother(const reco::Candidate *particle,
		int &ancestorPID, int &ancestorStatus){
	if( particle == 0 ){
		printf("SimplePhotonNtupler: ERROR! null candidate pointer, this should never happen\n");
		return;
	}
	// Is this the first non-photon parent? If yes, return, otherwise
	// go deeper into recursion
	if( abs(particle->pdgId()) == 22 ){
		findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
	}else{
		ancestorPID = particle->pdgId();
		ancestorStatus = particle->status();
	}
	return;
}
//------------------------------------
//------------------------------------
ZPKUTreeMaker::~ZPKUTreeMaker()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}
//------------------------------------
//------------------------------------
	void
ZPKUTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	edm::ParameterSet  vDefaults1;
	vDefaults1.addParameter<std::string>("useTrack", "tracker");
	vDefaults1.addParameter<std::string>("useState", "atVertex");
	vDefaults1.addParameter<bool>("useSimpleGeometry", true);
	vDefaults1.addParameter<bool>("useStation2",true );
	vDefaults1.addParameter<bool>("fallbackToME1", true);
	/*PropagateToMuon **/muPropagator2nd_ = new PropagateToMuon(vDefaults1);
	muPropagator2nd_->init(iSetup);
	setDummyValues(); //Initalize variables with dummy values
	nevent = iEvent.eventAuxiliary().event();
	run    = iEvent.eventAuxiliary().run();
	ls     = iEvent.eventAuxiliary().luminosityBlock();
	//events weight
	if(!RunOnMC_) hout->Fill(0.5);
	if (RunOnMC_){
		//L1 prefiring
		edm::Handle< double > theprefweight;
		iEvent.getByToken(prefweight_token, theprefweight ) ;
		prefiringweight =(*theprefweight);
		
		edm::Handle< double > theprefweightup;
		iEvent.getByToken(prefweightup_token, theprefweightup ) ;
		prefiringweightup =(*theprefweightup);
		
		edm::Handle< double > theprefweightdown;
		iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
		prefiringweightdown =(*theprefweightdown);

		edm::Handle<GenEventInfoProduct> genEvtInfo;
		iEvent.getByToken(GenToken_,genEvtInfo);
		theWeight = genEvtInfo->weight();
		if(theWeight>0) {nump = nump+1;hout->Fill(0.5);}
		if(theWeight<0) {numm = numm+1;hout->Fill(-0.5);}
	
		edm::Handle<LHEEventProduct> wgtsource;
		// cout the weight name 2018.07.03
                iEvent.getByToken(LheToken_, wgtsource);
//		std::cout<<"weight number "<<wgtsource->weights().size()<<std::endl;
//                for ( int i=0; i<703; ++i) {
//                     pweight[i]= wgtsource->weights()[i].wgt/wgtsource->originalXWGTUP();
//                }

		edm::Handle<std::vector<PileupSummaryInfo>>  PupInfo;
		iEvent.getByToken(PUToken_, PupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			nBX = PVI->getBunchCrossing();
			if(nBX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
				npT = PVI->getTrueNumInteractions();
				npIT = PVI->getPU_NumInteractions();
			}
		} 
	} 

	Handle<TriggerResults> trigRes;
	iEvent.getByToken(hltToken_, trigRes);
	int xtemp1=0;
	for (size_t i=0; i<elPaths1.size();i++) {
		xtemp1 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths1[i]));
		if(HLT_Ele1<xtemp1) HLT_Ele1=xtemp1;
	}
	int xtemp2=0;
	for (size_t i=0; i<elPaths2.size();i++) {
		xtemp2 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths2[i]));
		if(HLT_Ele2<xtemp2) HLT_Ele2=xtemp2;
	}
	int xtemp3=0;
        for (size_t i=0; i<elPaths3.size();i++) {
                xtemp3 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths3[i]));
                if(HLT_Ele3<xtemp3) HLT_Ele3=xtemp3;
        }
        int xtemp4=0;
        for (size_t i=0; i<elPaths4.size();i++) {
                xtemp4 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths4[i]));
                if(HLT_Ele4<xtemp4) HLT_Ele4=xtemp4;
        }
	int xtemp5=0;
        for (size_t i=0; i<elPaths5.size();i++) {
                xtemp5 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths5[i]));
                if(HLT_Ele5<xtemp5) HLT_Ele5=xtemp5;
        }
        int xtemp6=0;
        for (size_t i=0; i<elPaths6.size();i++) {
                xtemp6 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths6[i]));
                if(HLT_Ele6<xtemp6) HLT_Ele6=xtemp6;
        }
	int xtemp7=0;
        for (size_t i=0; i<elPaths7.size();i++) {
                xtemp7 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths7[i]));
                if(HLT_Ele7<xtemp7) HLT_Ele7=xtemp7;
        }
	int mtemp1=0;
	for (size_t i=0; i<muPaths1.size();i++) {
		mtemp1 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths1[i]));
		if(HLT_Mu1<mtemp1) HLT_Mu1=mtemp1;
	}
	int mtemp2=0;
	for (size_t i=0; i<muPaths2.size();i++) {
		mtemp2 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths2[i]));
		if(HLT_Mu2<mtemp2) HLT_Mu2=mtemp2;
	}
	int mtemp3=0;
	for (size_t i=0; i<muPaths3.size();i++) {
		mtemp3 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths3[i]));
		if(HLT_Mu3<mtemp3) HLT_Mu3=mtemp3;
	}
	int mtemp4=0;
	for (size_t i=0; i<muPaths4.size();i++) {
		mtemp4 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths4[i]));
		if(HLT_Mu4<mtemp4) HLT_Mu4=mtemp4;
	}
	int mtemp5=0;
	for (size_t i=0; i<muPaths5.size();i++) {
		mtemp5 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths5[i]));
		if(HLT_Mu5<mtemp5) HLT_Mu5=mtemp5;
	}
	int mtemp6=0;
	for (size_t i=0; i<muPaths6.size();i++) {
		mtemp6 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths6[i]));
		if(HLT_Mu6<mtemp6) HLT_Mu6=mtemp6;
	}
	int mtemp7=0;
	for (size_t i=0; i<muPaths7.size();i++) {
		mtemp7 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths7[i]));
		if(HLT_Mu7<mtemp7) HLT_Mu7=mtemp7;
	}
	int mtemp8=0;
	for (size_t i=0; i<muPaths8.size();i++) {
		mtemp8 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths8[i]));
		if(HLT_Mu8<mtemp8) HLT_Mu8=mtemp8;
	}
	int mtemp9=0;
        for (size_t i=0; i<muPaths9.size();i++) {
                mtemp9 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths9[i]));
                if(HLT_Mu9<mtemp9) HLT_Mu9=mtemp9;
        }
        int mtemp10=0;
        for (size_t i=0; i<muPaths10.size();i++) {
                mtemp10 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths10[i]));
                if(HLT_Mu10<mtemp10) HLT_Mu10=mtemp10;
        }
        int mtemp11=0;
        for (size_t i=0; i<muPaths11.size();i++) {
                mtemp11 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths11[i]));
                if(HLT_Mu11<mtemp11) HLT_Mu11=mtemp11;
        }

	edm::Handle<edm::View<reco::GenParticle> > genParticles; 
	iEvent.getByToken(genSrc_, genParticles);
	//   iEvent.getByLabel(InputTag("packedGenParticles"), genParticles);

	if (RunOnMC_){
		int ipp=0, imm=0, iee=0;
		for(size_t i=0; i<genParticles->size();i++){
			const reco::Candidate *particle = &(*genParticles)[i];
			if( abs(particle->pdgId())== 22 && particle->status()== 1 && (*genParticles)[i].isPromptFinalState()>0 && ipp<6 ) {
				genphoton_pt[ipp]=particle->pt();
				genphoton_eta[ipp]=particle->eta();
				genphoton_phi[ipp]=particle->phi();
				ipp++;
			}
			if( abs(particle->pdgId())== 13 && particle->status()== 1 && (*genParticles)[i].isPromptFinalState()>0 && imm<6 ) {
				genmuon_pt[imm]=particle->pt();
				genmuon_eta[imm]=particle->eta();
				genmuon_phi[imm]=particle->phi();
				genmuon_pid[imm]=particle->pdgId();
				imm++;
			}
			if( abs(particle->pdgId())== 11 && particle->status()== 1 && (*genParticles)[i].isPromptFinalState()>0 && iee<6 ) {
				genelectron_pt[iee]=particle->pt();
				genelectron_eta[iee]=particle->eta();
				genelectron_phi[iee]=particle->phi();
				iee++;
			}
		}
	}

	if(RunOnMC_){
		int ijj=0;
		edm::Handle<reco::GenJetCollection> genJets;
		iEvent.getByToken(genJet_,genJets);
		reco::GenJetCollection::const_iterator i_jet;
//		std::cout<<"nevent "<<nevent<<std::endl;//gen jet debug
		for( i_jet=genJets->begin(); i_jet != genJets->end();i_jet++){
			if(ijj>5) continue;
			genjet_pt[ijj]= i_jet->pt();
//			std::cout<<"pt "<<genjet_pt[ijj]<<std::endl;//gen jet debug
			genjet_eta[ijj]= i_jet->eta();
//			std::cout<<"eta "<<genjet_eta[ijj]<<std::endl;//gen jet debug
			genjet_phi[ijj]=i_jet->phi();
//			std::cout<<"phi "<<genjet_phi[ijj]<<std::endl;//gen jet debug
			genjet_e[ijj] = i_jet->energy();
//			std::cout<<"e "<<genjet_e[ijj]<<std::endl;//gen jet debug
			ijj++;
		}
	}

	edm::Handle<edm::View<reco::Candidate> > leptonicVs;
        iEvent.getByToken(leptonicVSrc_, leptonicVs);
        if (leptonicVs->empty()) { /* outTree_->Fill();*/ return;  }

        iEvent.getByToken(rhoToken_      , rho_     );
        double fastJetRho = *(rho_.product());
        useless = fastJetRho;

        edm::Handle<edm::View<pat::Jet> > ak4jets;
        iEvent.getByToken(ak4jetsSrc_, ak4jets);

        edm::Handle<edm::View<pat::Photon> > photons;
        iEvent.getByToken(photonSrc_, photons);

	edm::Handle<edm::View<pat::Muon>> loosemus;
	iEvent.getByToken(loosemuonToken_,loosemus); 
	edm::Handle<edm::View<pat::Muon>> goodmus;
	iEvent.getByToken(goodmuonToken_,goodmus);
	edm::Handle<edm::View<pat::Electron>> looseeles;
	iEvent.getByToken(looseelectronToken_,looseeles); 
	edm::Handle<edm::View<reco::Candidate> > metHandle; 
	iEvent.getByToken(metSrc_, metHandle); 
	//filter
	iEvent.getByToken(noiseFilterToken_, noiseFilterBits_);
	const edm::TriggerNames &names = iEvent.triggerNames(*noiseFilterBits_);
	for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
		if (names.triggerName(i) == HBHENoiseFilter_Selector_)
			passFilter_HBHE_ = noiseFilterBits_->accept(i);  
		if (names.triggerName(i) == HBHENoiseIsoFilter_Selector_)
			passFilter_HBHEIso_ = noiseFilterBits_->accept(i);  
		if (names.triggerName(i) ==ECALDeadCellNoiseFilter_Selector_)
			passFilter_globalTightHalo_ = noiseFilterBits_->accept(i);  
		if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
			passFilter_ECALDeadCell_ = noiseFilterBits_->accept(i);  
		if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_)
			passFilter_GoodVtx_ = noiseFilterBits_->accept(i);  
		if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
			passFilter_EEBadSc_ = noiseFilterBits_->accept(i);  
		// Meng
		if (names.triggerName(i) == badMuonFilter_Selector_)
			passFilter_MetbadMuon_ = noiseFilterBits_->accept(i);
		if (names.triggerName(i) == duplicateMuonFilter_Selector_)
			passFilter_duplicateMuon_ = noiseFilterBits_->accept(i);
	}
	edm::Handle<bool> badMuonResultHandle;
	edm::Handle<bool> badChargedHadronResultHandle;
	iEvent.getByToken(badMuon_Selector_, badMuonResultHandle);
	iEvent.getByToken(badChargedHadron_Selector_, badChargedHadronResultHandle);
	passFilter_badMuon_ = *badMuonResultHandle;
	passFilter_badChargedHadron_ = *badChargedHadronResultHandle;

	const reco::Candidate& leptonicV = leptonicVs->at(0);
	const reco::Candidate& metCand = metHandle->at(0);

	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(VertexToken_, vertices);
	if (vertices->empty()) { outTree_->Fill(); return;} // skip the event if no PV found
	nVtx = vertices->size();
	reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
	for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
		// Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
		// Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
		if (  // !vtx->isFake() &&
				!(vtx->chi2()==0 && vtx->ndof()==0) 
				&&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
				&& fabs(vtx->position().Z())<=24.0) {
			firstGoodVertex = vtx;
			break;
		}           
	}
	if ( firstGoodVertex==vertices->end() ) {outTree_->Fill();  return;} // skip event if there are no good PVs


	// ************************* MET ********************** //
	edm::Handle<pat::METCollection>  METs_;
	bool defaultMET = iEvent.getByToken(metInputToken_ , METs_ );
	if(RunOnMC_){
		const pat::MET &xmet = METs_->front();
		genMET=xmet.genMET()->pt();
	}
	if(defaultMET){
		addTypeICorr(iEvent);
//		addTypeICorr_user(iEvent);
		for (const pat::MET &met : *METs_) {
			//         const float  rawPt    = met.shiftedPt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
			//         const float  rawPhi   = met.shiftedPhi(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
			//         const float  rawSumEt = met.shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
			const float rawPt = met.uncorPt();
			const float rawPhi = met.uncorPhi();
			const float rawSumEt = met.uncorSumEt();
			TVector2 rawMET_;
			rawMET_.SetMagPhi (rawPt, rawPhi );
			Double_t rawPx = rawMET_.Px();
			Double_t rawPy = rawMET_.Py();
			Double_t rawEt = std::hypot(rawPx,rawPy);
			METraw_et = rawEt;
			METraw_phi = rawPhi;
			METraw_sumEt = rawSumEt;

			double pxcorr = rawPx+TypeICorrMap_["corrEx"];
			double pycorr = rawPy+TypeICorrMap_["corrEy"];
			double et     = std::hypot(pxcorr,pycorr);
			double sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];
			TLorentzVector corrmet;
			corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
			useless = sumEtcorr;
			useless = rawEt;
			MET_et = et;
                        MET_phi = corrmet.Phi();
                        MET_sumEt = sumEtcorr;
                        MET_corrPx = TypeICorrMap_["corrEx"];
                        MET_corrPy = TypeICorrMap_["corrEy"];
		}
	}
	//------------------------------------
	/// For the time being, set these to 1
	triggerWeight=1.0;
	pileupWeight=1.0;
	double targetEvents = targetLumiInvPb_*crossSectionPb_;
	lumiWeight = targetEvents/originalNEvents_;
	lep          = std::max(abs(leptonicV.daughter(0)->pdgId()), abs(leptonicV.daughter(1)->pdgId()));
	ptVlep       = leptonicV.pt();
	yVlep        = leptonicV.eta();
	phiVlep      = leptonicV.phi();
	massVlep     = leptonicV.mass();
	// muon station2 retrieve, L1 issue, Meng 2017/3/26
	if(goodmus->size()>1){
		lep1_etaphi_ = EtaPhiAtME2X(&((*goodmus)[0]), muPropagator2nd_);
		lep1_eta_station2 = lep1_etaphi_.first;
		lep1_phi_station2 = lep1_etaphi_.second;
		lep1_sign = leptonicV.daughter(0)->pdgId();
		lep2_etaphi_ = EtaPhiAtME2X(&((*goodmus)[1]), muPropagator2nd_); 
		lep2_eta_station2 = lep2_etaphi_.first;
		lep2_phi_station2 = lep2_etaphi_.second;
		lep2_sign = leptonicV.daughter(1)->pdgId();
	}
	// Lu
	ptlep1       = leptonicV.daughter(0)->pt();
//	std::cout<<"ptlep1 "<<ptlep1<<std::endl;
	etalep1      = leptonicV.daughter(0)->eta();
	philep1      = leptonicV.daughter(0)->phi(); 
	ptlep2       = leptonicV.daughter(1)->pt();
//	std::cout<<"ptlep2 "<<ptlep2<<std::endl;
	etalep2      = leptonicV.daughter(1)->eta();
	philep2      = leptonicV.daughter(1)->phi();
	// for muon rochester correction
	if(goodmus->size()>1){
		muon1_trackerLayers      = (*goodmus)[0].innerTrack()->hitPattern().trackerLayersWithMeasurement(); 
		muon2_trackerLayers      = (*goodmus)[1].innerTrack()->hitPattern().trackerLayersWithMeasurement(); 
	}
	if(lep==13)
	{	
		double dr_temp=1e2;
		double matchmuon_pt = -1e2;
		for (int i=0;i<6;i++)
		{
			if(lep1_sign==genmuon_pid[i]) 
			{
				double dr_mugenmu=getDR(etalep1,philep1,genmuon_eta[i],genmuon_phi[i]);
				if(dr_mugenmu<dr_temp) {matchmuon_pt=genmuon_pt[i];dr_temp=dr_mugenmu;}
			}
		}
		if(dr_temp<0.3) matchedgenMu1_pt=matchmuon_pt;
	}
//	std::cout<<"matchedgenMu1_pt "<<matchedgenMu1_pt<<std::endl;
	if(lep==13)
        {
                double dr_temp=1e2;
                double matchmuon_pt = -1e2;
                for (int i=0;i<6;i++)
                {
                        if(lep2_sign==genmuon_pid[i])
                        {
                                double dr_mugenmu=getDR(etalep2,philep2,genmuon_eta[i],genmuon_phi[i]);
                                if(dr_mugenmu<dr_temp) {matchmuon_pt=genmuon_pt[i];dr_temp=dr_mugenmu;}
                        }
                }
                if(dr_temp<0.3) matchedgenMu2_pt=matchmuon_pt;
        }
//	std::cout<<"matchedgenMu2_pt "<<matchedgenMu2_pt<<std::endl;
	// for muon rochester correction
	double energylep1     = leptonicV.daughter(0)->energy();
	double energylep2     = leptonicV.daughter(1)->energy();
	met          = metCand.pt();
	metPhi       = metCand.phi();
	nlooseeles = looseeles->size(); 
	nloosemus = loosemus->size(); 
	ngoodmus = goodmus->size();
	TLorentzVector  glepton;
	glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);
	TLorentzVector  glepton2;
	glepton2.SetPtEtaPhiE(ptlep2, etalep2, philep2, energylep2);

	// ************************* Photon Jets Information****************** //
	// *************************************************************//
	double rhoVal_;
	rhoVal_=-99.;
	rhoVal_ = *rho_;
	edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
	iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
	edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
	iEvent.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
	edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
	iEvent.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
	edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
	iEvent.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);

	photonet=-100.; photonet_f=-100.;  iphoton=-1; iphoton_f=-1;
	for (size_t ip=0; ip<photons->size();ip++)
	{
		const auto pho = photons->ptrAt(ip);

		double phosc_eta=pho->superCluster()->eta();
//		double phosc_phi=pho->superCluster()->phi();
		double pho_ieie=(*photons)[ip].full5x5_sigmaIetaIeta();//(*full5x5SigmaIEtaIEtaMap)[pho];
		double chIso1 = (*photons)[ip].chargedHadronIso();// (*phoChargedIsolationMap)[pho];
		double nhIso1 = (*photons)[ip].neutralHadronIso();// (*phoNeutralHadronIsolationMap)[pho];
		double phIso1 = (*photons)[ip].photonIso();//(*phoPhotonIsolationMap)[pho];
		double chiso=std::max(0.0, chIso1 - rhoVal_*EAch(fabs((*photons)[ip].eta()))); //effAreaChHadrons_.getEffectiveArea(fabs(phosc_eta)));
		//            double chiso=std::max((*photons)[ip].chargedHadronIso()-rhoVal_*EAch(fabs((*photons)[ip].eta())),0.0);
		double nhiso=std::max(0.0, nhIso1 - rhoVal_*EAnh(fabs((*photons)[ip].eta()))); //effAreaNeuHadrons_.getEffectiveArea(fabs(phosc_eta)));
		//            double nhiso=std::max((*photons)[ip].neutralHadronIso()-rhoVal_*EAnh(fabs((*photons)[ip].eta())),0.0);
		double phoiso=std::max(0.0, phIso1 - rhoVal_*EApho(fabs((*photons)[ip].eta()))); //effAreaPhotons_.getEffectiveArea(fabs(phosc_eta)));
		//            double phoiso=std::max((*photons)[ip].photonIso()-rhoVal_*EApho(fabs((*photons)[ip].eta())),0.0);

		int ismedium_photon=0;
		int ismedium_photon_f=0;
		edm::Handle<edm::View<pat::Electron> > electrons;
		iEvent.getByToken(electronToken_, electrons);
		edm::Handle<reco::BeamSpot> beamSpot;
		iEvent.getByToken(beamSpotToken_,beamSpot);
		edm::Handle<std::vector<reco::Conversion> > conversions;
		iEvent.getByToken(conversionsToken_,conversions);
		passEleVeto = (!hasMatchedPromptElectron((*photons)[ip].superCluster(),electrons, conversions, beamSpot->position() ) );

		passEleVetonew=(*photons)[ip].passElectronVeto();
		passPixelSeedVeto=(*photons)[ip].hasPixelSeed();


		if(ip<6)  {
			photon_pt[ip] = (*photons)[ip].pt();
			photon_eta[ip] = (*photons)[ip].eta();
			photon_sc_eta[ip] =phosc_eta ;
			photon_phi[ip] = (*photons)[ip].phi();
			photon_e[ip] = (*photons)[ip].energy();
			photon_pev[ip]=passEleVeto;
			photon_pevnew[ip]=passEleVetonew;
			photon_ppsv[ip]=passPixelSeedVeto;
			photon_iseb[ip]=(*photons)[ip].isEB();
			photon_isee[ip]=(*photons)[ip].isEE();
			photon_hoe[ip]=(*photons)[ip].hadTowOverEm();
			photon_sieie[ip]=pho_ieie;//(*photons)[ip].sigmaIetaIeta();
			photon_sieie2[ip]=(*photons)[ip].sigmaIetaIeta();
			photon_chiso[ip]=chiso;
			photon_nhiso[ip]=nhiso;
			photon_phoiso[ip]=phoiso;
			if(RunOnMC_ && photon_pt[ip]>0){
				const auto pho = photons->ptrAt(ip);
				photon_istrue[ip]=matchToTruth(*pho, genParticles, ISRPho, dR_, photon_isprompt[ip]);
			}
			photon_drla[ip]=deltaR(photon_eta[ip],photon_phi[ip],etalep1,philep1);
			photon_drla2[ip]=deltaR(photon_eta[ip],photon_phi[ip],etalep2,philep2);
			TLorentzVector tp4;
			tp4.SetPtEtaPhiE(photon_pt[ip],photon_eta[ip],photon_phi[ip],photon_e[ip]);
			photon_mla[ip]=(tp4+glepton).M();
			photon_mla2[ip]=(tp4+glepton2).M();
			photon_mva[ip]=(tp4+glepton+glepton2).M();
		}


		if(passEleVetonew && fabs(phosc_eta)<1.4442 && (*photons)[ip].hadTowOverEm()<0.0396 && photon_sieie[ip]<0.01022 && chiso<0.441 && nhiso<(2.725 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(2.571+0.0047*(*photons)[ip].pt())) {ismedium_photon=1;}
		if(passEleVetonew && fabs(phosc_eta)<2.5 &&fabs(phosc_eta)>1.566 && (*photons)[ip].hadTowOverEm()<0.0219 && photon_sieie[ip]<0.03001 && chiso<0.442 && nhiso<(1.715 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(3.863+0.0034*(*photons)[ip].pt())) {ismedium_photon=1;}

		if(ismedium_photon==1 && deltaR(photon_eta[ip],photon_phi[ip],etalep1,philep1) > 0.5 && deltaR(photon_eta[ip],photon_phi[ip],etalep2,philep2) > 0.5) { 
			if(ip==0) {photonet=(*photons)[ip].pt(); iphoton=ip;}
			if((*photons)[ip].pt()>photonet)
			{
				photonet=(*photons)[ip].pt(); iphoton=ip;
			}
		}

		//////////////////////////////////for fake photon study, store photon without sieie cut
		////Inverting loose ID
		if(passEleVetonew && fabs(phosc_eta)<1.4442 && (*photons)[ip].hadTowOverEm()<0.0597 && chiso<std::min(0.2*(*photons)[ip].pt(), 5.*1.295) && nhiso<std::min(0.2*(*photons)[ip].pt(), 5.*(10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())) ) && phoiso<std::min(0.2*(*photons)[ip].pt(), 5.*(3.630+0.0047*(*photons)[ip].pt()))  && !(photon_sieie[ip]<0.01031&&chiso<1.295 && nhiso<(10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(3.630+0.0047*(*photons)[ip].pt()))) {ismedium_photon_f=1;}
		if(passEleVetonew && fabs(phosc_eta)<2.5 &&fabs(phosc_eta)>1.566 && (*photons)[ip].hadTowOverEm()<0.0481 && chiso<std::min(0.2*(*photons)[ip].pt(), 5.*1.011) && nhiso<std::min(0.2*(*photons)[ip].pt(), 5.*(5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())) ) && phoiso<std::min(0.2*(*photons)[ip].pt(), 5.*(6.641+0.0034*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.03013 && chiso<1.011 && nhiso<(5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(6.641+0.0034*(*photons)[ip].pt()))) {ismedium_photon_f=1;}
		if(ismedium_photon_f==1 && deltaR(photon_eta[ip],photon_phi[ip],etalep1,philep1) > 0.5 && deltaR(photon_eta[ip],photon_phi[ip],etalep2,philep2) > 0.5) { 
			if(ip==0) {photonet_f=(*photons)[ip].pt(); iphoton_f=ip;}
			if((*photons)[ip].pt()>photonet_f) {
				photonet_f=(*photons)[ip].pt(); iphoton_f=ip;
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////
	}

	//Gen photon matching
	if(RunOnMC_ && iphoton>-1){
		const auto pho1 = photons->ptrAt(iphoton);
		isTrue_= matchToTruth(*pho1, genParticles, ISRPho, dR_, isprompt_);
	}

	if(iphoton>-1 && iphoton<6) {
		photonet=photon_pt[iphoton];//(*photons)[iphoton].pt();
		photoneta=photon_eta[iphoton];//(*photons)[iphoton].eta();
		photoneta_sc=photon_sc_eta[iphoton];//(*photons)[iphoton].eta();
		photonphi=photon_phi[iphoton];//(*photons)[iphoton].phi();
		photone=photon_e[iphoton];//(*photons)[iphoton].energy();
		photonsieie=photon_sieie[iphoton];//(*photons)[iphoton].sigmaIetaIeta();
		photonphoiso=photon_phoiso[iphoton];//std::max((*photons)[iphoton].photonIso()-rhoVal_*EApho(fabs((*photons)[iphoton].eta())),0.0);
		photonchiso=photon_chiso[iphoton];//std::max((*photons)[iphoton].chargedHadronIso()-rhoVal_*EAch(fabs((*photons)[iphoton].eta())),0.0);
		photonnhiso=photon_nhiso[iphoton];//std::max((*photons)[iphoton].neutralHadronIso()-rhoVal_*EAnh(fabs((*photons)[iphoton].eta())),0.0);
		drla=deltaR(photon_eta[iphoton],photon_phi[iphoton],etalep1,philep1);
		drla2=deltaR(photon_eta[iphoton],photon_phi[iphoton],etalep2,philep2);
		TLorentzVector photonp4;
		photonp4.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
		Mla=(photonp4+glepton).M();
		Mla2=(photonp4+glepton2).M();
		TLorentzVector wp4;
		wp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
		Mva=(photonp4+wp4).M();
	}

	if(iphoton_f>-1 && iphoton_f<6) {
		photonet_f=photon_pt[iphoton_f];
		photoneta_f=photon_eta[iphoton_f];
		photoneta_sc_f=photon_sc_eta[iphoton_f];
		photonphi_f=photon_phi[iphoton_f];
		photone_f=photon_e[iphoton_f];
		photonsieie_f=photon_sieie[iphoton_f];
		photonphoiso_f=photon_phoiso[iphoton_f];
		photonchiso_f=photon_chiso[iphoton_f];
		photonnhiso_f=photon_nhiso[iphoton_f];
		drla_f=deltaR(photon_eta[iphoton_f],photon_phi[iphoton_f],etalep1,philep1);
		drla2_f=deltaR(photon_eta[iphoton_f],photon_phi[iphoton_f],etalep2,philep2);
		TLorentzVector photonp4_f;
		photonp4_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
		Mla_f=(photonp4_f+glepton).M();
		Mla2_f=(photonp4_f+glepton2).M();
		TLorentzVector wp4_f;
		wp4_f.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
		Mva_f=(photonp4_f+wp4_f).M();
	}

	// ************************* AK4 Jets Information****************** //
	// ***********************************************************//
	Int_t jetindexphoton12[2] = {-1,-1}; 
	Int_t jetindexphoton12_f[2] = {-1,-1};

	std::vector<JetCorrectorParameters> vPar;
	for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4Labels_.begin(), payloadEnd = jecAK4Labels_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters pars(*ipayload);
		vPar.push_back(pars);    }
	jecAK4_ = new FactorizedJetCorrector(vPar);
	vPar.clear();

	int nujets=0 ;
	double tmpjetptcut=20.0;
	std::vector<TLorentzVector*> jets;

	//################Jet Correction##########################
	for (size_t ik=0; ik<ak4jets->size();ik++)
	{
		reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
		jecAK4_->setJetEta( uncorrJet.eta() );
		jecAK4_->setJetPt ( uncorrJet.pt() );
		jecAK4_->setJetE ( uncorrJet.energy() );
		jecAK4_->setRho ( rhoVal_ );
		jecAK4_->setNPV ( vertices->size() );
		jecAK4_->setJetA ( (*ak4jets)[ik].jetArea() );
		double corr = jecAK4_->getCorrection();
		if(corr*uncorrJet.pt()>tmpjetptcut) {
                        TLorentzVector *dummy = new TLorentzVector(0,0,0,0);
                        dummy->SetPtEtaPhiE(corr*uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), corr*uncorrJet.energy());
                        jets.push_back(dummy);
                        ++nujets;
                }
                if(ik<6)  {
                        ak4jet_pt[ik] =  corr*uncorrJet.pt();
                        ak4jet_eta[ik] = (*ak4jets)[ik].eta();
                        ak4jet_phi[ik] = (*ak4jets)[ik].phi();
                        ak4jet_e[ik] =   corr*uncorrJet.energy();
                        ak4jet_csv[ik] = (*ak4jets)[ik].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
                        ak4jet_icsv[ik] = (*ak4jets)[ik].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");   }
}
	sort (jets.begin (), jets.end (), ZmysortPt);
	for (size_t i=0;i<jets.size();i++) {
		if(iphoton>-1) {
			double drtmp1=deltaR(jets.at(i)->Eta(), jets.at(i)->Phi(), photoneta,photonphi);
			if(drtmp1>0.5 && jetindexphoton12[0]==-1&&jetindexphoton12[1]==-1) {
				jetindexphoton12[0] = i;
				continue;
			}
			if(drtmp1>0.5 && jetindexphoton12[0]!=-1&&jetindexphoton12[1]==-1) {
				jetindexphoton12[1] = i;
				continue;
			}
		}
	}

	for (size_t i=0;i<jets.size();i++) {
		if(iphoton_f>-1) {    
			double drtmp1_f=deltaR(jets.at(i)->Eta(), jets.at(i)->Phi(), photoneta_f,photonphi_f);
			if(drtmp1_f>0.5 && jetindexphoton12_f[0]==-1&&jetindexphoton12_f[1]==-1) {
				jetindexphoton12_f[0] = i;
				continue;  
			}
			if(drtmp1_f>0.5 && jetindexphoton12_f[0]!=-1&&jetindexphoton12_f[1]==-1) {
				jetindexphoton12_f[1] = i;
				continue; 
			}
		}		
	}

	if(jetindexphoton12[0]>-1 && jetindexphoton12[1]>-1) {
		jet1pt=jets[jetindexphoton12[0]]->Pt();
		jet1eta=jets[jetindexphoton12[0]]->Eta();
		jet1phi=jets[jetindexphoton12[0]]->Phi();
		jet1e=jets[jetindexphoton12[0]]->E();
		jet2pt=jets[jetindexphoton12[1]]->Pt();
		jet2eta=jets[jetindexphoton12[1]]->Eta();
		jet2phi=jets[jetindexphoton12[1]]->Phi();
		jet2e=jets[jetindexphoton12[1]]->E();
		jet1csv =(*ak4jets)[jetindexphoton12[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
		jet2csv =(*ak4jets)[jetindexphoton12[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
		jet1icsv =(*ak4jets)[jetindexphoton12[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		jet2icsv =(*ak4jets)[jetindexphoton12[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		drj1a=deltaR(jet1eta,jet1phi,photoneta,photonphi);
		drj2a=deltaR(jet2eta,jet2phi,photoneta,photonphi);
		drj1l=deltaR(jet1eta,jet1phi,etalep1,philep1);
		drj2l=deltaR(jet2eta,jet2phi,etalep1,philep1);
		drj1l2=deltaR(jet1eta,jet1phi,etalep2,philep2);
		drj2l2=deltaR(jet2eta,jet2phi,etalep2,philep2);
		TLorentzVector j1p4;
		j1p4.SetPtEtaPhiE(jet1pt, jet1eta, jet1phi, jet1e);
		TLorentzVector j2p4;
		j2p4.SetPtEtaPhiE(jet2pt, jet2eta, jet2phi, jet2e);
		TLorentzVector photonp42;
		photonp42.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
		TLorentzVector vp4;
		vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
		j1metPhi=fabs(jet1phi-MET_phi);
		if(j1metPhi>Pi) {j1metPhi=2.0*Pi-j1metPhi;}
		j2metPhi=fabs(jet2phi-MET_phi);
		if(j2metPhi>Pi) {j2metPhi=2.0*Pi-j2metPhi;}
		Mjj=(j1p4 + j2p4).M();
		deltaetajj = fabs(jet1eta - jet2eta);
		zepp = fabs((vp4+photonp42).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity())/ 2.0);
	}


	if(jetindexphoton12_f[0]>-1 && jetindexphoton12_f[1]>-1) {
		jet1pt_f=jets[jetindexphoton12_f[0]]->Pt();
		jet1eta_f=jets[jetindexphoton12_f[0]]->Eta();
		jet1phi_f=jets[jetindexphoton12_f[0]]->Phi();
		jet1e_f=jets[jetindexphoton12_f[0]]->E();
		jet2pt_f=jets[jetindexphoton12_f[1]]->Pt();
		jet2eta_f=jets[jetindexphoton12_f[1]]->Eta();
		jet2phi_f=jets[jetindexphoton12_f[1]]->Phi();
		jet2e_f=jets[jetindexphoton12_f[1]]->E();
		jet1csv_f =(*ak4jets)[jetindexphoton12_f[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
		jet2csv_f =(*ak4jets)[jetindexphoton12_f[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
		jet1icsv_f =(*ak4jets)[jetindexphoton12_f[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		jet2icsv_f =(*ak4jets)[jetindexphoton12_f[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		drj1a_f=deltaR(jet1eta_f,jet1phi_f,photoneta_f,photonphi_f);
		drj2a_f=deltaR(jet2eta_f,jet2phi_f,photoneta_f,photonphi_f);
		drj1l_f=deltaR(jet1eta_f,jet1phi_f,etalep1,philep1);
		drj2l_f=deltaR(jet2eta_f,jet2phi_f,etalep1,philep1);
		drj1l2_f=deltaR(jet1eta_f,jet1phi_f,etalep2,philep2);
		drj2l2_f=deltaR(jet2eta_f,jet2phi_f,etalep2,philep2);
		TLorentzVector j1p4_f;
		j1p4_f.SetPtEtaPhiE(jet1pt_f, jet1eta_f, jet1phi_f, jet1e_f);
		TLorentzVector j2p4_f;
		j2p4_f.SetPtEtaPhiE(jet2pt_f, jet2eta_f, jet2phi_f, jet2e_f);
		TLorentzVector photonp42_f;
		photonp42_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
		TLorentzVector vp4_f;
		vp4_f.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
		j1metPhi_f=fabs(jet1phi_f-MET_phi);
		if(j1metPhi_f>Pi) {j1metPhi_f=2.0*Pi-j1metPhi_f;}
		j2metPhi_f=fabs(jet2phi_f-MET_phi);
		if(j2metPhi_f>Pi) {j2metPhi_f=2.0*Pi-j2metPhi_f;}
		Mjj_f=(j1p4_f + j2p4_f).M();
		deltaetajj_f = fabs(jet1eta_f - jet2eta_f);
		zepp_f = fabs((vp4_f+photonp42_f).Rapidity() - (j1p4_f.Rapidity() + j2p4_f.Rapidity())/ 2.0);

	}

	outTree_->Fill();
//	std::cout<<"fill the outTree"<<std::endl;
	delete jecAK4_;
	jecAK4_=0;
}

//-------------------------------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------------------------------------------//

double ZPKUTreeMaker::getDR(double eta1, double phi1, double eta2, double phi2)
{
	double DR = -1e1;
	if (fabs(phi1-phi2)>3.1415926) DR=sqrt((eta1-eta2)*(eta1-eta2)+(2*3.1415926-fabs(phi1-phi2))*(2*3.1415926-fabs(phi1-phi2)));
	else DR=sqrt((eta1-eta2)*(eta1-eta2)+(fabs(phi1-phi2))*(fabs(phi1-phi2)));
	return DR;
}


void ZPKUTreeMaker::setDummyValues() {
	// muon station2 retrieve, L1 issue, Meng 2017/3/26
//	std::cout << "begin setDummyValues()..." << std::endl;
	lep1_sign = -1e2;
	lep2_sign = -1e2;
	lep1_etaphi_.first = -99.;
	lep1_etaphi_.second = -99.;
	lep2_etaphi_.first= -99.;
	lep2_etaphi_.second = -99.;
	lep1_eta_station2 = -99.;
	lep1_phi_station2 = -99.;
	lep2_eta_station2 = -99.;
	lep2_phi_station2 = -99.;
	//Lu
	npT=-1.;
	npIT=-1.;
	nBX=-1;
	nVtx           = -1e1;
	triggerWeight  = -1e1;
	pileupWeight   = -1e1;
	lumiWeight     = -1e1;
	theWeight = -99;
	lep            = -1e1;
	nlooseeles=-1e1;
	nloosemus=-1e1;
	ngoodmus = -1e1;
	ptVlep         = -1e1;
	yVlep          = -1e1;
	phiVlep        = -1e1;
	massVlep       = -1e1;
	Mla=-1e1;		Mla_f=-1e1;
	Mla2=-1e1;         Mla2_f=-1e1;
	Mva=-1e1;		Mva_f=-1e1;     
	ptlep1         = -1e1;
	etalep1        = -1e1;
	philep1        = -1e1;
	ptlep2         = -1e1;
	etalep2        = -1e1;
	philep2        = -1e1;
	// for muon rochester correction
	muon1_trackerLayers =-1e1;
	muon2_trackerLayers =-1e1;
	matchedgenMu1_pt =-1e2;
	matchedgenMu2_pt =-1e2;
	// for muon rochester correction
	met            = -1e1;
	metPhi         = -1e1;
	j1metPhi         = -1e1; j1metPhi_f         = -1e1;
	j2metPhi         = -1e1; j2metPhi_f         = -1e1;
	METraw_et = -99;
	METraw_phi = -99;
	METraw_sumEt = -99;
	genMET=-99;
	MET_et = -99;
	MET_phi = -99;
	MET_sumEt = -99;
//	for(int j=0; j<703; j++){
//		pweight[j]=0.0;
//	}
	for(int i=0; i<6;i++){
		genjet_pt[i]  = -1e1;
		genjet_eta[i]  = -1e1;
		genjet_phi[i]  = -1e1;
		genjet_e[i]  = -1e1;
	}      
	for(int i=0; i<6; i++) {
		genphoton_pt[i] = -1e1;
		genphoton_eta[i] = -1e1;
		genphoton_phi[i] = -1e1;
		genmuon_pt[i] = -1e1;
		genmuon_eta[i] = -1e1;
		genmuon_phi[i] = -1e1;
		genmuon_pid[i] = -1e2;
		genelectron_pt[i] = -1e1;
		genelectron_eta[i] = -1e1;
		genelectron_phi[i] = -1e1;
		photon_pt[i] = -1e1;
		photon_eta[i] = -1e1;
		photon_phi[i] = -1e1;
		photon_e[i] = -1e1;
		photon_pev[i] = false;
		photon_pevnew[i] = false;
		photon_ppsv[i] = false;
		photon_iseb[i] = false;
		photon_isee[i] = false;
		photon_hoe[i] = -1e1;
		photon_sieie[i] = -1e1;
		photon_sieie2[i] = -1e1;
		photon_chiso[i] = -1e1;
		photon_nhiso[i] = -1e1;
		photon_phoiso[i] = -1e1;
		photon_istrue[i] = -1;
		photon_isprompt[i] = -1;
		photon_drla[i] = 1e1;
		photon_drla2[i] = 1e1;
		photon_mla[i]=-1e1;
		photon_mla2[i]=-1e1;
		photon_mva[i]=-1e1;
		ak4jet_pt[i]=-1e1;
		ak4jet_eta[i]=-1e1;
		ak4jet_phi[i]=-1e1;
		ak4jet_e[i]=-1e1;
		ak4jet_csv[i]=-1e1;
		ak4jet_icsv[i]=-1e1;
	}

	photonet=-1e1;	 photonet_f=-1e1;
	photoneta=-1e1;  photoneta_f=-1e1;
	photonphi=-1e1;  photonphi_f=-1e1;
	photone=-1e1;   photone_f=-1e1;
	photonsieie=-1e1;  photonsieie_f=-1e1;
	photonphoiso=-1e1;  photonphoiso_f=-1e1;
	photonchiso=-1e1;  photonchiso_f=-1e1;
	photonnhiso=-1e1;  photonnhiso_f=-1e1;
	iphoton=-1;	 iphoton_f=-1;
	drla=1e1; 		 drla_f=1e1;
	drla2=1e1; 	 drla2_f=1e1;
	passEleVeto=false;
	passEleVetonew=false;
	passPixelSeedVeto=false;


	ISRPho = false;
	dR_ = 999;
	isTrue_=-1;
	isprompt_=-1; 
	jet1pt=-1e1;  jet1pt_f=-1e1;
	jet1eta=-1e1;  jet1eta_f=-1e1;
	jet1phi=-1e1;  jet1phi_f=-1e1;
	jet1e=-1e1;  jet1e_f=-1e1;
	jet1csv=-1e1;  jet1csv_f=-1e1;
	jet1icsv=-1e1;  jet1icsv_f=-1e1;
	jet2pt=-1e1;  jet2pt_f=-1e1;
	jet2eta=-1e1;  jet2eta_f=-1e1;
	jet2phi=-1e1;  jet2phi_f=-1e1;
	jet2e=-1e1;  jet2e_f=-1e1;
	jet2csv=-1e1;  jet2csv_f=-1e1;
	jet2icsv=-1e1;  jet2icsv_f=-1e1;
	drj1a=1e1;  drj1a_f=1e1;
	drj2a=1e1;  drj2a_f=1e1;
	drj1l=1e1;  drj1l_f=1e1;
	drj2l=1e1;  drj2l_f=1e1;
	drj1l2=1e1;  drj1l2_f=1e1;
	drj2l2=1e1;  drj2l2_f=1e1;
	Mjj=-1e1;   Mjj_f=-1e1;
	deltaetajj=-1e1;   deltaetajj_f=-1e1;
	zepp=-1e1;   zepp_f=-1e1;
	
	//L1 prefiring
	prefiringweight=-1.;
	prefiringweightup=-1.;
	prefiringweightdown=-1.;

	HLT_Ele1=-99;
	HLT_Ele2=-99;
	HLT_Ele3=-99;
	HLT_Ele4=-99;
	HLT_Ele5=-99;
	HLT_Ele6=-99;
	HLT_Ele7=-99;
	HLT_Mu1=-99;
	HLT_Mu2=-99;
	HLT_Mu3=-99;
	HLT_Mu4=-99;
	HLT_Mu5=-99;
	HLT_Mu6=-99;
	HLT_Mu7=-99;
	HLT_Mu8=-99;
	HLT_Mu9=-99;
        HLT_Mu10=-99;
        HLT_Mu11=-99;


	passFilter_HBHE_                  = false;
	passFilter_HBHEIso_               = false;
	passFilter_globalTightHalo_       = false;
	passFilter_ECALDeadCell_          = false;
	passFilter_GoodVtx_               = false;
	passFilter_EEBadSc_               = false;
	passFilter_badMuon_               = false;
	passFilter_badChargedHadron_      = false; 
	// Meng
	passFilter_MetbadMuon_	       = false;
	passFilter_duplicateMuon_ 	       = false;
//	std::cout << "end setDummyValues()..." << std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
	void 
ZPKUTreeMaker::beginJob()
{
//	std::cout << "ZPKUTreeMaker beginJob()..." << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void ZPKUTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
//	std::cout << "ZPKUTreeMaker beginRun()..." << std::endl;

	elPaths1.clear();
	elPaths2.clear();
	elPaths3.clear();
	elPaths4.clear();
	elPaths5.clear();
	elPaths6.clear();
	elPaths7.clear();
	muPaths1.clear();
	muPaths2.clear();
	muPaths3.clear();
	muPaths4.clear();
	muPaths5.clear();
	muPaths6.clear();
	muPaths7.clear();
	muPaths8.clear();
	muPaths9.clear();
	muPaths10.clear();
	muPaths11.clear();
	bool changed;
	if ( !hltConfig.init(iRun, iSetup, "HLT", changed) ) {
		edm::LogError("HltAnalysis") << "Initialization of HLTConfigProvider failed!!";
		return;
	}

	for (size_t i = 0; i < elPaths1_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), elPaths1_[i] );
		while ( !foundPaths.empty() ){
			elPaths1.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < elPaths2_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), elPaths2_[i] );
		while ( !foundPaths.empty() ){
			elPaths2.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < elPaths3_.size(); i++) {
                std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), elPaths3_[i] );
                while ( !foundPaths.empty() ){
                        elPaths3.push_back( foundPaths.back() );
                        foundPaths.pop_back();
                }
        }
	for (size_t i = 0; i < elPaths4_.size(); i++) {
                std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), elPaths4_[i] );
                while ( !foundPaths.empty() ){
                        elPaths4.push_back( foundPaths.back() );
                        foundPaths.pop_back();
                }
        }
	for (size_t i = 0; i < elPaths5_.size(); i++) {
                std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), elPaths5_[i] );
                while ( !foundPaths.empty() ){
                        elPaths5.push_back( foundPaths.back() );
                        foundPaths.pop_back();
                }
        }
	for (size_t i = 0; i < elPaths6_.size(); i++) {
                std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), elPaths6_[i] );
                while ( !foundPaths.empty() ){
                        elPaths6.push_back( foundPaths.back() );
                        foundPaths.pop_back();
                }
        }
	for (size_t i = 0; i < elPaths7_.size(); i++) {
                std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), elPaths7_[i] );
                while ( !foundPaths.empty() ){
                        elPaths7.push_back( foundPaths.back() );
                        foundPaths.pop_back();
                }
        }
	for (size_t i = 0; i < muPaths1_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths1_[i] );
		while ( !foundPaths.empty() ){
			muPaths1.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < muPaths2_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths2_[i] );
		while ( !foundPaths.empty() ){
			muPaths2.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < muPaths3_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths3_[i] );
		while ( !foundPaths.empty() ){
			muPaths3.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < muPaths4_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths4_[i] );
		while ( !foundPaths.empty() ){
			muPaths4.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < muPaths5_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths5_[i] );
		while ( !foundPaths.empty() ){
			muPaths5.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < muPaths6_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths6_[i] );
		while ( !foundPaths.empty() ){
			muPaths6.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < muPaths7_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths7_[i] );
		while ( !foundPaths.empty() ){
			muPaths7.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	for (size_t i = 0; i < muPaths8_.size(); i++) {
		std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths8_[i] );
		while ( !foundPaths.empty() ){
			muPaths8.push_back( foundPaths.back() );
			foundPaths.pop_back();
		}
	}
	std::cout<<"\n************** HLT Information **************\n";
	for (size_t i=0; i < elPaths1.size(); i++) std::cout << "\n Electron paths: " << elPaths1[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < elPaths2.size(); i++) std::cout << "\n Electron paths: " << elPaths2[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < muPaths1.size(); i++) std::cout << "\n Muon paths1    : " << muPaths1[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < muPaths2.size(); i++) std::cout << "\n Muon paths2    : " << muPaths2[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < muPaths3.size(); i++) std::cout << "\n Muon paths3    : " << muPaths3[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < muPaths4.size(); i++) std::cout << "\n Muon paths4    : " << muPaths4[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < muPaths5.size(); i++) std::cout << "\n Muon paths5    : " << muPaths5[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < muPaths6.size(); i++) std::cout << "\n Muon paths6    : " << muPaths6[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < muPaths7.size(); i++) std::cout << "\n Muon paths7    : " << muPaths7[i].c_str() <<"\t"<< std::endl;
	for (size_t i=0; i < muPaths8.size(); i++) std::cout << "\n Muon paths8    : " << muPaths8[i].c_str() <<"\t"<< std::endl;
	std::cout<<"\n*********************************************\n\n";

}

void ZPKUTreeMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
	std::cout << "ZPKUTreeMaker endRun()..." << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZPKUTreeMaker::endJob() {
	std::cout << "ZPKUTreeMaker endJob()..." << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZPKUTreeMaker);
