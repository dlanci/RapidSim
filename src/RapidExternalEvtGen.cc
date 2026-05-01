#include "RapidExternalEvtGen.h"

#include <cstdlib>
#include <fstream>
#include <queue>

#include "TRandom.h"
#include "TSystem.h"

#ifdef RAPID_EVTGEN
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"

#include "EvtGenExternal/EvtExternalGenList.hh"
#endif

int B0Id( 511 ), B0barId( -511 );
int BsId( 531 ), BsbarId( -531 );


bool RapidExternalEvtGen::decay(std::vector<RapidParticle*>& parts) {
#ifdef RAPID_EVTGEN

	if(parts.size() < 1) {
		std::cout << "WARNING in RapidExternalEvtGen::decay : There are no particles to decay." << std::endl;
		return false;
	}

	// EvtParticle* baseParticle(0);
	EvtSpinDensity* spinDensity = 0;
	
	TString parentName = getEvtGenName(parts[0]->id());
	EvtId theId = EvtPDL::evtIdFromLundKC(parts[0]->id());

	int PDGId = EvtPDL::getStdHep( theId );

    if ( theId.getId() == -1 && theId.getAlias() == -1 ) {
        std::cout << "Error. Could not find valid EvtId for " << parentName << std::endl;
        return -1;
    }

	EvtSpinType::spintype baseSpin = EvtPDL::getSpinType( theId );

    if ( baseSpin == EvtSpinType::VECTOR ) {
        std::cout << "Setting spin density for vector particle " << parentName << std::endl;
        spinDensity = new EvtSpinDensity();
        spinDensity->setDiag( EvtSpinType::getSpinStates( EvtSpinType::VECTOR ) );
        spinDensity->set( 1, 1, EvtComplex( 0.0, 0.0 ) );
    }

	// RapidVertex* vtx(0);
	// vtx = parts[0]->getOriginVertex();
	
	// Parent particle XYZ-position
	// ROOT::Math::XYZPoint point = vtx->getVertex(true);
	// EvtVector4R origin(0.0, point.X(), point.Y(), point.Z());
	EvtVector4R origin(0.0, 0.0, 0.0, 0.0);

	// Parent particle 4-momentum
	TLorentzVector pIn = parts[0]->getP();
	EvtVector4R pInit(pIn.E(),pIn.Px(),pIn.Py(),pIn.Pz());


    EvtHepMCEvent* theEvent =
        evtGen_->generateDecay( PDGId, pInit, origin, spinDensity );

	// Retrieve the HepMC event information
    GenEvent* hepMCEvent = theEvent->getEvent();
    // hepMCEvent->print();

	std::list<GenVertexPtr> allVertices;
	
	double flightTime;

	int iVtx(0);
	int iPart(1);
	bool hasOsc(false);

	TLorentzVector p4TLV;
	FourVector FourMom;
	FourVector DecayVtx;
	FourVector OrigVtx;

	for ( auto theVertex : hepMCEvent->vertices() ) {
		if ( theVertex == 0 ) {
		            continue;
		}
		auto nin  = theVertex->particles_in_size();
		auto nout = theVertex->particles_out_size();

		
		if(nin==1 && nout ==1){

			auto inParticle  = theVertex->particles_in()[0];
			auto outParticle = theVertex->particles_out()[0];
			int inPDGId = inParticle->pdg_id();
			int PDGId = outParticle->pdg_id();
			// std::cout<<""<<std::endl;
			// std::cout<<"In nin==1 and nout==1"<<std::endl;
			// std::cout<<""<<std::endl;
			if ( inPDGId == B0Id || inPDGId == BsId || inPDGId == B0barId || inPDGId == BsbarId ) {


				// std::cout<<"inPDGId"<<inPDGId<<std::endl;
				// std::cout<<"PDGId"<<PDGId<<std::endl;
				
				if (PDGId!=inPDGId){
					hasOsc=true;
					// std::cout<<""<<std::endl;
					// std::cout<<" HAS OSCILLATED "<<std::endl;
					// std::cout<<" inPDGId: "<<inPDGId<<std::endl;
					// std::cout<<" PDGId: "<<PDGId<<std::endl;
					// std::cout<<""<<std::endl;
				}
				else{
					hasOsc=false;
					// std::cout<<""<<std::endl;
					// std::cout<<" HASN'T OSCILLATED "<<std::endl;
					// std::cout<<" inPDGId: "<<inPDGId<<std::endl;
					// std::cout<<" PDGId: "<<PDGId<<std::endl;
					// std::cout<<""<<std::endl;
				}
				
        	
			}

			continue;
			
		}
		else if(nin==1 && nout > 1){
			// For these, get the mother decay vertex position and the 4-momentum to calculate
			// the flight time.

			for ( auto inParticle : theVertex->particles_in() ) {
        	if ( inParticle == 0 ) {
        	    continue;
        	}
			
			FourMom = inParticle->momentum();
			DecayVtx = theVertex->position();
			
			flightTime = calcFlightTime( DecayVtx, FourMom );
			
			parts[iVtx]->setId(inParticle->pdg_id());
			parts[iVtx]->setDecaytime(flightTime);
			parts[iVtx]->getDecayVertex()->setXYZ(DecayVtx.x(), DecayVtx.y(), DecayVtx.x());
			
			int inPDGId = inParticle->pdg_id();
			// std::cout<<""<<std::endl;
			// std::cout<<"In nin==1 and nout>1"<<std::endl;
			// std::cout<<""<<std::endl;			
			if(iVtx==0){
        		if ( inPDGId == B0Id || inPDGId == BsId || inPDGId == B0barId || inPDGId == BsbarId ) {

					
					if (PDGId!=inPDGId){
						hasOsc=true;
						// std::cout<<""<<std::endl;
						// std::cout<<" HAS OSCILLATED "<<std::endl;
						// std::cout<<" inPDGId: "<<inPDGId<<std::endl;
						// std::cout<<" PDGId: "<<PDGId<<std::endl;
						// std::cout<<""<<std::endl;
					}
					else{
						hasOsc=false;
						// std::cout<<""<<std::endl;
						// std::cout<<" HASN'T OSCILLATED "<<std::endl;
						// std::cout<<" inPDGId: "<<inPDGId<<std::endl;
						// std::cout<<" PDGId: "<<PDGId<<std::endl;
						// std::cout<<""<<std::endl;
					}

					parts[iVtx]->setHasOsc(hasOsc);
        		}
			}
       		}  

			for ( auto outParticle : theVertex->particles_out() ) {
			
			if(outParticle->pdg_id() != 22){//discard PHOTOS info
				
				FourMom = outParticle->momentum();
				OrigVtx = theVertex->position();

				p4TLV.SetPxPyPzE(FourMom.px(),FourMom.py(),FourMom.pz(),FourMom.e());

				parts[iPart]->setP(p4TLV);

				parts[iPart]->getOriginVertex()->setXYZ(OrigVtx.x(),OrigVtx.y(),OrigVtx.z());
				parts[iPart]->setId(outParticle->pdg_id());

				++iPart;
				
			}
			}
			++iVtx;
			// std::cout<<"hey: "<<iVtx<<std::endl;
		}
	}

	delete hepMCEvent;

	return true;
#else
	if(!suppressWarning_) {
		std::cout << "WARNING in RapidExternalEvtGen::decay : EvtGen extension not compiled. Will not use EvtGen to decay " << parts[0]->name() << "." << std::endl;
		suppressWarning_=true;
	}

	return false;
#endif
}

// bool RapidExternalEvtGen::decay(std::vector<RapidParticle*>& parts) {
// #ifdef RAPID_EVTGEN

// 	if(parts.size() < 1) {
// 		std::cout << "WARNING in RapidExternalEvtGen::decay : There are no particles to decay." << std::endl;
// 		return false;
// 	}

// 	EvtParticle* baseParticle(0);

// 	EvtId theId = EvtPDL::evtIdFromLundKC(parts[0]->id());
	

// 	// Parent particle 4-momentum
// 	TLorentzVector pIn = parts[0]->getP();
// 	EvtVector4R pInit(pIn.E(),pIn.Px(),pIn.Py(),pIn.Pz());

// 	baseParticle = EvtParticleFactory::particleFactory(theId, pInit);
// 	if (baseParticle->getSpinStates() == 3) {baseParticle->setVectorSpinDensity();}

// 	// std::cout<<""<<std::endl;
// 	// std::cout<<""<<std::endl;
// 	// std::cout<<"before"<<std::endl;
// 	// baseParticle->printTree();
// 	int motherIDbefore(EvtPDL::getStdHep(baseParticle->getId()));
// 	// std::cout<<""<<std::endl;
// 	// Generate the event
// 	evtGen_->generateDecay(baseParticle);
// 	// std::cout<<""<<std::endl;
// 	// std::cout<<"after"<<std::endl;
// 	// baseParticle->printTree();
// 	int motherIDafter(EvtPDL::getStdHep(baseParticle->getId()));
// 	// std::cout<<""<<std::endl;
// 	// std::cout<<""<<std::endl;

// 	// Store particles to read in the order RapidSim stores them
// 	std::queue<EvtParticle*> evtParts;
// 	// Also store the number of children expected for each of these particles so we can remove PHOTOS photons
// 	std::queue<int> nExpectedChildren;
		
// 	evtParts.push(baseParticle);
// 	nExpectedChildren.push(parts[0]->nDaughters());

// 	EvtVector4R x4Evt;
// 	EvtVector4R p4Evt;
// 	TLorentzVector p4TLV;

// 	int iPart=1; // The momentum and origin vertex of the first particle are already set

// 	parts[0]->setId(EvtPDL::getStdHep(baseParticle->getId()));
	
// 	// bool isBMixed((isB0Mixed(baseParticle) || isBsMixed(baseParticle)));
// 	bool isBMixed(motherIDbefore!=motherIDafter);

// 	std::cout<<"isBMixed :"<<isBMixed<<std::endl;
// 	if (isBMixed){
// 		std::cout<<"motherIDbefore: "<<motherIDbefore <<std::endl;
// 		std::cout<<"motherIDafter:  "<<motherIDafter  <<std::endl;
// 	}

// 	//get the flight distance from EvtGen somehow
// 	//option 1) use (L*m)/(|p|*c)
// 	//option 2) use t = gamma*tau
// 	double gamma, tlab;

// 	p4Evt = baseParticle->getP4Lab();
// 	gamma = p4Evt.get(0)/baseParticle->mass();
// 	tlab = gamma*baseParticle->getLifetime();


// 	std::cout<<"Bmass: "<<baseParticle->mass()<<std::endl;
// 	std::cout<<"energy: "<<p4Evt.get(0)<<std::endl;
// 	std::cout<<"gamma: "<<gamma<<std::endl;
// 	std::cout<<"t"<<baseParticle->getLifetime()<<std::endl;
// 	std::cout<<"tlab"<<tlab<<std::endl;


// 	parts[0]->setDecaytime(tlab);

// 	while(!evtParts.empty()) {
		
// 		EvtParticle* theParticle = evtParts.front();

// 		// // B0 and Bs may mix in EvtGen - RapidSim ignores this step and only records the second state
// 		//if isBMixed?
// 		while(theParticle->getNDaug()==1) {
// 			theParticle = theParticle->getDaug(0);
// 		}		
	
// 		uint nChildren = nExpectedChildren.front();
// 		// Loop over the daughter tracks
// 		for (uint iChild = 0; iChild < nChildren; ++iChild) {

// 			EvtParticle* child = theParticle->getDaug(iChild);

// 			if (child != 0) {
// 				p4Evt = child->getP4Lab();
// 				x4Evt = child->get4Pos();
// 				p4TLV.SetPxPyPzE(p4Evt.get(1),p4Evt.get(2),p4Evt.get(3),p4Evt.get(0));
// 				if(parts.size() < iPart+1u) {
// 					std::cout << "WARNING in RapidExternalEvtGen::decay : EvtGen has produced too many particles." << std::endl;
// 					return false;
// 				}
// 				parts[iPart]->setP(p4TLV);
// 				parts[iPart]->getOriginVertex()->setXYZ(x4Evt.get(1),x4Evt.get(2),x4Evt.get(3));
// 				parts[iPart]->setId(EvtPDL::getStdHep(child->getId()));
// 				parts[iPart]->setDecaytime(child->getLifetime());
// 				evtParts.push(child);
// 				nExpectedChildren.push(parts[iPart]->nDaughters());
// 				++iPart;
// 			}
// 		}
// 		// Clean up any PHOTOS photons
// 		for (uint iChild = nChildren; iChild < theParticle->getNDaug(); ++iChild) {
// 			delete theParticle->getDaug(iChild);
// 		}
// 		delete theParticle;
// 		evtParts.pop();
// 		nExpectedChildren.pop();
// 	}
// 	return true;
// #else
// 	if(!suppressWarning_) {
// 		std::cout << "WARNING in RapidExternalEvtGen::decay : EvtGen extension not compiled. Will not use EvtGen to decay " << parts[0]->name() << "." << std::endl;
// 		suppressWarning_=true;
// 	}

// 	return false;
// #endif
// }

bool RapidExternalEvtGen::setup() {
#ifdef RAPID_EVTGEN
	std::cout << "INFO in RapidExternalEvtGen::setup : Setting decay for external EvtGen generator." << std::endl;
	if(!evtGen_) setupGenerator();
	evtGen_->readUDecay(decFileName_.Data());
	return true;
#else
	std::cout << "WARNING in RapidExternalEvtGen::setup : EvtGen extension not compiled." << std::endl;
	return false;
#endif
}

bool RapidExternalEvtGen::setupGenerator() {
#ifdef RAPID_EVTGEN
	std::cout << "INFO in RapidExternalEvtGen::setupGenerator : Setting up external EvtGen generator." << std::endl;
	EvtRandomEngine* randomEngine = 0;
	EvtAbsRadCorr* radCorrEngine = 0;
	std::list<EvtDecayBase*> extraModels;

	// Define the random number generator
	uint seed = gRandom->TRandom::GetSeed();
	randomEngine = new EvtMTRandomEngine(seed);

	bool useEvtGenRandom(false);
	EvtExternalGenList genList(true, "", "gamma", useEvtGenRandom);
	radCorrEngine = genList.getPhotosModel();
	extraModels = genList.getListOfModels();

	TString evtPDLPath;
	evtPDLPath += getenv("EVTGEN_ROOT");
	evtPDLPath += "/evt.pdl";

	bool foundDec=false;

	TString decPath;
	decPath += getenv("RAPIDSIM_CONFIG");
	if(decPath!="") {
		decPath += "/config/evtgen/DECAY.DEC";
		if(!gSystem->AccessPathName(decPath)) foundDec=true;
	}

	// We want to initialise EvtGen before we define our DEC file so we can use EvtPDL
	// To do this pass an empty DEC file as the main decay file and pass our file later as a user file
	if(!foundDec) {
		decPath += getenv("RAPIDSIM_ROOT");
		decPath += "/config/evtgen/DECAY.DEC";
	}

	int mixingType = EvtCPUtil::Incoherent;
	// int mixingType = 1;

	evtGen_ = new EvtGen(decPath.Data(), evtPDLPath.Data(), randomEngine,
			radCorrEngine, &extraModels, mixingType);

	return true;
#else
	std::cout << "WARNING in RapidExternalEvtGen::setup : EvtGen extension not compiled." << std::endl;
	return false;
#endif
}

void RapidExternalEvtGen::writeDecFile(TString fname, std::vector<RapidParticle*>& parts, bool usePhotos) {
#ifdef RAPID_EVTGEN
	if(!evtGen_) setupGenerator();

	decFileName_ = fname+".DEC";
	std::cout << "INFO in RapidExternalEvtGen::writeDecFile : Writing EvtGen DEC file : " << decFileName_ << std::endl;

	std::ofstream fout;
	fout.open(decFileName_, std::ofstream::out);

	if(usePhotos) {
		fout << "yesPhotos\n" << std::endl;
	} else {
		fout << "noPhotos\n" << std::endl;
	}

	// Loop over all particles and write out Decay rule for each
	for(unsigned int iPart=0; iPart<parts.size(); ++iPart) {
		unsigned int nChildren = parts[iPart]->nDaughters();
		if(nChildren>0) {
			int id = parts[iPart]->id();
			fout << "Decay " << getEvtGenName(id) << "\n1.00\t";
			if ( !(parts[iPart]->evtGenDecayModel()).Contains("TAUOLA") ) {
				for(unsigned int iChild=0; iChild<nChildren; ++iChild) {
					fout << getEvtGenName(parts[iPart]->daughter(iChild)->id()) << "\t";
				}
			}
			fout << parts[iPart]->evtGenDecayModel() << ";" << std::endl;
			fout <<"Enddecay" << std::endl;

			// Workaround to deal with mixing of B0 and Bs
			if(TMath::Abs(id)==531||TMath::Abs(id)==511) fout <<"CDecay " << getEvtGenConjName(id) << std::endl << std::endl;
		}
	}
	fout <<"End\n" << std::endl;
	fout.close();
#else
	std::cout << "WARNING in RapidExternalEvtGen::writeDecFile : EvtGen extension not compiled. Cannot write DEC file "
		  << fname << " for " << parts.size() << "particles with usePhotos=" << usePhotos << "." << std::endl;
#endif
}

TString RapidExternalEvtGen::getEvtGenName(int id) {
#ifdef RAPID_EVTGEN
	EvtId evtId = EvtPDL::evtIdFromStdHep(id);
	TString name = EvtPDL::name(evtId);
	return name;
#else
	std::cout << "WARNING in RapidExternalEvtGen::getEvtGenName : EvtGen extension not compiled. Cannot lookup name for particle ID " << id << "." << std::endl;
	return "";
#endif
}

TString RapidExternalEvtGen::getEvtGenConjName(int id) {
#ifdef RAPID_EVTGEN
	EvtId evtId = EvtPDL::evtIdFromStdHep(id);
	EvtId evtConjId = EvtPDL::chargeConj(evtId);
	TString name = EvtPDL::name(evtConjId);
	return name;
#else
	std::cout << "WARNING in RapidExternalEvtGen::getEvtGenConjName : EvtGen extension not compiled. Cannot lookup conjugate name for particle ID " << id << "." << std::endl;
	return "";
#endif
}

bool RapidExternalEvtGen::isB0Mixed( EvtParticle* p )
{
	if ( !( p->getDaug(0) ) )
        return false;

    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );
 
    if ( ( p->getId() != B0 ) && ( p->getId() != B0B ) ){
		
        return false;
	}

    if ( ( p->getDaug(0)->getId() == B0 ) || ( p->getDaug(0)->getId() == B0B ) ){
    	
		return true;
	}

    return false;
}
 
bool RapidExternalEvtGen::isBsMixed( EvtParticle* p )
{

	if ( !( p->getDaug(0) ) )
        return false;

    static EvtId BS0 = EvtPDL::getId( "B_s0" );
    static EvtId BSB = EvtPDL::getId( "anti-B_s0" );

    if ( ( p->getId() != BS0 ) && ( p->getId() != BSB ) )
	{
        return false;
	}
    if ( ( p->getDaug(0)->getId() == BS0 ) || ( p->getDaug(0)->getId() == BSB ) )
	{
        return true;
	}
    return false;
}

double RapidExternalEvtGen::calcFlightTime( FourVector& DecayVtx, FourVector& P4mtm )
{
    double flightTime( 0.0 );

#ifdef EVTGEN_HEPMC3
    double distance = DecayVtx.length() * 1e-3;    // metres
    double momentum = P4mtm.length();               // GeV/c
	double PMass = P4mtm.m();
#else
    double distance = DecayVtx.rho() * 1e-3;    // metres
    double momentum = P4mtm.rho();  
	double PMass = P4mtm.m();
#endif
	// std::cout<<"mass: "<<PMass<<std::endl;
	// std::cout<<"energy: "<<P4mtm.e()<<std::endl;
	// std::cout<<"distance: "<<distance<<std::endl;
	// std::cout<<"momentum: "<<momentum<<std::endl;
	// std::cout<<"gamma: "<<P4mtm.e()/PMass<<std::endl;

    
    double c0 = 299792458.0;    // m/s
	
    if ( momentum > 0.0 ) {
        flightTime = 1.0e12 * distance * PMass /
                     ( momentum * c0 );    // picoseconds
    }

    return flightTime;
}
