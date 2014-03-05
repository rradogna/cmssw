#include <Geometry/GEMGeometry/interface/GEMGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include "FWCore/Framework/interface/ESHandle.h"

#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHitCollection.h>
#include <RecoLocalMuon/GEMRecHit/interface/CSCSegtoGEM.h>

ObjectMapCSC* ObjectMapCSC::mapInstance = NULL;

ObjectMapCSC* ObjectMapCSC::GetInstance(const edm::EventSetup& iSetup){
  if (mapInstance == NULL){
    mapInstance = new ObjectMapCSC(iSetup);
  }
  return mapInstance;
}

ObjectMapCSC::ObjectMapCSC(const edm::EventSetup& iSetup){
  edm::ESHandle<GEMGeometry> gemGeo;
  edm::ESHandle<CSCGeometry> cscGeo;
  
  iSetup.get<MuonGeometryRecord>().get(gemGeo);
  iSetup.get<MuonGeometryRecord>().get(cscGeo);
  
  for (TrackingGeometry::DetContainer::const_iterator it=gemGeo->dets().begin();it<gemGeo->dets().end();it++){
    if(dynamic_cast< GEMChamber* >( *it ) != 0 ){
      GEMChamber* ch = dynamic_cast< GEMChamber* >( *it ); 
      std::vector< const GEMEtaPartition*> roles = (ch->etaPartitions());
      for(std::vector<const GEMEtaPartition*>::const_iterator r = roles.begin();r != roles.end(); ++r){
	GEMDetId gemId = (*r)->id();
	int region=gemId.region();
	if(region!=0){
	  int station=gemId.station();
          int ring=gemId.ring();
	  int gemchamber=gemId.chamber();
          int cscring=ring;
          int cscstation=station;
	  int cscchamber = gemchamber;
          
	  CSCStationIndex ind(region,cscstation,cscring,cscchamber);
          std::set<GEMDetId> myrolls;
	  if (rollstoreCSC.find(ind)!=rollstoreCSC.end()) myrolls=rollstoreCSC[ind];
	  myrolls.insert(gemId);
          rollstoreCSC[ind]=myrolls;
	}
      }
    }
  }
}
  
CSCSegtoGEM::CSCSegtoGEM(edm::Handle<CSCSegmentCollection> allCSCSegments, const edm::EventSetup& iSetup,const edm::Event& iEvent, bool debug, double eyr){
  
  edm::ESHandle<GEMGeometry> gemGeo;
  edm::ESHandle<CSCGeometry> cscGeo;
  
  iSetup.get<MuonGeometryRecord>().get(gemGeo);
  iSetup.get<MuonGeometryRecord>().get(cscGeo);
  

  if(debug) std::cout<<"CSC \t Number of CSC Segments in this event = "<<allCSCSegments->size()<<std::endl;

  _ThePoints = new GEMRecHitCollection();

  if(allCSCSegments->size()==0){
    if(debug) std::cout<<"CSC 0 segments skiping event"<<std::endl;
  }else {
    std::map<CSCDetId,int> CSCSegmentsCounter;
    CSCSegmentCollection::const_iterator segment;
      
    int segmentsInThisEventInTheEndcap=0;
      
    for (segment = allCSCSegments->begin();segment!=allCSCSegments->end(); ++segment){
      CSCSegmentsCounter[segment->cscDetId()]++;
      segmentsInThisEventInTheEndcap++;
    }    
      
    if(debug) std::cout<<"CSC \t loop over all the CSCSegments "<<std::endl;
    for (segment = allCSCSegments->begin();segment!=allCSCSegments->end(); ++segment){
      CSCDetId CSCId = segment->cscDetId();
	
      if(debug) std::cout<<"CSC \t \t This Segment is in Chamber id: "<<CSCId<<std::endl;
      if(debug) std::cout<<"CSC \t \t Number of segments in this CSC = "<<CSCSegmentsCounter[CSCId]<<std::endl;
      //if(debug) std::cout<<"CSC \t \t Is the only one in this CSC? is not ind the ring 1 or station 4? Are there more than 2 segments in the event?"<<std::endl;
       if(debug) std::cout<<"CSC segments diff to zero?"<<std::endl;

      //if(CSCSegmentsCounter[CSCId]==1 && CSCId.station()!=4 && CSCId.ring()!=1 && allCSCSegments->size()>=2){
      if(allCSCSegments->size()!=0){ //???????? inutile, ma che selezione dovrei mettere?
	if(debug) std::cout<<"CSC \t \t yes"<<std::endl;
	int cscEndCap = CSCId.endcap();
	int cscStation = CSCId.station();
	int cscRing = CSCId.ring();
	int gemRegion = 1; if(cscEndCap==2) gemRegion= -1;//Relacion entre las endcaps
	int gemRing = cscRing;
	if(cscRing==4)gemRing =1; ////csc ring=4 è me1A, mentre ring 1 è me1B, gemring=1
	int gemStation = cscStation;
	int gemChamber = CSCId.chamber();
	
	//LocalPoint segmentPosition= segment->localPosition();
	LocalVector segmentDirection=segment->localDirection();
	float dz=segmentDirection.z();
	
	if(debug) std::cout<<"CSC \t \t \t Information about the segment" 
			   <<"RecHits ="<<segment->nRecHits()
			   <<"Angle ="<<acos(dz)*180/3.1415926<<std::endl;
		      
	if(debug) std::cout<<"CSC \t \t Is a good Segment? dim = 4, 4 <= nRecHits <= 10 Incident angle int range 45 < "<<acos(dz)*180/3.1415926<<" < 135? "<<std::endl;

	if((segment->dimension()==4) && (segment->nRecHits()<=10 && segment->nRecHits()>=4)){ ///////??????????????
	  //&& acos(dz)*180/3.1415926 > 45. && acos(dz)*180/3.1415926 < 135.){ 
	  //&& segment->chi2()< ??)Add 3 segmentes in the endcaps???


	  if(debug) std::cout<<"CSC \t \t yes"<<std::endl;
	  if(debug) std::cout<<"CSC \t \t CSC Segment Dimension "<<segment->dimension()<<std::endl; 
	    
	  /*float Xo=segmentPosition.x();
	  float Yo=segmentPosition.y();
	  float Zo=segmentPosition.z();
	  float dx=segmentDirection.x();
	  float dy=segmentDirection.y();
	  float dz=segmentDirection.z();
	*/
	  if(debug)  std::cout<<"Calling to Object Map class"<<std::endl;
	  ObjectMapCSC* TheObjectCSC = ObjectMapCSC::GetInstance(iSetup);
	  if(debug) std::cout<<"Creating the CSCIndex"<<std::endl;
	  CSCStationIndex theindex(gemRegion,gemStation,gemRing,gemChamber);
	  if(debug) std::cout<<"Getting the Rolls for the given index"<<std::endl;
	
	  std::set<GEMDetId> rollsForThisCSC = TheObjectCSC->GetInstance(iSetup)->GetRolls(theindex);
		
	   
	  if(debug) std::cout<<"CSC \t \t Getting chamber from Geometry"<<std::endl;
	  const CSCChamber* TheChamber=cscGeo->chamber(CSCId); 
	  if(debug) std::cout<<"CSC \t \t Getting ID from Chamber"<<std::endl;
	  const CSCDetId TheId=TheChamber->id();

	  if(debug) std::cout<<"CSC \t \t Number of rolls for this CSC = "<<rollsForThisCSC.size()<<std::endl;

	  if(debug) std::cout<<"CSC \t \t Printing The Id"<<TheId<<std::endl;
	  
//???????????????
	/*  if(rpcRing!=1&&rpcStation!=4){//They don't exist! //////??????
	  
	    assert(rollsForThisCSC.size()>=1);

	    if(debug) std::cout<<"CSC \t \t Loop over all the rolls asociated to this CSC"<<std::endl;	    
	    for (std::set<GEMDetId>::iterator iteraRoll = rollsForThisCSC.begin();iteraRoll != rollsForThisCSC.end(); iteraRoll++){
	      const GEMEtaPartition* rollasociated = gemGeo->etaPartitions(*iteraRoll);
	      GEMDetId gemId = rollasociated->id();
		
	      if(debug) std::cout<<"CSC \t \t \t We are in the roll getting the surface"<<gemId<<std::endl;
	      const BoundPlane & GEMSurface = rollasociated->surface(); 

	      if(debug) std::cout<<"CSC \t \t \t RollID: "<<gemId<<std::endl;
		
	      if(debug) std::cout<<"CSC \t \t \t Doing the extrapolation to this roll"<<std::endl;
	      if(debug) std::cout<<"CSC \t \t \t CSC Segment Direction in CSCLocal "<<segmentDirection<<std::endl;
	      if(debug) std::cout<<"CSC \t \t \t CSC Segment Point in CSCLocal "<<segmentPosition<<std::endl;  
		
	      GlobalPoint CenterPointRollGlobal = GEMSurface.toGlobal(LocalPoint(0,0,0));
	      if(debug) std::cout<<"CSC \t \t \t Center (0,0,0) of the Roll in Global"<<CenterPointRollGlobal<<std::endl;
	      GlobalPoint CenterPointCSCGlobal = TheChamber->toGlobal(LocalPoint(0,0,0));
	      if(debug) std::cout<<"CSC \t \t \t Center (0,0,0) of the CSC in Global"<<CenterPointCSCGlobal<<std::endl;
	      GlobalPoint segmentPositionInGlobal=TheChamber->toGlobal(segmentPosition); //new way to convert to global
	      if(debug) std::cout<<"CSC \t \t \t Segment Position in Global"<<segmentPositionInGlobal<<std::endl;
	      LocalPoint CenterRollinCSCFrame = TheChamber->toLocal(CenterPointRollGlobal);
*/
	      /*if(debug){//to check CSC RPC phi relation!
		float rpcphi=0;
		float cscphi=0;
		  
		(CenterPointRollGlobal.barePhi()<0)? 
		  rpcphi = 2*3.141592+CenterPointRollGlobal.barePhi():rpcphi=CenterPointRollGlobal.barePhi();
		  
		(CenterPointCSCGlobal.barePhi()<0)? 
		  cscphi = 2*3.1415926536+CenterPointCSCGlobal.barePhi():cscphi=CenterPointCSCGlobal.barePhi();

		float df=fabs(cscphi-rpcphi); 
		float dr=fabs(CenterPointRollGlobal.perp()-CenterPointCSCGlobal.perp());
		float diffz=CenterPointRollGlobal.z()-CenterPointCSCGlobal.z();
		float dfg=df*180./3.14159265;

		if(debug) std::cout<<"CSC \t \t \t z of RPC="<<CenterPointRollGlobal.z()<<"z of CSC"<<CenterPointCSCGlobal.z()<<" dfg="<<dfg<<std::endl;
		  
		RPCGeomServ rpcsrv(rpcId);
		  
		if(dr>200.||fabs(dz)>55.||dfg>1.){ 
		  //if(rpcRegion==1&&dfg>1.&&dr>100.){  
		  if (debug) std::cout
		    <<"\t \t \t CSC Station= "<<CSCId.station()
		    <<" Ring= "<<CSCId.ring()
		    <<" Chamber= "<<CSCId.chamber()
		    <<" cscphi="<<cscphi*180/3.14159265
		    <<"\t RPC Station= "<<rpcId.station()
		    <<" ring= "<<rpcId.ring()
		    <<" segment =-> "<<rpcsrv.segment()
		    <<" rollphi="<<rpcphi*180/3.14159265
		    <<"\t dfg="<<dfg
		    <<" dz="<<diffz
		    <<" dr="<<dr
		    <<std::endl;
		    
		}
	      }// end if debug*/
	      
	    

	    //}
	  //}
	}
      }
    }
  }
}

CSCSegtoGEM::~CSCSegtoGEM(){

}

