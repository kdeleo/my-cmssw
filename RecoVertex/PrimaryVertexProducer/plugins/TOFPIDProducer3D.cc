#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

using namespace std;
using namespace edm;

class TOFPIDProducer3D : public edm::stream::EDProducer<> {
public:
  TOFPIDProducer3D(const ParameterSet& pset);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  template <class H, class T>
  void fillValueMap(edm::Event& iEvent,
                    const edm::Handle<H>& handle,
                    const std::vector<T>& vec,
                    const std::string& name) const;

  void produce(edm::Event& ev, const edm::EventSetup& es) final;

private:
  static constexpr char t0Name[] = "t0";
  static constexpr char sigmat0Name[] = "sigmat0";
  static constexpr char t0safeName[] = "t0safe";
  static constexpr char sigmat0safeName[] = "sigmat0safe";
  static constexpr char probPiName[] = "probPi";
  static constexpr char probKName[] = "probK";
  static constexpr char probPName[] = "probP";

  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofkToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofpToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxsToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeErrorToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDMomentumToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDPathLengthToken_;

  struct TrackInfo {
    reco::TrackRef trackreference;
    double trkWeight;
    double trkTime;
    double trkTimeError;
    double trkTimeHyp[3];
  };

  edm::ValueMap<float> trackMTDTimes_;
  edm::ValueMap<float> trackMTDTimeErrors_;
  edm::ValueMap<float> trackMTDTimeQualities_;
  edm::ValueMap<float> trackMTDMomenta_;
  edm::ValueMap<float> trackMTDPathLengths_;

  double const minTrackVtxWeight_;
  double const minTrackTimeQuality_;
  double const massPion_;
  double const massKaon_;
  double const massProton_;
  double const probPion_;
  double const probKaon_;
  double const probProton_;
  double const Tstart_;
  double const coolingFactor_;

  double fixedT0Error_;

  float trackTime(float const mass, float const mtdTime, float const mtdPathLength, float const mtdMomentum) const;

};

TOFPIDProducer3D::TOFPIDProducer3D(const ParameterSet& iConfig)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"))),
      t0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"))),
      sigmat0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"))),
      tofkToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofkSrc"))),
      tofpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofpSrc"))),
      vtxsToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxsSrc"))),
      trackMTDTimeToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMTDTimeVMapTag"))),
      trackMTDTimeErrorToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMTDTimeErrorVMapTag"))),
      trackMTDTimeQualityToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"))),
      trackMTDMomentumToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMTDMomentumVMapTag"))),
      trackMTDPathLengthToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMTDPathLengthVMapTag"))),
      minTrackVtxWeight_(iConfig.getParameter<double>("minTrackVtxWeight")),
      minTrackTimeQuality_(iConfig.getParameter<double>("minTrackTimeQuality")),
      massPion_(iConfig.getParameter<double>("massPion")),
      massKaon_(iConfig.getParameter<double>("massKaon")),
      massProton_(iConfig.getParameter<double>("massProton")),
      probPion_(iConfig.getParameter<double>("probPion")),
      probKaon_(iConfig.getParameter<double>("probKaon")),
      probProton_(iConfig.getParameter<double>("probProton")),
      Tstart_(iConfig.getParameter<double>("Tstart")),
      coolingFactor_(iConfig.getParameter<double>("coolingFactor")),
      fixedT0Error_(iConfig.getParameter<double>("fixedT0Error"))	
{
  produces<edm::ValueMap<float>>(t0Name);
  produces<edm::ValueMap<float>>(sigmat0Name);
  produces<edm::ValueMap<float>>(t0safeName);
  produces<edm::ValueMap<float>>(sigmat0safeName);
  produces<edm::ValueMap<float>>(probPiName);
  produces<edm::ValueMap<float>>(probKName);
  produces<edm::ValueMap<float>>(probPName);
}

// Configuration descriptions
void TOFPIDProducer3D::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"))->setComment("Input tracks collection");
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"))
      ->setComment("Input ValueMap for track time at beamline");
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"))
      ->setComment("Input ValueMap for track time uncertainty at beamline");
  desc.add<edm::InputTag>("tofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"))
      ->setComment("Input ValueMap for track tof as kaon");
  desc.add<edm::InputTag>("tofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"))
      ->setComment("Input ValueMap for track tof as proton");
  desc.add<edm::InputTag>("vtxsSrc", edm::InputTag("unsortedOfflinePrimaryVertices3Dt"))
      ->setComment("Input primary vertex collection");
  desc.add<edm::InputTag>("trackMTDTimeVMapTag", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"))
      ->setComment("");
  desc.add<edm::InputTag>("trackMTDTimeErrorVMapTag", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"))
      ->setComment("");
  desc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"))
      ->setComment("");
  desc.add<edm::InputTag>("trackMTDMomentumVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackp"))
      ->setComment("");
  desc.add<edm::InputTag>("trackMTDPathLengthVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"))
      ->setComment("");

  desc.add<double>("minTrackVtxWeight", 0.5)->setComment("");
  desc.add<double>("minTrackTimeQuality", 0.8)->setComment("");

  desc.add<double>("massPion", 0.139570)->setComment("");
  desc.add<double>("massKaon", 0.493677)->setComment("");
  desc.add<double>("massProton", 0.938272)->setComment("");

  desc.add<double>("probPion", 0.7)->setComment("");
  desc.add<double>("probKaon", 0.2)->setComment("");
  desc.add<double>("probProton", 0.1)->setComment("");

  desc.add<double>("Tstart", 256.)->setComment("");
  desc.add<double>("coolingFactor", 0.5)->setComment("");

  desc.add<double>("fixedT0Error", 0.)->setComment("Use a fixed T0 uncertainty [ns]");

  descriptions.add("tofPIDProducer3D", desc);
}

template <class H, class T>
void TOFPIDProducer3D::fillValueMap(edm::Event& iEvent,
                                  const edm::Handle<H>& handle,
                                  const std::vector<T>& vec,
                                  const std::string& name) const {
  auto out = std::make_unique<edm::ValueMap<T>>();
  typename edm::ValueMap<T>::Filler filler(*out);
  filler.insert(handle, vec.begin(), vec.end());
  filler.fill();
  iEvent.put(std::move(out), name);
}

float TOFPIDProducer3D::trackTime(float const mass,
                                                  float const mtdTime,
                                                  float const mtdPathLength,
                                                  float const mtdMomentum) const {
  // speed of light, c = 29.9792458 cm/ns
  return (mtdTime - mtdPathLength * std::sqrt(1.f + mass * mass / mtdMomentum / mtdMomentum) / 29.9792458f);
}


void TOFPIDProducer3D::produce(edm::Event& ev, const edm::EventSetup& es) {
  edm::Handle<reco::TrackCollection> tracksH;
  ev.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  const auto& t0In = ev.get(t0Token_);

  const auto& sigmat0In = ev.get(sigmat0Token_);

  const auto& tofkIn = ev.get(tofkToken_);

  const auto& tofpIn = ev.get(tofpToken_);

  const auto& vtxs = ev.get(vtxsToken_);

  trackMTDTimes_ = ev.get(trackMTDTimeToken_);
  trackMTDTimeErrors_ = ev.get(trackMTDTimeErrorToken_);
  trackMTDTimeQualities_ = ev.get(trackMTDTimeQualityToken_);
  trackMTDMomenta_ = ev.get(trackMTDMomentumToken_);
  trackMTDPathLengths_ = ev.get(trackMTDPathLengthToken_);

  std::vector<float> t0OutRaw;
  std::vector<float> sigmat0OutRaw;
  std::vector<float> t0safeOutRaw;
  std::vector<float> sigmat0safeOutRaw;
  std::vector<float> probPiOutRaw;
  std::vector<float> probKOutRaw;
  std::vector<float> probPOutRaw;

  for (unsigned int itrack = 0; itrack < tracks.size(); ++itrack) {
    const reco::TrackRef trackref(tracksH, itrack);
    //std::cout << "NEW TRACK n. " << itrack << " dt = " << trackMTDTimes_[trackref] << std::endl; 
 
    float t0 = t0In[trackref];
    float t0safe = t0;
    float sigmat0safe = sigmat0In[trackref];
    float sigmatmtd = (trackMTDTimeErrors_[trackref] > 0. && fixedT0Error_ > 0.) ? fixedT0Error_ : trackMTDTimeErrors_[trackref];
    float sigmat0 = sigmatmtd;

    double tmtd = trackMTDTimes_[trackref];
    double t0_k = tmtd - tofkIn[trackref];
    double t0_p = tmtd - tofpIn[trackref];

    //std::cout << "Track t0_pi " << t0 << " t0_k " << t0_k << " t0_p " << t0_p << " track t0 best " << t0 << " t0safe " << t0safe << " sigmat0 " << sigmat0 << " sigmat0safe " << sigmat0safe << std::endl;

    float track_prob_PID_pi = -1;
    float track_prob_PID_k  = -1;
    float track_prob_PID_p  = -1;
    float t0_best  = t0;
  
    for (unsigned int ivtx = 0; ivtx < vtxs.size(); ++ivtx) {
     const reco::Vertex& vtx = vtxs[ivtx];
     float w = vtx.trackWeight(trackref);
     if (w > 0.5) {
  
      if (vtx.tracksSize()==0) {
        continue;
      }
      //std::cout << "NEW VERTEX " << std::endl;
    
      auto vtxTime(0.f), vtxTimeError(-1.f);
    
      auto const vtxTime_init = vtxTime;
      auto const vtxTimeError_init = vtxTimeError;
    
      double tsum = 0;
      double wsum = 0;
      double w2sum = 0;
      float prob_PID_pi = -1;
      float prob_PID_k  = -1;
      float prob_PID_p  = -1;
    
      double const a[3] = {probPion_, probKaon_, probProton_};
    
      std::vector<TrackInfo> v_trackInfo;
      v_trackInfo.reserve(vtx.tracksSize());
    
      // initial guess
      for (const auto& trk : vtx.tracks()) {
	reco::TrackRef trackreference = trk.castTo<reco::TrackRef>();
        auto const trkWeight = vtx.trackWeight(trk);
        if (trkWeight > minTrackVtxWeight_) {
          auto const trkTimeQuality = trackMTDTimeQualities_[trk];
    
          if (trkTimeQuality >= minTrackTimeQuality_) {
            auto const trkTime = trackMTDTimes_[trk];
            auto const trkTimeError = trackMTDTimeErrors_[trk];
            auto const trkPathLength = trackMTDPathLengths_[trk];
            auto const trkMomentum = trackMTDMomenta_[trk];
    
            v_trackInfo.emplace_back();
            auto& trkInfo = v_trackInfo.back();
    
            trkInfo.trackreference = trackreference;
            trkInfo.trkWeight = trkWeight;
            trkInfo.trkTime = trkTime;
            trkInfo.trkTimeError = trkTimeError;
    
            if (trkPathLength > 0) {
              trkInfo.trkTimeHyp[0] = trackTime(massPion_, trkTime, trkPathLength, trkMomentum);
              trkInfo.trkTimeHyp[1] = trackTime(massKaon_, trkTime, trkPathLength, trkMomentum);
              trkInfo.trkTimeHyp[2] = trackTime(massProton_, trkTime, trkPathLength, trkMomentum);
            } else {
              trkInfo.trkTimeHyp[0] = 0.f;
              trkInfo.trkTimeHyp[1] = 0.f;
              trkInfo.trkTimeHyp[2] = 0.f;
            }
    
            auto const wgt = trkWeight / (trkTimeError * trkTimeError);
            wsum += wgt;
    
            for (uint j = 0; j < 3; ++j) {
              tsum += wgt * trkInfo.trkTimeHyp[j] * a[j];
            }
            //std::cout << "vertexTimeFromTracks:     track"
            //    << " pt=" << track.pt() << " eta=" << track.eta() << " phi=" << track.phi()
            //    << " vtxWeight=" << trkWeight << " time=" << trkTime << " timeError=" << trkTimeError
            //    << " timeQuality=" << trkTimeQuality << " pathLength=" << trkPathLength << " momentum=" << trkMomentum
            //    << " timeHyp[pion]=" << trkInfo.trkTimeHyp[0] << " timeHyp[kaon]=" << trkInfo.trkTimeHyp[1]
            //    << " timeHyp[proton]=" << trkInfo.trkTimeHyp[2] << std::endl;
          }
        }
      }
      if (wsum > 0) {
        //std::cout << "vertexTimeFromTracks:   wsum = " << wsum << " tsum = " << tsum << " t0 = " << (wsum > 0 ? tsum / wsum : 0)
        //    << " trec = " << vtx.t()<< std::endl;
    
        auto t0_vtx = tsum / wsum;
        auto beta = 1. / Tstart_;
        int nit = 0;
        while ((nit++) < 100) {
          tsum = 0;
          wsum = 0;
          w2sum = 0;
    
          for (auto const& trkInfo : v_trackInfo) {
            double dt = trkInfo.trkTimeError;
            double e[3] = {0, 0, 0};
            double Z = exp(-beta * 0.5 * 3. * 3.);
    	    for (unsigned int j = 0; j < 3; j++) {
              auto const tpull = (trkInfo.trkTimeHyp[j] - t0_vtx) / dt;
              e[j] = exp(-0.5 * beta * tpull * tpull);
              Z += a[j] * e[j];
    	      double raw_prob_PID_pi = a[0] * e[0];
    	      double raw_prob_PID_k = a[1] * e[1];
    	      double raw_prob_PID_p = a[2] * e[2];
	      double normprob = 1. / (raw_prob_PID_pi + raw_prob_PID_k + raw_prob_PID_p);
  	      prob_PID_pi = raw_prob_PID_pi * normprob;
  	      prob_PID_k = raw_prob_PID_k * normprob;
  	      prob_PID_p = raw_prob_PID_p * normprob;
            }
            double wsum_trk = 0;
            for (uint j = 0; j < 3; j++) {
              double wt = a[j] * e[j] / Z;
              double w = wt * trkInfo.trkWeight / (dt * dt);
              wsum_trk += w;
              tsum += w * trkInfo.trkTimeHyp[j];
            }
    
            wsum += wsum_trk;
            w2sum += wsum_trk * wsum_trk * (dt * dt) / trkInfo.trkWeight;
  	    if(trkInfo.trackreference == trackref){
               track_prob_PID_pi = prob_PID_pi;
               track_prob_PID_k = prob_PID_k;
               track_prob_PID_p = prob_PID_p;
               //std::cout << nit << " prob_pi " << track_prob_PID_pi << " prob_k " << track_prob_PID_k << " prob_p " << track_prob_PID_p << std::endl;
  	    }
          }
          
    
          if (wsum < 1e-10) {
          //  std::cout << "vertexTimeFromTracks:   failed while iterating"<< std::endl;
            continue;
          }
    
          vtxTime = tsum / wsum;
    
          //std::cout << "vertexTimeFromTracks:   iteration=" << nit << ", T= " << 1 / beta << ", t=" << vtxTime
          //    << ", t-t0=" << vtxTime - t0_vtx<< std::endl;
    
          if ((std::abs(vtxTime - t0_vtx) < 1e-4 / std::sqrt(beta)) and beta >= 1.) {
            vtxTimeError = std::sqrt(w2sum) / wsum;
    
            //std::cout << "vertexTimeFromTracks:   tfit = " << vtxTime << " +/- " << vtxTimeError << " trec = " << vtx.t()
            //    << ", iteration=" << nit<< std::endl;
    
            break;
          }
    
          if ((std::abs(vtxTime - t0_vtx) < 1e-3) and beta < 1.) {
            beta = std::min(1., beta / coolingFactor_);
          }
    
          t0_vtx = vtxTime;
        }
    
        //std::cout << "vertexTimeFromTracks: failed to converge"<< std::endl;
      } else {
        //std::cout << "vertexTimeFromTracks: has no track timing info"<< std::endl;
      }
    
      vtxTime = vtxTime_init;
      vtxTimeError = vtxTimeError_init;
    
    } // RECO VTX LOOP
   } // IF VTX-tracks 
                    
   sigmat0safe = sigmatmtd;

   if(track_prob_PID_k > track_prob_PID_pi && track_prob_PID_k > track_prob_PID_p){
      t0_best = t0_k;
      t0safe = t0_k;
   }
   if(track_prob_PID_p > track_prob_PID_pi && track_prob_PID_p > track_prob_PID_k){
      t0_best = t0_p;
      t0safe = t0_p;
   }
   if( (1. - track_prob_PID_pi ) > 0.75 ) t0 = t0_best;
   
   //std::cout << "-Final prob_pi " << track_prob_PID_pi << " prob_k " << track_prob_PID_k << " prob_p " << track_prob_PID_p << std::endl;
  
   //std::cout << "- Final Track t0 best reassign " << t0 << " t0_k " << t0_k << " t0_p " << t0_p << " t0safe max " << t0safe << " sigmat0 " << sigmat0 << " sigmat0safe " << sigmat0safe << std::endl;

   t0OutRaw.push_back(t0);
   sigmat0OutRaw.push_back(sigmat0);
   t0safeOutRaw.push_back(t0safe);
   sigmat0safeOutRaw.push_back(sigmat0safe); 
   probPiOutRaw.push_back(track_prob_PID_pi);
   probKOutRaw.push_back(track_prob_PID_k);
   probPOutRaw.push_back(track_prob_PID_p);
  
  } // TRACK LOOP

fillValueMap(ev, tracksH, t0OutRaw, t0Name);
fillValueMap(ev, tracksH, sigmat0OutRaw, sigmat0Name);
fillValueMap(ev, tracksH, t0safeOutRaw, t0safeName);
fillValueMap(ev, tracksH, sigmat0safeOutRaw, sigmat0safeName);
fillValueMap(ev, tracksH, probPiOutRaw, probPiName);
fillValueMap(ev, tracksH, probKOutRaw, probKName);
fillValueMap(ev, tracksH, probPOutRaw, probPName);

}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TOFPIDProducer3D);
