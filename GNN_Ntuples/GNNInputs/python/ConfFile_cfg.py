import FWCore.ParameterSet.Config as cms

process = cms.Process("GNNInputs")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("file.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                               #'file:/home/kdeleo/Work/MTD_14_1_0_pre5/CMSSW_14_1_0_pre5/work/TTbar_PU200/29834.0_TTbar_14TeV+2026D110PU/step3.root'
                               'file:/home/kdeleo/Work/MTD_14_1_0_pre5/CMSSW_14_1_0_pre5/work/TTbar_noPU/step3.root'
                            )
)

process.gnnInputs = cms.EDAnalyzer('GNNInputs',
    inputTracks = cms.InputTag('generalTracks'),
    SimTag = cms.InputTag('mix', 'MergedTrackTruth'),
    TPtoRecoTrackAssoc = cms.InputTag('trackingParticleRecoTrackAsssociation'),
    pathLengthSrc = cms.InputTag('trackExtenderWithMTD', 'generalTrackPathLength'),
    momentumSrc = cms.InputTag('trackExtenderWithMTD', 'generalTrackp'),
    offlineBS = cms.InputTag('offlineBeamSpot'),
    tmtd = cms.InputTag('trackExtenderWithMTD', 'generalTracktmtd'),
    timeSrc = cms.InputTag('trackExtenderWithMTD', 'generalTracktmtd'),
    sigmaSrc = cms.InputTag('trackExtenderWithMTD', 'generalTracksigmatmtd'),
    t0PID = cms.InputTag('tofPID', 't0'),
    t0SafePID = cms.InputTag('tofPID', 't0safe'),
    sigmat0SafePID = cms.InputTag('tofPID', 'sigmat0safe'),
    trackMVAQual = cms.InputTag('mtdTrackQualityMVA', 'mtdQualMVA'),
    tofPi = cms.InputTag('trackExtenderWithMTD', 'generalTrackTofPi'),
    tofK = cms.InputTag('trackExtenderWithMTD', 'generalTrackTofK'),
    tofP = cms.InputTag('trackExtenderWithMTD', 'generalTrackTofP'),
    probPi = cms.InputTag('tofPID', 'probPi'),
    probK = cms.InputTag('tofPID', 'probK'),
    probP = cms.InputTag('tofPID', 'probP'),
    sigmatofpiSrc = cms.InputTag('trackExtenderWithMTD', 'generalTrackSigmaTofPi'),
    sigmatofkSrc = cms.InputTag('trackExtenderWithMTD', 'generalTrackSigmaTofK'),
    sigmatofpSrc = cms.InputTag('trackExtenderWithMTD', 'generalTrackSigmaTofP'),
    btlMatchChi2Src = cms.InputTag('trackExtenderWithMTD', 'btlMatchChi2'),
    btlMatchTimeChi2Src = cms.InputTag('trackExtenderWithMTD', 'btlMatchTimeChi2'),
    etlMatchChi2Src = cms.InputTag('trackExtenderWithMTD', 'etlMatchChi2'),
    etlMatchTimeChi2Src = cms.InputTag('trackExtenderWithMTD', 'etlMatchTimeChi2'),
    npixBarrelSrc = cms.InputTag('trackExtenderWithMTD', 'npixBarrel'),
    npixEndcapSrc = cms.InputTag('trackExtenderWithMTD', 'npixEndcap'),
    maxD0Significance = cms.double(4),
    maxD0Error = cms.double(1),
    maxDzError = cms.double(1),
    minPt = cms.double(0),
    maxEta = cms.double(2.4),
    maxNormalizedChi2 = cms.double(10),
    minSiliconLayersWithHits = cms.int32(5),
    minPixelLayersWithHits = cms.int32(2),
)

process.p = cms.Path(process.gnnInputs)
