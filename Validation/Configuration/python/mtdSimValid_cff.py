import FWCore.ParameterSet.Config as cms

# MTD validation sequences
from Validation.MtdValidation.btlSimHitsValid_cfi import btlSimHitsValid
from Validation.MtdValidation.btlDigiHitsValid_cfi import btlDigiHitsValid
from Validation.MtdValidation.btlLocalRecoValid_cfi import btlLocalRecoValid
from Validation.MtdValidation.etlLocalRecoValid_cfi import etlLocalRecoValid
from Validation.MtdValidation.etlSimHitsValid_cfi import etlSimHitsValid
from Validation.MtdValidation.etlDigiHitsValid_cfi import etlDigiHitsValid
from Validation.MtdValidation.mtdTracksValid_cfi import mtdTracksValid
from Validation.MtdValidation.vertices4DValid_cfi import vertices4DValid

#mtdTracksValid=mtdTracksValid.clone(
#   inputTagV='offlinePrimaryVertices',
#   t0SafePID='tofPID3D:t0safe',
#   sigmat0SafePID='tofPID3D:sigmat0safe',
#   sigmat0PID='tofPID3D:sigmat0',
#   t0PID='tofPID3D:t0'
#)
#
#vertices4DValid=vertices4DValid.clone(
#   offline4DPV='offlinePrimaryVertices',
#   t0PID='tofPID3D:t0',
#   t0SafePID='tofPID3D:t0safe',
#   sigmat0SafePID='tofPID3D:sigmat0safe',
#   probPi='tofPID3D:probPi',
#   probK='tofPID3D:probK',
#   probP='tofPID3D:probP'
#)

mtdSimValid  = cms.Sequence(btlSimHitsValid  + etlSimHitsValid )
mtdDigiValid = cms.Sequence(btlDigiHitsValid + etlDigiHitsValid)
mtdRecoValid = cms.Sequence(btlLocalRecoValid  + etlLocalRecoValid + mtdTracksValid + vertices4DValid)

