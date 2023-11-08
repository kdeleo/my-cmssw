import FWCore.ParameterSet.Config as cms

def customise(process):
    process.mtdTracksValid.inputTagV = 'offlinePrimaryVertices'
    process.vertices4DValid.offline4DPV = 'offlinePrimaryVertices'
    return(process)
