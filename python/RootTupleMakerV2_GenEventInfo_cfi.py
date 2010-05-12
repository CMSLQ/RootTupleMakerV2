import FWCore.ParameterSet.Config as cms

rootTupleGenEventInfo = cms.EDProducer("RootTupleMakerV2_GenEventInfo",
    GenEventInfoInputTag = cms.InputTag('generator'),
    StorePDFWeights      = cms.bool(True),
    PDFWeightsInputTag   = cms.InputTag('pdfWeights','cteq65')
)
