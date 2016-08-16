import FWCore.ParameterSet.Config as cms

rootTupleGenEventInfo = cms.EDProducer("RootTupleMakerV2_GenEventInfo",
    GenEventInfoInputTag = cms.InputTag('generator'),
    StorePDFWeights      = cms.bool(True),
    PDFCTEQWeightsInputTag   = cms.InputTag('pdfWeights','CT10nlo'),
    PDFPDF4LHCWeightsInputTag   = cms.InputTag('pdfWeights','PDF4LHC15'),
    PDFMMTHWeightsInputTag   = cms.InputTag('pdfWeights','MMHT2014nlo68cl'),
    PDFNNPDFWeightsInputTag   = cms.InputTag('pdfWeights','NNPDF30'),
    pileupInfo           = cms.InputTag('slimmedAddPileupInfo'),
    LHEEventProductInputTag   = cms.InputTag('externalLHEProducer'),
    LHERunInfoProductInputTag = cms.InputTag('externalLHEProducer')
    
)
