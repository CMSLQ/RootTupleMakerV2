import FWCore.ParameterSet.Config as cms

rootTupleGenEventInfo = cms.EDProducer("RootTupleMakerV2_GenEventInfo",
    GenEventInfoInputTag = cms.InputTag('generator'),
    StorePDFWeights      = cms.bool(True),
    PDFCTEQWeightsInputTag   = cms.InputTag('pdfWeights','CT10'),
    #PDFMSTWWeightsInputTag   = cms.InputTag('pdfWeights','MSTW2008nlo68cl'),
    PDFMSTWWeightsInputTag   = cms.InputTag('pdfWeights','MMHT2014nlo68cl'),
    #PDFNNPDFWeightsInputTag   = cms.InputTag('pdfWeights','NNPDF20'),
    PDFNNPDFWeightsInputTag   = cms.InputTag('pdfWeights','NNPDF30'),
    pileupInfo           = cms.InputTag('slimmedAddPileupInfo')                                       
)
