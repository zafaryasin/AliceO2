#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TNtuple.h>

#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "ITStracking/Vertexer.h"

#if defined(__VERTEXER_ITS_GPU)
#include "GPUO2Interface.h"
#include "GPUReconstruction.h"
#include "GPUChainTracking.h"
#include "GPUChainITS.h"
#endif

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;
using namespace o2::gpu;

int run_primary_vertexer_ITS(const bool useGPU = false,
                             const bool useMCcheck = false,
                             const int inspEvt = -1,
                             const int numEvents = 1,
                             const std::string inputClustersITS = "o2clus_its.root",
                             const std::string inputGRP = "o2sim_grp.root",
                             const std::string simfilename = "o2sim.root",
                             const std::string paramfilename = "O2geometry.root",
                             const std::string path = "./")
{
  std::string outfile;

  if (useGPU) {
    if (useMCcheck) {
      outfile = "vertexer_gpu_data_MCCheck.root";
    } else {
      outfile = "vertexer_gpu_data.root";
    }
  } else {
    if (useMCcheck) {
      outfile = "vertexer_serial_data_MCCheck.root";
    } else {
      outfile = "vertexer_serial_data.root";
    }
  }

  const auto grp = o2::parameters::GRPObject::loadFrom(path + inputGRP);
  const bool isITS = grp->isDetReadOut(o2::detectors::DetID::ITS);
  const bool isContITS = grp->isDetContinuousReadOut(o2::detectors::DetID::ITS);
  std::cout << "ITS is in " << (isContITS ? "CONTINUOS" : "TRIGGERED") << " readout mode" << std::endl;
  TChain itsClusters("o2sim");
  itsClusters.AddFile((path + inputClustersITS).data());

  // Setup Runtime DB
  TFile paramFile((path + paramfilename).data());
  paramFile.Get("FAIRGeom");
  auto gman = o2::its::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
                                            o2::TransformType::L2G)); // request cached transforms

  // Get event header
  TChain mcHeaderTree("o2sim");
  mcHeaderTree.AddFile((path + simfilename).data());
  o2::dataformats::MCEventHeader* mcHeader = nullptr;
  if (!mcHeaderTree.GetBranch("MCEventHeader.")) {
    LOG(FATAL) << "Did not find MC event header in the input header file.";
  }
  mcHeaderTree.SetBranchAddress("MCEventHeader.", &mcHeader);

  // get clusters
  std::vector<o2::itsmft::Cluster>* clusters = nullptr;
  itsClusters.SetBranchAddress("ITSCluster", &clusters);

  TChain itsClustersROF("ITSClustersROF");
  itsClustersROF.AddFile((path + inputClustersITS).data());

  if (!itsClustersROF.GetBranch("ITSClustersROF")) {
    LOG(FATAL) << "Did not find ITS clusters branch ITSClustersROF in the input tree";
  }
  std::vector<o2::itsmft::ROFRecord>* rofs = nullptr;
  itsClustersROF.SetBranchAddress("ITSClustersROF", &rofs);
  itsClustersROF.GetEntry(0);

  // get labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels = nullptr;
  itsClusters.SetBranchAddress("ITSClusterMCTruth", &labels);

  TFile* outputfile = new TFile(outfile.data(), "recreate");

  TTree outTree("o2sim", "Vertexer Vertices");
  std::vector<o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>>* verticesITS =
    new std::vector<o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>>;
  outTree.Branch("ITSVertices", &verticesITS);

  // benchmarks
  TNtuple foundVerticesBenchmark("foundVerticesBenchmark", "Found vertices benchmark", "frameId:foundVertices");
  TNtuple timeBenchmark("timeBenchmark", "Time benchmarks", "init:trackletFinder:trackletMatcher:vertexFinder:total");
  // \benchmarks

  std::uint32_t roFrame = 0;

  //Settings
  o2::its::VertexingParameters parameters;
  // e.g. parameters.clusterContributorsCut = 5;
  //\Settings

  const int stopAt = (inspEvt == -1) ? rofs->size() : inspEvt + numEvents;
  const int startAt = (inspEvt == -1) ? 0 : inspEvt;

  std::unique_ptr<o2::its::VertexerTraits> traits;

#if defined(__VERTEXER_ITS_GPU)
  if (useGPU) {
    traits = o2::its::createVertexerTraitsGPU();
  }
#else
  if (!useGPU) {
    traits = std::make_unique<o2::its::VertexerTraits>();
  }
#endif

  o2::its::Vertexer vertexer(traits.get());
  vertexer.setParameters(parameters);
  itsClusters.GetEntry(0);
  mcHeaderTree.GetEntry(0);
  for (size_t iROfCount{static_cast<size_t>(startAt)}; iROfCount < static_cast<size_t>(stopAt); ++iROfCount) {
    auto& rof = (*rofs)[iROfCount];
    o2::its::ROframe frame(iROfCount); // to get meaningful roframeId
    std::cout << "ROframe: " << iROfCount << std::endl;
    int nclUsed = o2::its::ioutils::loadROFrameData(rof, frame, clusters, labels);

    std::array<float, 3> total{0.f, 0.f, 0.f};
    o2::its::ROframe* eventptr = &frame;

    // debug
    vertexer.setDebugTrackletSelection();
    vertexer.setDebugLines(); // Handle with care, takes very long
    vertexer.setDebugCombinatorics();
    vertexer.setDebugSummaryLines();
    // \debug

    total[0] = vertexer.evaluateTask(&o2::its::Vertexer::initialiseVertexer, "Vertexer initialisation", std::cout, eventptr);
    // total[1] = vertexer.evaluateTask(&o2::its::Vertexer::findTrivialMCTracklets, "Trivial Tracklet finding", std::cout); // If enable this, comment out the validateTracklets
    total[1] = vertexer.evaluateTask(&o2::its::Vertexer::findTracklets, "Tracklet finding", std::cout);
    if (useMCcheck) {
      vertexer.evaluateTask(&o2::its::Vertexer::filterMCTracklets, "MC tracklets filtering", std::cout);
    }
    total[2] = vertexer.evaluateTask(&o2::its::Vertexer::validateTracklets, "Adjacent tracklets validation", std::cout);
    total[3] = vertexer.evaluateTask(&o2::its::Vertexer::findVertices, "Vertex finding", std::cout);

    std::vector<Vertex> vertITS = vertexer.exportVertices();
    const size_t numVert = vertITS.size();
    foundVerticesBenchmark.Fill(static_cast<float>(iROfCount), static_cast<float>(numVert));
    verticesITS->swap(vertITS);
    // TODO: get vertexer postion form MC truth

    timeBenchmark.Fill(total[0], total[1], total[2], total[3], total[0] + total[1] + total[2] + total[3]);
    outTree.Fill();
  }
  traits.reset();
  outputfile->cd();
  outTree.Write();
  foundVerticesBenchmark.Write();
  timeBenchmark.Write();

  outputfile->Close();
  return 0;
}
#endif
