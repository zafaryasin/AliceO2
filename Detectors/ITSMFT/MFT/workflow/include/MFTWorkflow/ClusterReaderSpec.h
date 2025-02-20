// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   ClusterReaderSpec.h

#ifndef O2_MFT_CLUSTERREADER_H_
#define O2_MFT_CLUSTERREADER_H_

#include "TFile.h"

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"

namespace o2
{
namespace mft
{

class ClusterReader : public o2::framework::Task
{
 public:
  ClusterReader() = default;
  ~ClusterReader() override = default;
  void init(o2::framework::InitContext& ic) final;
  void run(o2::framework::ProcessingContext& pc) final;

 private:
  int mState = 0;
  std::unique_ptr<TFile> mFile = nullptr;
};

/// create a processor spec
/// read simulated MFT digits from a root file
o2::framework::DataProcessorSpec getClusterReaderSpec();

} // namespace mft
} // namespace o2

#endif /* O2_MFT_CLUSTERREADER */
