// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "EMCALCalib/BadChannelMap.h"
#include "EMCALCalib/TimeCalibrationParams.h"
#include "EMCALCalib/TempCalibrationParams.h"
#include "EMCALCalib/TimeCalibParamL1Phase.h"
#include "EMCALCalib/TempCalibParamSM.h"
#include "EMCALCalib/CalibDB.h"

using namespace o2::emcal;

CalibDB::CalibDB(const std::string_view server) : CalibDB()
{
  mCCDBServer = server;
}

void CalibDB::init()
{
  mCCDBManager.init(mCCDBServer);
  mInit = true;
}

void CalibDB::storeBadChannelMap(BadChannelMap* bcm, const std::map<std::string, std::string>& metadata, ULong_t rangestart, ULong_t rangeend)
{
  if (!mInit)
    init();
  mCCDBManager.store(new o2::TObjectWrapper<o2::emcal::BadChannelMap>(bcm), "EMC/BadChannelMap", metadata, rangestart, rangeend);
}

void CalibDB::storeTimeCalibParam(TimeCalibrationParams* tcp, const std::map<std::string, std::string>& metadata, ULong_t rangestart, ULong_t rangeend)
{
  if (!mInit)
    init();
  mCCDBManager.store(new o2::TObjectWrapper<o2::emcal::TimeCalibrationParams>(tcp), "EMC/TimeCalibParams", metadata, rangestart, rangeend);
}

void CalibDB::storeTimeCalibParamL1Phase(TimeCalibParamL1Phase* tcp, const std::map<std::string, std::string>& metadata, ULong_t rangestart, ULong_t rangeend)
{
  if (!mInit)
    init();
  mCCDBManager.store(new o2::TObjectWrapper<o2::emcal::TimeCalibParamL1Phase>(tcp), "EMC/TimeCalibParamsL1Phase", metadata, rangestart, rangeend);
}

void CalibDB::storeTempCalibParam(TempCalibrationParams* tcp, const std::map<std::string, std::string>& metadata, ULong_t rangestart, ULong_t rangeend)
{
  if (!mInit)
    init();
  mCCDBManager.store(new o2::TObjectWrapper<o2::emcal::TempCalibrationParams>(tcp), "EMC/TempCalibParams", metadata, rangestart, rangeend);
}

void CalibDB::storeTempCalibParamSM(TempCalibParamSM* tcp, const std::map<std::string, std::string>& metadata, ULong_t rangestart, ULong_t rangeend)
{
  if (!mInit)
    init();
  mCCDBManager.store(new o2::TObjectWrapper<o2::emcal::TempCalibParamSM>(tcp), "EMC/TempCalibParamsSM", metadata, rangestart, rangeend);
}

void CalibDB::storeGainCalibFactors(GainCalibrationFactors* gcf, const std::map<std::string, std::string>& metadata, ULong_t rangestart, ULong_t rangeend)
{
  if (!mInit)
    init();
  mCCDBManager.store(new o2::TObjectWrapper<o2::emcal::GainCalibrationFactors>(gcf), "EMC/GainCalibFactors", metadata, rangestart, rangeend);
}

BadChannelMap* CalibDB::readBadChannelMap(ULong_t timestamp, const std::map<std::string, std::string>& metadata)
{
  if (!mInit)
    init();
  auto result = mCCDBManager.retrieve("EMC/BadChannelMap", metadata, timestamp);
  if (!result)
    throw ObjectNotFoundException(mCCDBServer, "EMC/BadChannelMap", metadata, timestamp);
  if (result->IsA() != TObjectWrapper<o2::emcal::BadChannelMap>::Class())
    throw TypeMismatchException("TObjectWrapper<o2::emcal::BadChannelMap>", result->IsA()->GetName());
  auto wrap = dynamic_cast<TObjectWrapper<o2::emcal::BadChannelMap>*>(result);
  if (!wrap)
    throw TypeMismatchException("TObjectWrapper<o2::emcal::BadChannelMap>", result->IsA()->GetName()); // type checked before - should not enter here
  return wrap->getObj();
}

TimeCalibrationParams* CalibDB::readTimeCalibParam(ULong_t timestamp, const std::map<std::string, std::string>& metadata)
{
  if (!mInit)
    init();
  auto result = mCCDBManager.retrieve("EMC/TimeCalibParams", metadata, timestamp);
  if (!result)
    throw ObjectNotFoundException(mCCDBServer, "EMC/TimeCalibParams", metadata, timestamp);
  if (result->IsA() != TObjectWrapper<o2::emcal::TimeCalibrationParams>::Class())
    throw TypeMismatchException("TObjectWrapper<o2::emcal::TimeCalibrationParams>", result->IsA()->GetName());
  auto wrap = dynamic_cast<TObjectWrapper<o2::emcal::TimeCalibrationParams>*>(result);
  if (!wrap)
    throw TypeMismatchException("TObjectWrapper<o2::emcal::TimeCalibrationParams>", result->IsA()->GetName()); // type checked before - should not enter here
  return wrap->getObj();
}

TimeCalibParamL1Phase* CalibDB::readTimeCalibParamL1Phase(ULong_t timestamp, const std::map<std::string, std::string>& metadata)
{
  if (!mInit)
    init();
  auto result = mCCDBManager.retrieve("EMC/TimeCalibParamsL1Phase", metadata, timestamp);
  if (!result)
    throw ObjectNotFoundException(mCCDBServer, "EMC/TimeCalibParamsL1Phase", metadata, timestamp);
  if (result->IsA() != TObjectWrapper<o2::emcal::TimeCalibParamL1Phase>::Class())
    throw TypeMismatchException("TObjectWrapper<o2::emcal::TimeCalibParamL1Phase>", result->IsA()->GetName());
  auto wrap = dynamic_cast<TObjectWrapper<o2::emcal::TimeCalibParamL1Phase>*>(result);
  if (!wrap)
    throw TypeMismatchException("TObjectWrapper<o2::emcal::TimeCalibParamL1Phase>", result->IsA()->GetName()); // type checked before - should not enter here
  return wrap->getObj();
}

TempCalibrationParams* CalibDB::readTempCalibParam(ULong_t timestamp, const std::map<std::string, std::string>& metadata)
{
  if (!mInit)
    init();
  auto result = mCCDBManager.retrieve("EMC/TempCalibParams", metadata, timestamp);
  if (!result)
    throw ObjectNotFoundException(mCCDBServer, "EMC/TempCalibParams", metadata, timestamp);
  if (result->IsA() != TObjectWrapper<o2::emcal::TempCalibrationParams>::Class())
    throw TypeMismatchException("TObjectWrapper<o2::emcal::TempCalibrationParams>", result->IsA()->GetName());
  auto wrap = dynamic_cast<TObjectWrapper<o2::emcal::TempCalibrationParams>*>(result);
  if (!wrap)
    throw TypeMismatchException("TObjectWrapper<o2::emcal::TempCalibrationParams>", result->IsA()->GetName()); // type checked before - should not enter here
  return wrap->getObj();
}

TempCalibParamSM* CalibDB::readTempCalibParamSM(ULong_t timestamp, const std::map<std::string, std::string>& metadata)
{
  if (!mInit)
    init();
  auto result = mCCDBManager.retrieve("EMC/TempCalibParamsSM", metadata, timestamp);
  if (!result)
    throw ObjectNotFoundException(mCCDBServer, "EMC/TempCalibParamsSM", metadata, timestamp);
  if (result->IsA() != TObjectWrapper<o2::emcal::TempCalibParamSM>::Class())
    throw TypeMismatchException("TObjectWrapper<o2::emcal::TempCalibParamSM>", result->IsA()->GetName());
  auto wrap = dynamic_cast<TObjectWrapper<o2::emcal::TempCalibParamSM>*>(result);
  if (!wrap)
    throw TypeMismatchException("TObjectWrapper<o2::emcal::TempCalibParamSM>", result->IsA()->GetName()); // type checked before - should not enter here
  return wrap->getObj();
}

GainCalibrationFactors* CalibDB::readGainCalibFactors(ULong_t timestamp, const std::map<std::string, std::string>& metadata)
{
  if (!mInit)
    init();
  auto result = mCCDBManager.retrieve("EMC/GainCalibFactors", metadata, timestamp);
  if (!result)
    throw ObjectNotFoundException(mCCDBServer, "EMC/GainCalibFactors", metadata, timestamp);
  if (result->IsA() != TObjectWrapper<o2::emcal::GainCalibrationFactors>::Class())
    throw TypeMismatchException("TObjectWrapper<o2::emcal::GainCalibrationFactors>", result->IsA()->GetName());
  auto wrap = dynamic_cast<TObjectWrapper<o2::emcal::GainCalibrationFactors>*>(result);
  if (!wrap)
    throw TypeMismatchException("TObjectWrapper<o2::emcal::GainCalibrationFactors>", result->IsA()->GetName()); // type checked before - should not enter here
  return wrap->getObj();
}
