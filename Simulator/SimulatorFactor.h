#pragma once
#include "FiniteElementMethod\FemSimulator.h"
#include "mipc\MipcSimulator.h"

static BaseSimulator* createSimulator(std::string& simName, RunPlatform platform)
{
	if (simName == "fem_simulator")
	{
		return new FiniteElementMethod::FemSimulator(simName, platform);
	}
	else if (simName == "mipc_simulator")
	{
		return new MIPC::MipcSimulator(simName, platform);
	}
	else
	{
		return new BaseSimulator("base_simulator", platform);
	}
}

static bool readSimulatorFromConfigFile(BaseSimulator* sim, const std::string filename, TiXmlElement* item)
{
	std::string simName = sim->getSimulatorName();
	if (simName == "fem_simulator")
	{
		return dynamic_cast<FiniteElementMethod::FemSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "mipc_simulator")
	{
		return dynamic_cast<MIPC::MipcSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "base_simulator")
	{
		return sim->readSimulatorFromConfigFile(filename, item);
	}

	return false;
}
