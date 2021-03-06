#include <OpenSim/OpenSim.h>
#include "SimTKsimbody.h"
#include "OpenSim/Tools/InverseKinematicsTool.h"
#include <OpenSim/Simulation/CoordinateReference.h>
#include <OpenSim/Simulation/InverseKinematicsSolver.h>
//#include <OpenSim/Simulation/Model/CoordinateLimitForce.h>
using namespace OpenSim;
using namespace SimTK;
using namespace std;

// my simulation tools
class SimTools {
public:
	SimTools();
	double Calc_rms_VERSION2(Storage& data_trc, double& initialTime, double& finalTime, Array<Array<double>>& pk_data);
	Array<double> ik_constrain_torso(OpenSim::Model& osimModel, string& marker_filename, string& ik_setup_filename, double& ti, double& tf);
	Array<Array<double>> RunSimulation_LIMITSTOP(OpenSim::Model& osimModel, Vector& PARAMS, double& initialTime, double& finalTime, Array<double>& ICs, const bool& save_states, string& fd, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h);
	double RunSimulation_wRMS(Storage& data_trc, OpenSim::Model& osimModel, Vector& PARAMS, double& initialTime, double& finalTime, Array<double>& ICs, const bool& save_states, string& fd, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h);
	void ApplyForce(Storage& force_storage,Model& osimModel);
	void SimTools::AddPKAtest(OpenSim::Model& osimModel);	
	double SimTools::RunSimulation_FlexExt(OpenSim::Model& osimModel, Vector& PARAMS, double& TiFlex, double& TiExt, double& TfFlex, double& TfExt, Array<double>& ICs_flex, Array<double>& ICs_ext, const bool& save_states, string& fd, Storage& force_storage_flex, Storage& force_storage_ext, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h,Storage& data_trc_Flex,Storage& data_trc_Ext);


Array<Array<double>> SimTools::DoIntegration(string& fd, Model& osimModel,Array<double>& ICs, double& initialTime, double& finalTime, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h);
private:
};

