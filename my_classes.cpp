#include "my_classes.h"

extern double bestSoFar;
extern int stepCount;


/* ################################################## */
// SIMULATION TOOLS OBJECTS
SimTools::SimTools(){}


double SimTools::Calc_rms_VERSION2(Storage& data_trc, double& initialTime, double& finalTime, Array<Array<double>>& pk_data)
{


	// Read trc storage object and insert data into 2D array for comparison with point kinematic data
	Array<Array<double>> trc_data(0,18);
	// NOTE: only comparing the markers on the head so the index is offset by 12 to account for the first
	// 4 markers in the trc file being the prescribed spine markers.
	for(int i =12;i<(12+3*6);i++)
	{
		data_trc.getDataColumn(i,trc_data[i-12]);
		//cout<<"\nDATA:    "<<trc_data[i-12]; //-> checked that data only to 5 decimal points (i.e. 100 microns) */
	}

	//cout<<"\n\npk data: "<<pk_data;
	//cout<<"\n\ntrc data: "<<trc_data;

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	// WE NOW HAVE 2 2D ARRAYS 
	//PK_DATA = DATA FROM SIMULATION (POINT KINEMATICS)
	//TRC_DATA = DATA FROM MEASURED TRC FILE.
	//NB. THEY SHOULD HAVE THE SAME NUMBER OF ARRAYS BUT ARRAYS IN EACH WILL BE OF DIFFERENT LENGTH
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //


	// Specify time values for each data set -> necessary for interploation
	Array<double> TRC_TIME, SIM_TIME; // array of time values from trc and simulation
	
	data_trc.getTimeColumn(TRC_TIME);
	SIM_TIME = pk_data[0];

	//cout<<"\n\npk data: "<<SIM_TIME;
	//cout<<"\n\ntrc data: "<<TRC_TIME;
	
	double rms = 0;
	double rms2 = 0;
	double rms_total = 0;

	// number of markers
	int no_markers = trc_data.getSize()/3;
	double ed = 0;
	
	Array<double> a(0.0,6);//);

	//cout<<"\n\ntrc data length:"<<trc_data[0].getSize();
	Array<Array<double>> ed_array(a,trc_data[0].getSize());

	// Create splines of all the simulation results so can compare with same position in trc file	
	OpenSim::GCVSplineSet *spline_set = new OpenSim::GCVSplineSet();
	GCVSpline m1x(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*0 + 0 + 1][0]);
	GCVSpline m1y(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*0 + 1 + 1][0]);
	GCVSpline m1z(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*0 + 2 + 1][0]);
	
	GCVSpline m2x(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*1 + 0 + 1][0]);
	GCVSpline m2y(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*1 + 1 + 1][0]);
	GCVSpline m2z(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*1 + 2 + 1][0]);
	
	GCVSpline m3x(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*2 + 0 + 1][0]);
	GCVSpline m3y(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*2 + 1 + 1][0]);
	GCVSpline m3z(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*2 + 2 + 1][0]);
	
	GCVSpline m4x(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*3 + 0 + 1][0]);
	GCVSpline m4y(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*3 + 1 + 1][0]);
	GCVSpline m4z(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*3 + 2 + 1][0]);
	
	GCVSpline m5x(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*4 + 0 + 1][0]);
	GCVSpline m5y(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*4 + 1 + 1][0]);
	GCVSpline m5z(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*4 + 2 + 1][0]);
	
	GCVSpline m6x(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*5 + 0 + 1][0]);
	GCVSpline m6y(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*5 + 1 + 1][0]);
	GCVSpline m6z(3, SIM_TIME.getSize(), &SIM_TIME[0], &pk_data[3*5 + 2 + 1][0]);

	spline_set->insert(0,m1x);
	spline_set->insert(1,m1y);
	spline_set->insert(2,m1z);
	spline_set->insert(3,m2x);
	spline_set->insert(4,m2y);
	spline_set->insert(5,m2z);
	spline_set->insert(6,m3x);
	spline_set->insert(7,m3y);
	spline_set->insert(8,m3z);
	spline_set->insert(9,m4x);
	spline_set->insert(10,m4y);
	spline_set->insert(11,m4z);
	spline_set->insert(12,m5x);
	spline_set->insert(13,m5y);
	spline_set->insert(14,m5z);
	spline_set->insert(15,m6x);
	spline_set->insert(16,m6y);
	spline_set->insert(17,m6z);
	
	// create temporary spline to use for comparison
	GCVSpline * spline = new GCVSpline();
	

	// for each time step 
	for (int i = 0; i < trc_data[0].getSize(); i++) {

		// for each marker at each time step
		for (int k=0;k<6;k++) {
			ed = 0; // initialise eucledian distance variable

			// for each coord for marker (i.e. x, y, z)
			for (int j = 0; j<3; j++) {
				
				// get spline from spline set which contains correct marker in simulation simulation
				spline = spline_set->getGCVSpline(3*k + j);
				
				// Obtain current time step in trc file -> note: spline.calcvalue requres vector input
				SimTK::Vector inputTime(1,TRC_TIME[i]);
				
				// Caclulate the square of the eucledican distance (i.e (x-x1)^2 + (y1-y2)^2 + (z1-z2)^2 )				
				ed += (trc_data[3*k + j][i] - spline->calcValue(inputTime))*(trc_data[3*k + j][i] - spline->calcValue(inputTime));

			}
			
			// return eucliedian distance between measurements and model i.e ed = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
			// Store in opensim array structure
			//if (isNaN(ed)){
			//	cout<<"ERROR ed is negative";
			//	break;
			//}
			ed_array[i][k] = sqrt(ed);
		
		}
	}

	//cout<<"\n\ned_array: "<<ed_array[0];//.getSize()<<"\t\t"<<ed_array.getSize();

	// Each column of ed_array consists of distance between markeri and model_markeri throughout time. 
	// I want to then get the rms of each column and store in rms_list (note: size 6 as 6 markers being included)
	Array<double> rms_list(0,6);

	// Loop through each marker column
	for (int j = 0; j<6; j++) {
		
		//cout<<"\nEUCLEDIAN ARRAY: ";	
		for (int i = 0; i<trc_data[0].getSize(); i++) {
			//cout<<ed_array[i][j]<<"\t";

			// accumatively sum the distance squared (sq.err)
			rms_list[j]+=pow(ed_array[i][j],(double)2);
		}
	}

	//cout<<"\n\nrms_list: "<<rms_list;
	// divide by length of trc file to get mean squared error and sum to return total mean squared error (rms_total)
	for (int i=0;i<rms_list.getSize();i++){
		rms_list[i] = rms_list[i]/(trc_data[0].getSize());
		rms_total += rms_list[i];
	}

	// Destruct pointers
	//spline->~GCVSpline();
	//spline_set->~GCVSplineSet();
	
	return rms_total;	

}

Array<double> SimTools::ik_constrain_torso(OpenSim::Model& osimModel, string& marker_filename, string& ik_setup_filename, double& ti, double& tf)
{

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	
	//// Do I need to include???
	//osimModel.invalidateSystem();


	// Initialise model for analysis
	SimTK::State& s = osimModel.initSystem();
	osimModel.updMultibodySystem().updDefaultState();
	

	//Convert old Tasks to references for assembly and tracking
	MarkersReference markersReference;
	Set<MarkerWeight> markerWeights;
	SimTK::Array_<CoordinateReference> coordinateReferences;

	// Create inverse kinematics tool using setup file
	OpenSim::InverseKinematicsTool* ik_tool = new InverseKinematicsTool(ik_setup_filename); 
	// only want to track @ initial time to provide IC's for simulation
	ik_tool->setStartTime(ti);
	ik_tool->setEndTime(ti+0.05);

	// Define reporter for output
	Kinematics kinematicsReporter(&osimModel);
	kinematicsReporter.setRecordAccelerations(false);
	kinematicsReporter.setInDegrees(false);

	// Loop through old "IKTaskSet" and assign weights to the coordinate and marker references
	// For coordinates, create the functions for coordinate reference values
	OpenSim::IKTaskSet ik_tasks = ik_tool->getIKTaskSet();
	int index = 0;
	for(int i=0; i < ik_tasks.getSize(); i++){
		if (!ik_tasks[i].getApply()) continue;

		else if(IKMarkerTask *markerTask = dynamic_cast<IKMarkerTask *>(&ik_tasks[i])){
			MarkerWeight *markerWeight = new MarkerWeight(markerTask->getName(), markerTask->getWeight());
			//cout<<"\n"<<markerTask->getWeight()<<"\n";
			//ERROR IN 3.0 //markerWeights..append(markerWeight)
			markerWeights.adoptAndAppend(markerWeight);
		}
	}

	//Set the weights for markers
	markersReference.setMarkerWeightSet(markerWeights);
	//Load the makers
	markersReference.loadMarkersFile(marker_filename);
	// marker names
	const SimTK::Array_<std::string>& markerNames =  markersReference.getNames();

	double start_time = ik_tool->getStartTime();
	double final_time = ik_tool->getEndTime();

	double accuracy = 0.00001;
	// create the solver given the input data
	double constraint_weight = SimTK::Infinity;
	InverseKinematicsSolver ikSolver(osimModel, markersReference, coordinateReferences, constraint_weight);
	ikSolver.setAccuracy(accuracy);
	s.updTime() = start_time;
	ikSolver.assemble(s);
	kinematicsReporter.begin(s);

	const clock_t start = clock();
	double dt = 1.0/markersReference.getSamplingFrequency();
	int Nframes = int((final_time-start_time)/dt)+1;
	AnalysisSet& analysisSet = osimModel.updAnalysisSet();
	analysisSet.begin(s);

	// number of markers
	int nm = markerWeights.getSize();
	// Initialise arrays of marker errors and marker xyz locations
	SimTK::Array_<double> squaredMarkerErrors(nm, 0.0);
	SimTK::Array_<Vec3> markerLocations(nm, Vec3(0));


	bool reportErrors = true;

	// Loop through each frame tracking the markers
	for (int i = 1; i <= Nframes; i++) {
		s.updTime() = start_time + (i-1)*dt;
		ikSolver.track(s);

		// Display errors
		if(reportErrors){
			double totalSquaredMarkerError = 0.0;
			double maxSquaredMarkerError = 0.0;
			int worst = -1;
			Array<double> error_report(0.0, nm+1);

			ikSolver.computeCurrentSquaredMarkerErrors(squaredMarkerErrors);

			for(int j=0; j<nm; ++j){
				totalSquaredMarkerError += squaredMarkerErrors[j];
				error_report.set(j,sqrt(squaredMarkerErrors[j]));

				if(squaredMarkerErrors[j] > maxSquaredMarkerError){
					maxSquaredMarkerError = squaredMarkerErrors[j];
					worst = j;
				}
			}
			error_report.set(error_report.getSize()-1,totalSquaredMarkerError);
			cout << "Frame " << i << " (t=" << s.getTime() << "):\t";
			cout << "total squared error = " << totalSquaredMarkerError;
			cout << ", marker error: RMS=" << sqrt(totalSquaredMarkerError/nm);
			cout << ", max=" << sqrt(maxSquaredMarkerError) << " (" << markerNames[worst] << ")" << endl;


			//TrackingErrors->append(s.getTime(), nm+1, &error_report[0]);

		}

		//if(reportMarkerLocations){
		//	ikSolver.computeCurrentMarkerLocations(markerLocations);
		//	Array<double> locations(0.0, 3*nm);
		//	for(int j=0; j<nm; ++j){
		//		for(int k=0; k<3; ++k)
		//			locations.set(3*j+k, markerLocations[j][k]);
		//	}

		//	modelMarkerLocations->append(s.getTime(), 3*nm, &locations[0]);

		//}

		kinematicsReporter.step(s, i);
		analysisSet.step(s, i);
	}	

	// Create storage object with results of tracking and set to be in radians
	Storage* k_storage = kinematicsReporter.getPositionStorage();
	k_storage->setInDegrees(false);

	//k_storage->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version3/kinematics.mot");

	// Open coord set to access names so can create IC array with results attributed to correct coordinate
	OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
	Array<string> rNames;
	CS.getNames(rNames);
	//cout<<"\n\nCS names: "<<rNames;
	//cout<<"\n\nsize of CS = "<<CS.getSize();


	Array<double> state_at_ti(0,CS.getSize());
	for(int i = 0; i<CS.getSize(); i++)
	{
		Array<double> a;
		k_storage->getDataColumn(rNames[i],a);
		//cout<<"\n\na: "<<a;
		
		for (int j = 0; j< a.getSize(); j++){
			state_at_ti[i] += a[j];
		}
		
		state_at_ti[i] = state_at_ti[i]/a.getSize();
		//cout<<"\n\nstate_at_ti[i]: "<<state_at_ti[i];
	} 

	// Sets torso coordinates to be set @ value specified by ik analysis (assuming torso does not move....)
	for(int i = 0; i<6; i++)
	{
		OpenSim::Function* function_ik = new OpenSim::Constant(state_at_ti[i]);
		// First 6 coordinates in cooridnate set correspond to rotx,roty,rotz,tx,ty,tz
		CS[i].setPrescribedFunction(*function_ik);
		CS[i].setDefaultValue(state_at_ti[i]);
		CS[i].setDefaultIsPrescribed(true);
		

	} 
	
	
	//osimModel.invalidateSystem();

	// return array of IC's
	return (state_at_ti);
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
}

Array<Array<double>> SimTools::RunSimulation_LIMITSTOP(OpenSim::Model& osimModel, Vector& PARAMS, double& initialTime, double& finalTime, Array<double>& ICs, const bool& save_states, string& fd, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h)//, Array<double>& times)
	//const double& ti, const double& tf, const OpenSim::Storage& st)
{

	// in validate model so model can be changed
	osimModel.invalidateSystem();

	
	double ThetaStar = PARAMS[0];
	Vec3 k1(PARAMS[1]);
	double k2 = PARAMS[2]*SimTK::Pi/180;
	Vec3 Rdamping(PARAMS[3]);

	Vec3 Tk1(0); // translational dof stiffness
	Vec3 Tdamping(0);//PARAMS[3]); // translational dof damping
	

	//// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	//// ADD BUSHING FORCES WITH SPECIFIED PROPERTIES
	//// Define position in body 1 (p1), position in body 2 (p2), orientation in body 1 (o1), orientation in  body 2 (o2)
	//// These are position and orientation os joints wrt bodies. Can set to be constant throughout...
	Vec3 p1bushing(0),p2bushing(0),o1(0),o2(0);
	double transition = 2.0;

	double vert_mass = 0.018;
	double head_mass = 0.5;
	
	////Remove forces from osimModel so correct number of forces remain in the model
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);

	
	OpenSim::Force *limit1 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("t1Tot2_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit1);
	OpenSim::Force *limit2 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c7Tot1_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit2);
	OpenSim::Force *limit3 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c6Toc7_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit3);
	OpenSim::Force *limit4 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c5Toc6_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit4);
	OpenSim::Force *limit5 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c4Toc5_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit5);
	OpenSim::Force *limit6 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c3Toc4_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit6);
	OpenSim::Force *limit7 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c2Toc3_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit7);
	OpenSim::Force *limit8 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c1Toc2_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit8);
	OpenSim::Force *limit9 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("skToc1_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit9);
	
	//k2 = 50000;
	//transition = 1e-6;
	//ThetaStar = 10;
	//Rdamping[0] = 1000;
	//double k2a = 10;
	//OpenSim::Force *limit1s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("t1Tot2_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit1s);
	//OpenSim::Force *limit2s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c7Tot1_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit2s);
	//OpenSim::Force *limit3s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c6Toc7_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit3s);
	//OpenSim::Force *limit4s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c5Toc6_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit4s);
	//OpenSim::Force *limit5s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c4Toc5_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit5s);
	//OpenSim::Force *limit6s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c3Toc4_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit6s);
	//OpenSim::Force *limit7s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c2Toc3_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit7s);
	//OpenSim::Force *limit8s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c1Toc2_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit8s);
	//OpenSim::Force *limit9s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("skToc1_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit9s);

	//Create new bushing forces and add to model
	//t2-t1 joint
	p1bushing = osimModel.updJointSet().get("t1Tot2").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("t1Tot2").upd_location();
	o1 = osimModel.updJointSet().get("t1Tot2").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("t1Tot2").upd_orientation();
	OpenSim::BushingForce* t2t1 = new OpenSim::BushingForce("t2",p1bushing,o1,"t1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	t2t1->setName("t2t1_bushing");
	osimModel.addForce(t2t1);
	
	p1bushing = osimModel.updJointSet().get("c7Tot1").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("c7Tot1").upd_location();
	o1 = osimModel.updJointSet().get("c7Tot1").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("c7Tot1").upd_orientation();
	//t1-c7 joint
	OpenSim::BushingForce* t1c7 = new OpenSim::BushingForce("t1",p1bushing,o1,"c7",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	t1c7->setName("t1c7_bushing");
	osimModel.addForce(t1c7);
	
	p1bushing = osimModel.updJointSet().get("c6Toc7").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("c6Toc7").upd_location();
	o1 = osimModel.updJointSet().get("c6Toc7").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("c6Toc7").upd_orientation();
	//c7-c6 joint
	OpenSim::BushingForce* c7c6 = new OpenSim::BushingForce("c7",p1bushing,o1,"c6",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c7c6->setName("c7c6_bushing");
	osimModel.addForce(c7c6);

	p1bushing = osimModel.updJointSet().get("c5Toc6").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("c5Toc6").upd_location();
	o1 = osimModel.updJointSet().get("c5Toc6").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("c5Toc6").upd_orientation();
	//c6-c5 joint
	OpenSim::BushingForce* c6c5 = new OpenSim::BushingForce("c6",p1bushing,o1,"c5",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c6c5->setName("c6c5_bushing");
	osimModel.addForce(c6c5);

	p1bushing = osimModel.updJointSet().get("c4Toc5").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("c4Toc5").upd_location();
	o1 = osimModel.updJointSet().get("c4Toc5").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("c4Toc5").upd_orientation();
	//c5-c4 joint
	OpenSim::BushingForce* c5c4 = new OpenSim::BushingForce("c5",p1bushing,o1,"c4",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c5c4->setName("c5c4_bushing");
	osimModel.addForce(c5c4);

	p1bushing = osimModel.updJointSet().get("c3Toc4").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("c3Toc4").upd_location();
	o1 = osimModel.updJointSet().get("c3Toc4").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("c3Toc4").upd_orientation();
	//c4-c3 joint
	OpenSim::BushingForce* c4c3 = new OpenSim::BushingForce("c4",p1bushing,o1,"c3",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c4c3->setName("c4c3_bushing");
	osimModel.addForce(c4c3);

	p1bushing = osimModel.updJointSet().get("c2Toc3").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("c2Toc3").upd_location();
	o1 = osimModel.updJointSet().get("c2Toc3").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("c2Toc3").upd_orientation();
	//c3-c2 joint
	OpenSim::BushingForce* c3c2 = new OpenSim::BushingForce("c3",p1bushing,o1,"c2",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c3c2->setName("c3c2_bushing");
	osimModel.addForce(c3c2);

	p1bushing = osimModel.updJointSet().get("c1Toc2").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("c1Toc2").upd_location();
	o1 = osimModel.updJointSet().get("c1Toc2").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("c1Toc2").upd_orientation();
	//c2-c1 joint
	OpenSim::BushingForce* c2c1 = new OpenSim::BushingForce("c2",p1bushing,o1,"c1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c2c1->setName("c2c1_bushing");
	osimModel.addForce(c2c1);

	p1bushing = osimModel.updJointSet().get("skToc1").upd_location_in_parent();
	p2bushing = osimModel.updJointSet().get("skToc1").upd_location();
	o1 = osimModel.updJointSet().get("skToc1").upd_orientation_in_parent();
	o2 = osimModel.updJointSet().get("skToc1").upd_orientation();
	//c1-sk joint
	OpenSim::BushingForce* c1sk = new OpenSim::BushingForce("c1",p1bushing,o1,"sk",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c1sk->setName("c1sk_bushing");
	osimModel.addForce(c1sk);		


	// mass properties
	osimModel.updBodySet().get("sk").setMass(head_mass);
	osimModel.updBodySet().get("t2").setMass(vert_mass);
	osimModel.updBodySet().get("t1").setMass(vert_mass);
	osimModel.updBodySet().get("c7").setMass(vert_mass);
	osimModel.updBodySet().get("c6").setMass(vert_mass);
	osimModel.updBodySet().get("c5").setMass(vert_mass);
	osimModel.updBodySet().get("c4").setMass(vert_mass);
	osimModel.updBodySet().get("c3").setMass(vert_mass);
	osimModel.updBodySet().get("c2").setMass(vert_mass);
	osimModel.updBodySet().get("c1").setMass(vert_mass);


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version8/Version8_constrainedwFORCES.osim");
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	
	// Initialise and run simulation
	SimTK::State& s = osimModel.initSystem();
	osimModel.updMultibodySystem().updDefaultState();



	// Create the integrator, force reporter, and manager for the simulation.
	SimTK::RungeKuttaMersonIntegrator integrator(osimModel.getMultibodySystem());
	// default settings
	integrator.setAccuracy(0.00001);
	integrator.setMaximumStepSize(1);
	integrator.setMinimumStepSize(1e-8);
	
	
	// Create the force reporter
	// ForceReporter* reporter = new ForceReporter(&osimModel);
	// osimModel.addAnalysis(reporter);

	// Create the manager
	Manager manager(osimModel,  integrator);

	// Print out details of the model
	//osimModel.printDetailedInfo(s, std::cout);

	// Define non-zero (defaults are 0) states for the free joint
	// States are from initial ik step. States are passed via ICs
	OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
	//cout<<"\n\nCS = "<<CS;

	for (int i = 0; i<CS.getSize(); i++)//
	{
		//nb first 6 coords are prescribed (torso assumed to be fixed)
		CS[i].setValue(s,ICs[i]);
		CS[i].setSpeedValue(s,0);
		
	}

	// Print out the initial position and velocity states
	// s.getQ().dump("Initial q's"); // block positions
	// s.getU().dump("Initial u's"); // block velocities
	// std::cout << "Initial time: " << s.getTime() << std::endl;
	

	// Integrate from initial time to final time
	manager.setInitialTime(initialTime);
	manager.setFinalTime(finalTime);
	std::cout<<"\n\nIntegrating from "<<initialTime<<" to "<<finalTime<<std::endl;

	manager.integrate(s);

	cout<<"\n\nintegrating done\n\n";
	
	// return pk analysis results
	int no_markers = osimModel.getMarkerSet().getSize();

	// Declare 2D array with all pk analysis results in it. Can then directly compare with trc data
	Array<Array<double>> a(0,3*no_markers+1);
	double resample_dt = 0.01;

	// Store marker 1 simulation results
	Storage* temp = m1h.getPositionStorage();
	Array<string> col_labels("",4);// = new Array<string>();
	col_labels[0] = "time";//,"m2","m3");		
	col_labels[1] = "m_x";
	col_labels[2] = "m_y";
	col_labels[3] = "m_z";
	temp->setColumnLabels(col_labels);

	temp->getTimeColumn(a[0]);
	
	temp->getDataColumn("m_x",a[1]);
	temp->getDataColumn("m_y",a[2]);
	temp->getDataColumn("m_z",a[3]);

	//cout<<"\n\na = "<<a[0];
	//cout<<"\n\na = "<<a[1];
	//cout<<"\n\na = "<<a[2];
	//cout<<"\n\na = "<<a[3];

	// Store marker 2 simulation results
	temp = m2h.getPositionStorage();
	//temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[4]);
	temp->getDataColumn("m_y",a[5]);
	temp->getDataColumn("m_z",a[6]);

	//cout<<"\n\na = "<<a[4];
	//cout<<"\n\na = "<<a[5];
	//cout<<"\n\na = "<<a[6];

	// Store marker 3 simulation results
	temp = m3h.getPositionStorage();
	//tempm3->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[7]);
	temp->getDataColumn("m_y",a[8]);
	temp->getDataColumn("m_z",a[9]);

	//cout<<"\n\na = "<<a[7];
	//cout<<"\n\na = "<<a[8];
	//cout<<"\n\na = "<<a[9];

	// Store marker 4 simulation results
	temp = m4h.getPositionStorage();
	//tempm3->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[10]);
	temp->getDataColumn("m_y",a[11]);
	temp->getDataColumn("m_z",a[12]);

	//cout<<"\n\na = "<<a[10];
	//cout<<"\n\na = "<<a[11];
	//cout<<"\n\na = "<<a[12];

	// Store marker 5 simulation results
	temp = m5h.getPositionStorage();
	//tempm3->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[13]);
	temp->getDataColumn("m_y",a[14]);
	temp->getDataColumn("m_z",a[15]);

	//cout<<"\n\na = "<<a[13];
	//cout<<"\n\na = "<<a[14];
	//cout<<"\n\na = "<<a[15];

	// Store marker 6 simulation results
	temp = m6h.getPositionStorage();
	/*temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version6/Model/Marker2_positionv2.sto");*/
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[16]);
	temp->getDataColumn("m_y",a[17]);
	temp->getDataColumn("m_z",a[18]);

	//cout<<"\n\na = "<<a[16];
	//cout<<"\n\na = "<<a[17];
	//cout<<"\n\na = "<<a[18];

	// Save simulation results to file if specified
	if (save_states)
	{
		// Save the states
		Storage statesDegrees(manager.getStateStorage());
		statesDegrees.resampleLinear(resample_dt);

		statesDegrees.print(fd);
		//statesDegrees.print(fd.append("Fwd_simulation_states.sto"));
		//osimModel.updSimbodyEngine().convertRadiansToDegrees(statesDegrees);
		//statesDegrees.setWriteSIMMHeader(true);
		//statesDegrees.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Fwd_simulation_states_degrees.mot");	
	}

	// Return 2D array
	return (a);

}

double SimTools::RunSimulation_wRMS(Storage& data_trc, OpenSim::Model& osimModel, Vector& PARAMS, double& initialTime, double& finalTime, Array<double>& ICs, const bool& save_states, string& fd, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h)//, Array<double>& times)
	//const double& ti, const double& tf, const OpenSim::Storage& st)
{

	// in validate model so model can be changed
	osimModel.invalidateSystem();

	
	double ThetaStar = PARAMS[0];
	Vec3 k1(PARAMS[1]);
	double k2 = PARAMS[2]*SimTK::Pi/180;
	Vec3 Rdamping(10);//PARAMS[3]);

	//double ThetaStar_sk = PARAMS[4];
	//Vec3 k1_sk(PARAMS[5]);
	//double k2_sk = PARAMS[6]*SimTK::Pi/180;
	

	Vec3 Tk1(0); // translational dof stiffness
	Vec3 Tdamping(0);//PARAMS[3]); // translational dof damping
	

	//// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	//// ADD BUSHING FORCES WITH SPECIFIED PROPERTIES
	//// Define position in body 1 (p1), position in body 2 (p2), orientation in body 1 (o1), orientation in  body 2 (o2)
	//// These are position and orientation os joints wrt bodies. Can set to be constant throughout...
	Vec3 p1bushing(0),p2bushing(0),o1(0),o2(0);
	double transition = 2.0;

	double vert_mass = 0.018;
	double head_mass = 0.45;//PARAMS[4];
	
	////Remove forces from osimModel so correct number of forces remain in the model
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);
	osimModel.updForceSet().remove(1);

	
	OpenSim::Force *limit1 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("t1Tot2_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit1);
	OpenSim::Force *limit2 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c7Tot1_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit2);
	OpenSim::Force *limit3 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c6Toc7_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit3);
	OpenSim::Force *limit4 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c5Toc6_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit4);
	OpenSim::Force *limit5 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c4Toc5_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit5);
	OpenSim::Force *limit6 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c3Toc4_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit6);
	OpenSim::Force *limit7 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c2Toc3_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit7);
	OpenSim::Force *limit8 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c1Toc2_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit8);
	OpenSim::Force *limit9 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("skToc1_FE").getName(),ThetaStar,k2,-ThetaStar,k2,0,transition);
	osimModel.addForce(limit9);
	
	//k2 = 50000;
	//transition = 1e-6;
	//ThetaStar = 10;
	//Rdamping[0] = 1000;
	//double k2a = 10;
	//OpenSim::Force *limit1s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("t1Tot2_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit1s);
	//OpenSim::Force *limit2s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c7Tot1_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit2s);
	//OpenSim::Force *limit3s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c6Toc7_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit3s);
	//OpenSim::Force *limit4s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c5Toc6_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit4s);
	//OpenSim::Force *limit5s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c4Toc5_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit5s);
	//OpenSim::Force *limit6s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c3Toc4_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit6s);
	//OpenSim::Force *limit7s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c2Toc3_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit7s);
	//OpenSim::Force *limit8s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c1Toc2_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit8s);
	//OpenSim::Force *limit9s = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("skToc1_stretch").getName(),0,k2,-ThetaStar,k2a,Rdamping[0],transition);
	//osimModel.addForce(limit9s);


	////////////////////////////////////////////////////////////////////////////////////////////
	//// CORRECT BUSHING PROPERTIES
	////Create new bushing forces and add to model
	////t2-t1 joint
	//p1bushing = osimModel.updJointSet().get("t1Tot2").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("t1Tot2").upd_location();
	//o1 = osimModel.updJointSet().get("t1Tot2").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("t1Tot2").upd_orientation();
	//OpenSim::BushingForce* t2t1 = new OpenSim::BushingForce("t2",p1bushing,o1,"t1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//t2t1->setName("t2t1_bushing");
	//osimModel.addForce(t2t1);
	//
	//p1bushing = osimModel.updJointSet().get("c7Tot1").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c7Tot1").upd_location();
	//o1 = osimModel.updJointSet().get("c7Tot1").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("c7Tot1").upd_orientation();
	////t1-c7 joint
	//OpenSim::BushingForce* t1c7 = new OpenSim::BushingForce("t1",p1bushing,o1,"c7",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//t1c7->setName("t1c7_bushing");
	//osimModel.addForce(t1c7);
	//
	//p1bushing = osimModel.updJointSet().get("c6Toc7").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c6Toc7").upd_location();
	//o1 = osimModel.updJointSet().get("c6Toc7").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("c6Toc7").upd_orientation();
	////c7-c6 joint
	//OpenSim::BushingForce* c7c6 = new OpenSim::BushingForce("c7",p1bushing,o1,"c6",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//c7c6->setName("c7c6_bushing");
	//osimModel.addForce(c7c6);

	//p1bushing = osimModel.updJointSet().get("c5Toc6").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c5Toc6").upd_location();
	//o1 = osimModel.updJointSet().get("c5Toc6").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("c5Toc6").upd_orientation();
	////c6-c5 joint
	//OpenSim::BushingForce* c6c5 = new OpenSim::BushingForce("c6",p1bushing,o1,"c5",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//c6c5->setName("c6c5_bushing");
	//osimModel.addForce(c6c5);

	//p1bushing = osimModel.updJointSet().get("c4Toc5").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c4Toc5").upd_location();
	//o1 = osimModel.updJointSet().get("c4Toc5").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("c4Toc5").upd_orientation();
	////c5-c4 joint
	//OpenSim::BushingForce* c5c4 = new OpenSim::BushingForce("c5",p1bushing,o1,"c4",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//c5c4->setName("c5c4_bushing");
	//osimModel.addForce(c5c4);

	//p1bushing = osimModel.updJointSet().get("c3Toc4").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c3Toc4").upd_location();
	//o1 = osimModel.updJointSet().get("c3Toc4").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("c3Toc4").upd_orientation();
	////c4-c3 joint
	//OpenSim::BushingForce* c4c3 = new OpenSim::BushingForce("c4",p1bushing,o1,"c3",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//c4c3->setName("c4c3_bushing");
	//osimModel.addForce(c4c3);

	//p1bushing = osimModel.updJointSet().get("c2Toc3").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c2Toc3").upd_location();
	//o1 = osimModel.updJointSet().get("c2Toc3").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("c2Toc3").upd_orientation();
	////c3-c2 joint
	//OpenSim::BushingForce* c3c2 = new OpenSim::BushingForce("c3",p1bushing,o1,"c2",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//c3c2->setName("c3c2_bushing");
	//osimModel.addForce(c3c2);

	//p1bushing = osimModel.updJointSet().get("c1Toc2").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c1Toc2").upd_location();
	//o1 = osimModel.updJointSet().get("c1Toc2").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("c1Toc2").upd_orientation();
	////c2-c1 joint
	//OpenSim::BushingForce* c2c1 = new OpenSim::BushingForce("c2",p1bushing,o1,"c1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//c2c1->setName("c2c1_bushing");
	//osimModel.addForce(c2c1);

	//p1bushing = osimModel.updJointSet().get("skToc1").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("skToc1").upd_location();
	//o1 = osimModel.updJointSet().get("skToc1").upd_orientation_in_parent();
	//o2 = osimModel.updJointSet().get("skToc1").upd_orientation();
	////c1-sk joint
	//OpenSim::BushingForce* c1sk = new OpenSim::BushingForce("c1",p1bushing,o1,"sk",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	//c1sk->setName("c1sk_bushing");
	//osimModel.addForce(c1sk);		
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Create new bushing forces and add to model
	o1 = Vec3(0,0,PARAMS[3]);

	//t2-t1 joint
	//p1bushing = osimModel.updJointSet().get("t1Tot2").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("t1Tot2").upd_location();
	//o1 = -osimModel.updJointSet().get("t1Tot2").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("t1Tot2").upd_orientation();
	OpenSim::BushingForce* t2t1 = new OpenSim::BushingForce("t2",p1bushing,o1,"t1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	t2t1->setName("t2t1_bushing");
	osimModel.addForce(t2t1);
	
	//p1bushing = osimModel.updJointSet().get("c7Tot1").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c7Tot1").upd_location();
	////o1 = -osimModel.updJointSet().get("c7Tot1").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("c7Tot1").upd_orientation();
	//t1-c7 joint
	OpenSim::BushingForce* t1c7 = new OpenSim::BushingForce("t1",p1bushing,o1,"c7",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	t1c7->setName("t1c7_bushing");
	osimModel.addForce(t1c7);
	
	//p1bushing = osimModel.updJointSet().get("c6Toc7").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c6Toc7").upd_location();
	////o1 = -osimModel.updJointSet().get("c6Toc7").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("c6Toc7").upd_orientation();
	//c7-c6 joint
	OpenSim::BushingForce* c7c6 = new OpenSim::BushingForce("c7",p1bushing,o1,"c6",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c7c6->setName("c7c6_bushing");
	osimModel.addForce(c7c6);

	//p1bushing = osimModel.updJointSet().get("c5Toc6").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c5Toc6").upd_location();
	////o1 = -osimModel.updJointSet().get("c5Toc6").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("c5Toc6").upd_orientation();
	//c6-c5 joint
	OpenSim::BushingForce* c6c5 = new OpenSim::BushingForce("c6",p1bushing,o1,"c5",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c6c5->setName("c6c5_bushing");
	osimModel.addForce(c6c5);

	//p1bushing = osimModel.updJointSet().get("c4Toc5").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c4Toc5").upd_location();
	////o1 = -osimModel.updJointSet().get("c4Toc5").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("c4Toc5").upd_orientation();
	//c5-c4 joint
	OpenSim::BushingForce* c5c4 = new OpenSim::BushingForce("c5",p1bushing,o1,"c4",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c5c4->setName("c5c4_bushing");
	osimModel.addForce(c5c4);

	//p1bushing = osimModel.updJointSet().get("c3Toc4").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c3Toc4").upd_location();
	////o1 = -osimModel.updJointSet().get("c3Toc4").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("c3Toc4").upd_orientation();
	//c4-c3 joint
	OpenSim::BushingForce* c4c3 = new OpenSim::BushingForce("c4",p1bushing,o1,"c3",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c4c3->setName("c4c3_bushing");
	osimModel.addForce(c4c3);

	//p1bushing = osimModel.updJointSet().get("c2Toc3").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c2Toc3").upd_location();
	////o1 = -osimModel.updJointSet().get("c2Toc3").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("c2Toc3").upd_orientation();
	//c3-c2 joint
	OpenSim::BushingForce* c3c2 = new OpenSim::BushingForce("c3",p1bushing,o1,"c2",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c3c2->setName("c3c2_bushing");
	osimModel.addForce(c3c2);

	//p1bushing = osimModel.updJointSet().get("c1Toc2").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("c1Toc2").upd_location();
	////o1 = -osimModel.updJointSet().get("c1Toc2").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("c1Toc2").upd_orientation();
	//c2-c1 joint
	OpenSim::BushingForce* c2c1 = new OpenSim::BushingForce("c2",p1bushing,o1,"c1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c2c1->setName("c2c1_bushing");
	osimModel.addForce(c2c1);

	//p1bushing = osimModel.updJointSet().get("skToc1").upd_location_in_parent();
	//p2bushing = osimModel.updJointSet().get("skToc1").upd_location();
	////o1 = -osimModel.updJointSet().get("skToc1").upd_orientation_in_parent();
	//o2 = -osimModel.updJointSet().get("skToc1").upd_orientation();
	//c1-sk joint
	OpenSim::BushingForce* c1sk = new OpenSim::BushingForce("c1",p1bushing,o1,"sk",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c1sk->setName("c1sk_bushing");
	osimModel.addForce(c1sk);		


	// mass properties
	osimModel.updBodySet().get("sk").setMass(head_mass);
	osimModel.updBodySet().get("t2").setMass(vert_mass);
	osimModel.updBodySet().get("t1").setMass(vert_mass);
	osimModel.updBodySet().get("c7").setMass(vert_mass);
	osimModel.updBodySet().get("c6").setMass(vert_mass);
	osimModel.updBodySet().get("c5").setMass(vert_mass);
	osimModel.updBodySet().get("c4").setMass(vert_mass);
	osimModel.updBodySet().get("c3").setMass(vert_mass);
	osimModel.updBodySet().get("c2").setMass(vert_mass);
	osimModel.updBodySet().get("c1").setMass(vert_mass);


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version8/Version8_constrainedwFORCES_FLEX.osim");
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	
	// Initialise and run simulation
	SimTK::State& s = osimModel.initSystem();
	osimModel.updMultibodySystem().updDefaultState();



	// Create the integrator, force reporter, and manager for the simulation.
	SimTK::RungeKuttaMersonIntegrator integrator(osimModel.getMultibodySystem());
	// default settings
	integrator.setAccuracy(0.00001);
	integrator.setMaximumStepSize(1);
	integrator.setMinimumStepSize(1e-8);
	
	
	// Create the force reporter
	// ForceReporter* reporter = new ForceReporter(&osimModel);
	// osimModel.addAnalysis(reporter);

	// Create the manager
	Manager manager(osimModel,  integrator);

	// Print out details of the model
	//osimModel.printDetailedInfo(s, std::cout);

	// Define non-zero (defaults are 0) states for the free joint
	// States are from initial ik step. States are passed via ICs
	OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
	//cout<<"\n\nCS = "<<CS;

	for (int i = 0; i<CS.getSize(); i++)//
	{
		//nb first 6 coords are prescribed (torso assumed to be fixed)
		CS[i].setValue(s,ICs[i]);
		CS[i].setSpeedValue(s,0);
		
	}

	// Print out the initial position and velocity states
	// s.getQ().dump("Initial q's"); // block positions
	// s.getU().dump("Initial u's"); // block velocities
	// std::cout << "Initial time: " << s.getTime() << std::endl;
	

	// Integrate from initial time to final time
	manager.setInitialTime(initialTime);
	manager.setFinalTime(finalTime);
	std::cout<<"\n\nIntegrating from "<<initialTime<<" to "<<finalTime<<std::endl;

	manager.integrate(s);

	cout<<"\n\nintegrating done\n\n";
	
	// return pk analysis results
	int no_markers = osimModel.getMarkerSet().getSize();

	// Declare 2D array with all pk analysis results in it. Can then directly compare with trc data
	Array<Array<double>> a(0,3*no_markers+1);
	double resample_dt = 0.01;

	// Store marker 1 simulation results
	Storage* temp = m1h.getPositionStorage();
	Array<string> col_labels("",4);// = new Array<string>();
	col_labels[0] = "time";//,"m2","m3");		
	col_labels[1] = "m_x";
	col_labels[2] = "m_y";
	col_labels[3] = "m_z";
	temp->setColumnLabels(col_labels);

	temp->getTimeColumn(a[0]);
	
	temp->getDataColumn("m_x",a[1]);
	temp->getDataColumn("m_y",a[2]);
	temp->getDataColumn("m_z",a[3]);

	//cout<<"\n\na = "<<a[0];
	//cout<<"\n\na = "<<a[1];
	//cout<<"\n\na = "<<a[2];
	//cout<<"\n\na = "<<a[3];

	// Store marker 2 simulation results
	temp = m2h.getPositionStorage();
	//temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[4]);
	temp->getDataColumn("m_y",a[5]);
	temp->getDataColumn("m_z",a[6]);

	//cout<<"\n\na = "<<a[4];
	//cout<<"\n\na = "<<a[5];
	//cout<<"\n\na = "<<a[6];

	// Store marker 3 simulation results
	temp = m3h.getPositionStorage();
	//tempm3->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[7]);
	temp->getDataColumn("m_y",a[8]);
	temp->getDataColumn("m_z",a[9]);

	//cout<<"\n\na = "<<a[7];
	//cout<<"\n\na = "<<a[8];
	//cout<<"\n\na = "<<a[9];

	// Store marker 4 simulation results
	temp = m4h.getPositionStorage();
	//tempm3->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[10]);
	temp->getDataColumn("m_y",a[11]);
	temp->getDataColumn("m_z",a[12]);

	//cout<<"\n\na = "<<a[10];
	//cout<<"\n\na = "<<a[11];
	//cout<<"\n\na = "<<a[12];

	// Store marker 5 simulation results
	temp = m5h.getPositionStorage();
	//tempm3->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[13]);
	temp->getDataColumn("m_y",a[14]);
	temp->getDataColumn("m_z",a[15]);

	//cout<<"\n\na = "<<a[13];
	//cout<<"\n\na = "<<a[14];
	//cout<<"\n\na = "<<a[15];

	// Store marker 6 simulation results
	temp = m6h.getPositionStorage();
	/*temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version6/Model/Marker2_positionv2.sto");*/
	temp->setColumnLabels(col_labels);

	temp->getDataColumn("m_x",a[16]);
	temp->getDataColumn("m_y",a[17]);
	temp->getDataColumn("m_z",a[18]);

	//cout<<"\n\na = "<<a[16];
	//cout<<"\n\na = "<<a[17];
	//cout<<"\n\na = "<<a[18];


	// Read trc storage object and insert data into 2D array for comparison with point kinematic data
	Array<Array<double>> trc_data(0,18);
	// NOTE: only comparing the markers on the head so the index is offset by 12 to account for the first
	// 4 markers in the trc file being the prescribed spine markers.
	for(int i =12;i<(12+3*6);i++)
	{
		data_trc.getDataColumn(i,trc_data[i-12]);
		//cout<<"\nDATA:    "<<trc_data[i-12]; //-> checked that data only to 5 decimal points (i.e. 100 microns) */
	}

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	// WE NOW HAVE 2 2D ARRAYS 
	//PK_DATA = DATA FROM SIMULATION (POINT KINEMATICS)
	//TRC_DATA = DATA FROM MEASURED TRC FILE.
	//NB. THEY SHOULD HAVE THE SAME NUMBER OF ARRAYS BUT ARRAYS IN EACH WILL BE OF DIFFERENT LENGTH
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //


	// Specify time values for each data set -> necessary for interploation
	Array<double> TRC_TIME, SIM_TIME; // array of time values from trc and simulation
	
	data_trc.getTimeColumn(TRC_TIME);
	SIM_TIME = a[0];

	double rms = 0;
	double rms2 = 0;
	double rms_total = 0;

	// number of markers
	//int no_markers = trc_data.getSize()/3;
	double ed = 0;

	Array<double> init(0.0,6);//);

	//cout<<"\n\ntrc data length:"<<trc_data[0].getSize();
	Array<Array<double>> ed_array(init,trc_data[0].getSize());

	// Create splines of all the simulation results so can compare with same position in trc file	
	OpenSim::GCVSplineSet *spline_set = new OpenSim::GCVSplineSet();
	GCVSpline m1x(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*0 + 0 + 1][0]);
	GCVSpline m1y(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*0 + 1 + 1][0]);
	GCVSpline m1z(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*0 + 2 + 1][0]);
	
	GCVSpline m2x(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*1 + 0 + 1][0]);
	GCVSpline m2y(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*1 + 1 + 1][0]);
	GCVSpline m2z(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*1 + 2 + 1][0]);
	
	GCVSpline m3x(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*2 + 0 + 1][0]);
	GCVSpline m3y(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*2 + 1 + 1][0]);
	GCVSpline m3z(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*2 + 2 + 1][0]);
	
	GCVSpline m4x(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*3 + 0 + 1][0]);
	GCVSpline m4y(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*3 + 1 + 1][0]);
	GCVSpline m4z(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*3 + 2 + 1][0]);
	
	GCVSpline m5x(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*4 + 0 + 1][0]);
	GCVSpline m5y(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*4 + 1 + 1][0]);
	GCVSpline m5z(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*4 + 2 + 1][0]);
	
	GCVSpline m6x(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*5 + 0 + 1][0]);
	GCVSpline m6y(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*5 + 1 + 1][0]);
	GCVSpline m6z(3, SIM_TIME.getSize(), &SIM_TIME[0], &a[3*5 + 2 + 1][0]);

	spline_set->insert(0,m1x);
	spline_set->insert(1,m1y);
	spline_set->insert(2,m1z);
	spline_set->insert(3,m2x);
	spline_set->insert(4,m2y);
	spline_set->insert(5,m2z);
	spline_set->insert(6,m3x);
	spline_set->insert(7,m3y);
	spline_set->insert(8,m3z);
	spline_set->insert(9,m4x);
	spline_set->insert(10,m4y);
	spline_set->insert(11,m4z);
	spline_set->insert(12,m5x);
	spline_set->insert(13,m5y);
	spline_set->insert(14,m5z);
	spline_set->insert(15,m6x);
	spline_set->insert(16,m6y);
	spline_set->insert(17,m6z);
	
	// create temporary spline to use for comparison
	GCVSpline * spline = new GCVSpline();

	// for each time step 
	for (int i = 0; i < trc_data[0].getSize(); i++) {

		// for each marker at each time step
		for (int k=0;k<6;k++) {
			ed = 0; // initialise eucledian distance variable

			// for each coord for marker (i.e. x, y, z)
			for (int j = 0; j<3; j++) {
				
				// get spline from spline set which contains correct marker in simulation simulation
				spline = spline_set->getGCVSpline(3*k + j);
				
				// Obtain current time step in trc file -> note: spline.calcvalue requres vector input
				SimTK::Vector inputTime(1,TRC_TIME[i]);
				
				// Caclulate the square of the eucledican distance (i.e (x-x1)^2 + (y1-y2)^2 + (z1-z2)^2 )				
				ed += (trc_data[3*k + j][i] - spline->calcValue(inputTime))*(trc_data[3*k + j][i] - spline->calcValue(inputTime));

			}
			
			// return eucliedian distance between measurements and model i.e ed = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
			// Store in opensim array structure
			//if (isNaN(ed)){
			//	cout<<"ERROR ed is negative";
			//	break;
			//}
			ed_array[i][k] = sqrt(ed);
		
		}
	}

	// Each column of ed_array consists of distance between markeri and model_markeri throughout time. 
	// I want to then get the rms of each column and store in rms_list (note: size 6 as 6 markers being included)
	Array<double> rms_list(0,6);

	// Loop through each marker column
	for (int j = 0; j<6; j++) {
		
		//cout<<"\nEUCLEDIAN ARRAY: ";	
		for (int i = 0; i<trc_data[0].getSize(); i++) {
			//cout<<ed_array[i][j]<<"\t";

			// accumatively sum the distance squared (sq.err)
			rms_list[j]+=pow(ed_array[i][j],(double)2);
		}
	}

	//cout<<"\n\nrms_list: "<<rms_list;
	// divide by length of trc file to get mean squared error and sum to return total mean squared error (rms_total)
	for (int i=0;i<rms_list.getSize();i++){
		rms_list[i] = rms_list[i]/(trc_data[0].getSize());
		rms_total += rms_list[i];
	}

	// Destruct pointers
	//spline->~GCVSpline();
	//spline_set->~GCVSplineSet();
	cout<<"\n\nrms total: "<<rms_total;
	cout<<"\nbestsoFar: "<<bestSoFar;

	//// Use an if statement to only store and print the results of an 
	////  optimization step if it is better than a previous result.
	if( rms_total < bestSoFar){
		// Save the states
		Storage statesDegrees(manager.getStateStorage());
		statesDegrees.resampleLinear(resample_dt);
		statesDegrees.print(fd);
		bestSoFar = rms_total;
		cout << "\nobjective evaluation #: " << stepCount << "  PARAMS = " << PARAMS <<  " bestSoFar = " << rms_total << std::endl;


	}	



	return rms_total;	

}

void SimTools::ApplyForce(Storage& force_storage,Model& osimModel)
{
		// USE force file to prescribe forces to head_markers body

		/////////////////////////////////////
		// NOTE: There seemed to be problem invoking "setdatacolumn" and saving result directly to double. Done a work around where save result as array, ...
		// and then convert to double for input into GCVSpline
		Array<double> time = Array<double>();
		Array<double> fx = Array<double>();
		Array<double> fy = Array<double>();
		Array<double> fz = Array<double>();
		Array<double> tx = Array<double>();
		Array<double> ty = Array<double>();
		Array<double> tz = Array<double>();

		force_storage.getTimeColumn(time,-1);
		//cout<<"\n\n"<<time;

		force_storage.getDataColumn("f1_vx",fx);
		force_storage.getDataColumn("f1_vy",fy);
		force_storage.getDataColumn("f1_vz",fz);
		force_storage.getDataColumn("t1_x",tx);
		force_storage.getDataColumn("t1_y",ty);
		force_storage.getDataColumn("t1_z",tz);
		//cout<<"\n\n"<<fx;

		

		vector<double> fx_v,fy_v,fz_v,tx_v,ty_v,tz_v;
		vector<double> time_v;

		for(int j = 0; j<force_storage.getSize(); j++)
		{
			time_v.push_back(time.get(j));
			fx_v.push_back(fx.get(j));
			fy_v.push_back(fy.get(j));
			fz_v.push_back(fz.get(j));
			tx_v.push_back(tx.get(j));
			ty_v.push_back(ty.get(j));
			tz_v.push_back(tz.get(j));
		} 

		double* t = &time_v[0];
		double* f1x = &fx_v[0];
		double* f1y = &fy_v[0];
		double* f1z = &fz_v[0];
		double* t1x = &tx_v[0];
		double* t1y = &ty_v[0];
		double* t1z = &tz_v[0];
		//cout<<"\n\n"<<*(t+10);
		//cout<<"\n\n"<<*(f1z+10);

		// Create function for all forces and torques
		OpenSim::Function* spline_fx = new OpenSim::GCVSpline(3,force_storage.getSize(),t,f1x,"fx",0);
		OpenSim::Function* spline_fy = new OpenSim::GCVSpline(3,force_storage.getSize(),t,f1y,"fy",0);
		OpenSim::Function* spline_fz = new OpenSim::GCVSpline(3,force_storage.getSize(),t,f1z,"fz",0);
		OpenSim::Function* spline_tx = new OpenSim::GCVSpline(3,force_storage.getSize(),t,t1x,"tx",0);
		OpenSim::Function* spline_ty = new OpenSim::GCVSpline(3,force_storage.getSize(),t,t1y,"ty",0);
		OpenSim::Function* spline_tz = new OpenSim::GCVSpline(3,force_storage.getSize(),t,t1z,"tz",0);
		/////////////////////////////////////

		OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
		OpenSim::Function* function1 = new OpenSim::Constant(0);

		// Create pointer to BodySet inside osimModel
		BodySet& bodyset = osimModel.updBodySet();

		// Add external prescribed force to head_markers body using spline functions created. 
		OpenSim::Body * head_markers = &bodyset.get("head_markers");
		PrescribedForce *force1 = new PrescribedForce(head_markers); 

		// Prescribe functions for specifying the point that force is applied to (function 1 = constant value of 0)
		force1->setPointFunctions(function1,function1,function1);

		// Create linear function
		// Array<double> function_coeffs(0,2);
		// function_coeffs[0] = 5;
		// function_coeffs[1] = 0;
		// OpenSim::Function* function2 = new OpenSim::LinearFunction(function_coeffs);
		// Create constant function
		// OpenSim::Function* function3 = new OpenSim::Constant(8);

		// Prescribe force and torque functions using splines from force measurement file
		force1->setForceFunctions(spline_fx,spline_fy,spline_fz);
		force1->setTorqueFunctions(spline_tx,spline_ty,spline_tz);

		// Specify that force and point are NOT in global frame (i.e. expressed in frame of head_markers body)
		force1->setForceIsInGlobalFrame(0);
		force1->setPointIsInGlobalFrame(0);
		force1->setName("External_Load");
		// add force to osim model
		osimModel.addForce(force1); 
}

void SimTools::AddPKAtest(OpenSim::Model& osimModel)
{
//		 Create list of point kinematics reporters which match marker positions from osim model
		PointKinematics* m1h = new PointKinematics(&osimModel);
		PointKinematics* m2h = new PointKinematics(&osimModel);
		PointKinematics* m3h = new PointKinematics(&osimModel);
		PointKinematics* m4h = new PointKinematics(&osimModel);
		PointKinematics* m5h = new PointKinematics(&osimModel);
		PointKinematics* m6h = new PointKinematics(&osimModel);
		
		// points in head_marker ref frame specified by inspection
		Vec3 p1(0.071700000,0.00000000,-0.054772000);
		Vec3 p2(-0.071700000,0.00000000,-0.038524000);
		Vec3 p3(0.071700000,0.00000000,0.005228);
		Vec3 p4(-0.0300000,0.00000000,0.026928);
		Vec3 p5(0.0300000,0.00000000,0.026928);
		Vec3 p6(-0.071700000,0.00000000,-0.008524000);

		m1h->setBodyPoint("head_markers",p1);
		m2h->setBodyPoint("head_markers",p2);
		m3h->setBodyPoint("head_markers",p3);
		m4h->setBodyPoint("head_markers",p4);
		m5h->setBodyPoint("head_markers",p5);
		m6h->setBodyPoint("head_markers",p6);
		
		m1h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m2h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m3h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m4h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m5h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m6h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		
		m1h->setName("m1");
		m2h->setName("m2");
		m3h->setName("m3");
		m4h->setName("m4");
		m5h->setName("m5");
		m6h->setName("m6");

		osimModel.addAnalysis(m1h);
		osimModel.addAnalysis(m2h);
		osimModel.addAnalysis(m3h);
		osimModel.addAnalysis(m4h);
		osimModel.addAnalysis(m5h);
		osimModel.addAnalysis(m6h);

		//AnalysisSet pk_set;
		//pk_set.adoptAndAppend(m1h);
		//pk_set.adoptAndAppend(m2h);
		//pk_set.adoptAndAppend(m3h);
		//pk_set.adoptAndAppend(m4h);
		//pk_set.adoptAndAppend(m5h);
		//pk_set.adoptAndAppend(m6h);
		//pk_set.setName("pk_set");
		//osimModel.addAnalysis(&pk_set[0]);

		//return(pk_set);
}

double SimTools::RunSimulation_FlexExt(OpenSim::Model& osimModel, Vector& PARAMS, double& TiFlex, double& TiExt, double& TfFlex, double& TfExt, Array<double>& ICs_flex, Array<double>& ICs_ext, const bool& save_states, string& fd, Storage& force_storage_flex, Storage& force_storage_ext, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h,Storage& data_trc_Flex,Storage& data_trc_Ext)
{
	/////////////////////////////////////////////////////////////////////	
	/////////////////////////////////////////////////////////////////////
	// Define all parameters 
	double ThetaStarF = PARAMS[0];
	double ThetaStarE = PARAMS[1];
	Vec3 k1(PARAMS[2]);
	double k2 = PARAMS[3]*SimTK::Pi/180;
	Vec3 Rdamping(10);//PARAMS[3]);

	Vec3 Tk1(0); // translational dof stiffness
	Vec3 Tdamping(0);//PARAMS[3]); // translational dof damping
	
	Vec3 p1bushing(0),p2bushing(0),o1(0),o2(0); // bushing properties
	double transition = 2.0; // transition region of limit stop

	double vert_mass = 0.018;
	double head_mass = 0.35;
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	// remove all forces
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0);// osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0); //osimModelE.updForceSet().remove(0); 
	osimModel.updForceSet().remove(0); //osimModelE.updForceSet().remove(0); 
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	// Re-apply forces
	
	// coordinte limit forces
	OpenSim::Force *limit1 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("t1Tot2_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit1); //osimModelE.addForce(limit1);
	OpenSim::Force *limit2 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c7Tot1_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit2); //osimModelE.addForce(limit2);
	OpenSim::Force *limit3 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c6Toc7_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit3); //osimModelE.addForce(limit3);
	OpenSim::Force *limit4 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c5Toc6_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit4); //osimModelE.addForce(limit4);
	OpenSim::Force *limit5 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c4Toc5_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit5); //osimModelE.addForce(limit5);
	OpenSim::Force *limit6 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c3Toc4_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit6); //osimModelE.addForce(limit6);
	OpenSim::Force *limit7 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c2Toc3_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit7); //osimModelE.addForce(limit7);
	OpenSim::Force *limit8 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("c1Toc2_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit8); //osimModelE.addForce(limit8);
	OpenSim::Force *limit9 = new OpenSim::CoordinateLimitForce(osimModel.updCoordinateSet().get("skToc1_FE").getName(),ThetaStarE,k2,-ThetaStarF,k2,0,transition);
	osimModel.addForce(limit9); //osimModelE.addForce(limit9);
	
	

	//Create new bushing forces and add to model
	o1 = Vec3(0,0,PARAMS[4]);

	OpenSim::BushingForce* t2t1 = new OpenSim::BushingForce("t2",p1bushing,o1,"t1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	t2t1->setName("t2t1_bushing");
	osimModel.addForce(t2t1);
	
	OpenSim::BushingForce* t1c7 = new OpenSim::BushingForce("t1",p1bushing,o1,"c7",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	t1c7->setName("t1c7_bushing");
	osimModel.addForce(t1c7);
	
	OpenSim::BushingForce* c7c6 = new OpenSim::BushingForce("c7",p1bushing,o1,"c6",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c7c6->setName("c7c6_bushing");
	osimModel.addForce(c7c6);

	OpenSim::BushingForce* c6c5 = new OpenSim::BushingForce("c6",p1bushing,o1,"c5",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c6c5->setName("c6c5_bushing");
	osimModel.addForce(c6c5);

	OpenSim::BushingForce* c5c4 = new OpenSim::BushingForce("c5",p1bushing,o1,"c4",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c5c4->setName("c5c4_bushing");
	osimModel.addForce(c5c4);

	OpenSim::BushingForce* c4c3 = new OpenSim::BushingForce("c4",p1bushing,o1,"c3",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c4c3->setName("c4c3_bushing");
	osimModel.addForce(c4c3);

	OpenSim::BushingForce* c3c2 = new OpenSim::BushingForce("c3",p1bushing,o1,"c2",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c3c2->setName("c3c2_bushing");
	osimModel.addForce(c3c2);

	OpenSim::BushingForce* c2c1 = new OpenSim::BushingForce("c2",p1bushing,o1,"c1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c2c1->setName("c2c1_bushing");
	osimModel.addForce(c2c1);

	OpenSim::BushingForce* c1sk = new OpenSim::BushingForce("c1",p1bushing,o1,"sk",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);
	c1sk->setName("c1sk_bushing");
	osimModel.addForce(c1sk);		
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	// mass properties
	osimModel.updBodySet().get("sk").setMass(head_mass);
	osimModel.updBodySet().get("t2").setMass(vert_mass);
	osimModel.updBodySet().get("t1").setMass(vert_mass);
	osimModel.updBodySet().get("c7").setMass(vert_mass);
	osimModel.updBodySet().get("c6").setMass(vert_mass);
	osimModel.updBodySet().get("c5").setMass(vert_mass);
	osimModel.updBodySet().get("c4").setMass(vert_mass);
	osimModel.updBodySet().get("c3").setMass(vert_mass);
	osimModel.updBodySet().get("c2").setMass(vert_mass);
	osimModel.updBodySet().get("c1").setMass(vert_mass);
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	
	
	Array<Array<double>> a_flex(0,3*6+1);
	Array<Array<double>> a_ext(0,3*6+1);
	SimTools* simtools = new SimTools();
	OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();

	// Sets torso coordinates to be set @ value specified by ik analysis (assuming torso does not move....)
	for(int i = 0; i<6; i++)
	{
		OpenSim::Function* function_ik = new OpenSim::Constant(ICs_flex[i]);
		// First 6 coordinates in cooridnate set correspond to rotx,roty,rotz,tx,ty,tz
		CS[i].setPrescribedFunction(*function_ik);
		CS[i].setDefaultValue(ICs_flex[i]);
		CS[i].setDefaultIsPrescribed(true);
		

	} 

	
	simtools->ApplyForce(force_storage_flex,osimModel);
	osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/Version9_constrainedwFORCES_FLEX.osim");
	string fd_flex = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/flex.sto";	
	a_flex = DoIntegration(fd_flex, osimModel, ICs_flex,TiFlex,TfFlex,m1h,m2h,m3h,m4h,m5h,m6h);

	for(int i = 0; i<6; i++)
	{
		OpenSim::Function* function_ik = new OpenSim::Constant(ICs_ext[i]);
		// First 6 coordinates in cooridnate set correspond to rotx,roty,rotz,tx,ty,tz
		CS[i].setPrescribedFunction(*function_ik);
		CS[i].setDefaultValue(ICs_ext[i]);
		CS[i].setDefaultIsPrescribed(true);
		

	} 

	//cout<<"\nb = "<<b;

	osimModel.updForceSet().remove(osimModel.updForceSet().getSize()-1);
	string fd_ext = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/ext.sto";	
	simtools->ApplyForce(force_storage_ext,osimModel);
	osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/Version9_constrainedwFORCES_ext.osim");
	 
	a_ext = DoIntegration(fd_ext, osimModel, ICs_ext,TiExt,TfExt,m1h,m2h,m3h,m4h,m5h,m6h);
//	cout<<"\nb = "<<b;
	//cout<<"\nflex size: "<<a_flex.getSize();
	//cout<<"\next size: "<<a_ext.getSize();

	a_ext.setSize(3*6 +1);
	a_flex.setSize(3*6+1);

	//cout<<"\nflex size: "<<a_flex[0];
	//cout<<"\n\next size: "<<a_ext[0];

	//cout<<"\n\nflex size: "<<a_flex[1];
	//cout<<"\n\next size: "<<a_ext[1];

	//// Remove flexion pk analysis... unsure why includeD?
	//for (int i = 0; i<a_flex.getSize(); i++){
	//	for (int j = 0; j<a_flex[0].getSize(); j++){	
	//		a_ext[i].remove(0);
	//	}
	//}



	// -> 19 array of length n representing t,x,y,z... for 6 markers for both flex and ext
	// want to interpolate this to same dt as trc.

////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	// Read trc storage object and insert data into 2D array for comparison with point kinematic data
	Array<Array<double>> trc_data_flex(0,18),trc_data_ext(0,18);

	// NOTE: only comparing the markers on the head so the index is offset by 12 to account for the first
	// 4 markers in the trc file being the prescribed spine markers.
	for(int i =12;i<(12+3*6);i++)
	{
		data_trc_Flex.getDataColumn(i,trc_data_flex[i-12]);
		data_trc_Ext.getDataColumn(i,trc_data_ext[i-12]);

		//trc_data_flex[i-12].append(trc_data_ext[i-12].getSize(),&trc_data_ext[i-12][0]);
		//trc_data_flex.append(trc_data_ext[0].getSize(),trc_data_ext);
		//cout<<"\nTRC DATA:    "<<trc_data_flex[i-12]; //-> checked that data only to 5 decimal points (i.e. 100 microns) */
		//cout<<"\nDATA:    "<<trc_data_ext[i-12];
		//cout<<"\nDATA zerp:    "<<trc_data_ext[i-12][0];
	}

	//cout<<"\ntrc FLEX: "<<trc_data_flex<<endl;
	//cout<<"\ntrc EXT: "<<trc_data_ext<<endl;
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	// WE NOW HAVE 2 2D ARRAYS 
	//PK_DATA = DATA FROM SIMULATION (POINT KINEMATICS)
	//TRC_DATA = DATA FROM MEASURED TRC FILE.
	//NB. THEY SHOULD HAVE THE SAME NUMBER OF ARRAYS BUT ARRAYS IN EACH WILL BE OF DIFFERENT LENGTH
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //


	//// Specify time values for each data set -> necessary for interploation
	Array<double> TRC_TIME_flex, TRC_TIME_ext, SIM_TIME_flex, SIM_TIME_ext; // array of time values from trc and simulation
	//
	data_trc_Flex.getTimeColumn(TRC_TIME_flex);
	data_trc_Ext.getTimeColumn(TRC_TIME_ext);

	//cout<<"\ntrc FLEX time"<<TRC_TIME_flex<<endl;
	//cout<<"\ntrc EXT time: "<<TRC_TIME_ext<<endl;


	//double* test = &TRC_TIME_ext[0];
	//TRC_TIME_flex.append(TRC_TIME_ext.getSize(),test);
	//cout<<"\nTRC TIME FLEX: "<<TRC_TIME_flex<<endl;
	////cout<<"TRC TIME EXT: "<<TRC_TIME_ext<<endl;

	SIM_TIME_ext = a_ext[0];
	SIM_TIME_flex = a_flex[0];
	cout<<"\nSIM TIME EXT: "<<SIM_TIME_ext<<endl;
	cout<<"\nSIM TIME flex: "<<SIM_TIME_flex<<endl;	

	double rms = 0;
	double rms2 = 0;
	double rms_total = 0;

	//// number of markers
	int no_markers = trc_data_flex.getSize()/3;
	double ed = 0;

	// initialise eucledian distance array
	Array<double> init(0.0,no_markers);
	// Create array of arrays -> initialised to 6 zeros (1 for each marker) for each time instance in trc file
	Array<Array<double>> ed_array(init,trc_data_flex[0].getSize());

	// create spline set for both flexion data and extension data...
	// Create splines of all the simulation results so can compare with same position in trc file	
	OpenSim::GCVSplineSet *Extspline_set = new OpenSim::GCVSplineSet();
	OpenSim::GCVSplineSet *Flexspline_set = new OpenSim::GCVSplineSet();
	GCVSpline m1xE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*0 + 0 + 1][0]); GCVSpline m1xF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*0 + 0 + 1][0]); 
	GCVSpline m1yE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*0 + 1 + 1][0]); GCVSpline m1yF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*0 + 1 + 1][0]); 
	GCVSpline m1zE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*0 + 2 + 1][0]); GCVSpline m1zF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*0 + 2 + 1][0]); 
	
	GCVSpline m2xE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*1 + 0 + 1][0]); GCVSpline m2xF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*1 + 0 + 1][0]); 
	GCVSpline m2yE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*1 + 1 + 1][0]); GCVSpline m2yF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*1 + 1 + 1][0]); 
	GCVSpline m2zE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*1 + 2 + 1][0]); GCVSpline m2zF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*1 + 2 + 1][0]); 
	
	GCVSpline m3xE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*2 + 0 + 1][0]); GCVSpline m3xF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*2 + 0 + 1][0]); 
	GCVSpline m3yE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*2 + 1 + 1][0]); GCVSpline m3yF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*2 + 1 + 1][0]); 
	GCVSpline m3zE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*2 + 2 + 1][0]); GCVSpline m3zF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*2 + 2 + 1][0]); 

	GCVSpline m4xE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*3 + 0 + 1][0]); GCVSpline m4xF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*3 + 0 + 1][0]); 
	GCVSpline m4yE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*3 + 1 + 1][0]); GCVSpline m4yF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*3 + 1 + 1][0]); 
	GCVSpline m4zE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*3 + 2 + 1][0]); GCVSpline m4zF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*3 + 2 + 1][0]); 

	GCVSpline m5xE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*4 + 0 + 1][0]); GCVSpline m5xF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*4 + 0 + 1][0]); 
	GCVSpline m5yE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*4 + 1 + 1][0]); GCVSpline m5yF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*4 + 1 + 1][0]); 
	GCVSpline m5zE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*4 + 2 + 1][0]); GCVSpline m5zF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*4 + 2 + 1][0]); 

	GCVSpline m6xE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*5 + 0 + 1][0]); GCVSpline m6xF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*5 + 0 + 1][0]); 
	GCVSpline m6yE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*5 + 1 + 1][0]); GCVSpline m6yF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*5 + 1 + 1][0]); 
	GCVSpline m6zE(3, SIM_TIME_ext.getSize(), &SIM_TIME_ext[0], &a_ext[3*5 + 2 + 1][0]); GCVSpline m6zF(3, SIM_TIME_flex.getSize(), &SIM_TIME_flex[0], &a_flex[3*5 + 2 + 1][0]); 

	// Insert splines into spline sets for easier access (indexing...)
	Extspline_set->insert(0,m1xE); Flexspline_set->insert(0,m1xF);
	Extspline_set->insert(1,m1yE); Flexspline_set->insert(1,m1yF);
	Extspline_set->insert(2,m1zE); Flexspline_set->insert(2,m1zF);
	Extspline_set->insert(3,m2xE); Flexspline_set->insert(3,m2xF);
	Extspline_set->insert(4,m2yE); Flexspline_set->insert(4,m2yF);
	Extspline_set->insert(5,m2zE); Flexspline_set->insert(5,m2zF);
	Extspline_set->insert(6,m3xE); Flexspline_set->insert(6,m3xF);
	Extspline_set->insert(7,m3yE); Flexspline_set->insert(7,m3yF);
	Extspline_set->insert(8,m3zE); Flexspline_set->insert(8,m3zF);
	Extspline_set->insert(9,m4xE); Flexspline_set->insert(9,m4xF);
	Extspline_set->insert(10,m4yE); Flexspline_set->insert(10,m4yF);
	Extspline_set->insert(11,m4zE); Flexspline_set->insert(11,m4zF);
	Extspline_set->insert(12,m5xE); Flexspline_set->insert(12,m5xF);
	Extspline_set->insert(13,m5yE); Flexspline_set->insert(13,m5yF);
	Extspline_set->insert(14,m5zE); Flexspline_set->insert(14,m5zF);
	Extspline_set->insert(15,m6xE); Flexspline_set->insert(15,m6xF);
	Extspline_set->insert(16,m6yE); Flexspline_set->insert(16,m6yF);
	Extspline_set->insert(17,m6zE); Flexspline_set->insert(17,m6zF);
	
	// create temporary spline to use for comparison
	GCVSpline * spline = new GCVSpline();

	// for each time step 
	for (int i = 0; i < trc_data_flex[0].getSize(); i++) {

		// for each marker at each time step
		for (int k=0;k<6;k++) {
			ed = 0; // initialise eucledian distance variable

			// for each coord for marker (i.e. x, y, z)
			for (int j = 0; j<3; j++) {
				
				// get spline from spline set which contains correct marker in simulation simulation
				spline = Flexspline_set->getGCVSpline(3*k + j);
				
				// Obtain current time step in trc file -> note: spline.calcvalue requres vector input
				SimTK::Vector inputTime(1,TRC_TIME_flex[i]);
				
				// Caclulate the square of the eucledican distance (i.e (x-x1)^2 + (y1-y2)^2 + (z1-z2)^2 )
				////cout<<"\ninput time: "<<spline->;
				//cout<<"\ntrc time: "<<trc_data_flex[3*k + j][i];
				//cout<<"\nsimulation value: "<<spline->calcValue(inputTime);
				ed += (trc_data_flex[3*k + j][i] - spline->calcValue(inputTime))*(trc_data_flex[3*k + j][i] - spline->calcValue(inputTime));

			}
			
			// return eucliedian distance between measurements and model i.e ed = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
			// Store in opensim array structure
			//if (isNaN(ed)){
			//	cout<<"ERROR ed is negative";
			//	break;
			//}
			ed_array[i][k] = sqrt(ed);
		
		}
	//cout<<"\ned array flex: "<<ed_array[i]<<endl;
	}

	Array<Array<double>> ed_array_flex = ed_array;


	// for each time step 
	for (int i = 0; i < trc_data_ext[0].getSize(); i++) {

		// for each marker at each time step
		for (int k=0;k<6;k++) {
			ed = 0; // initialise eucledian distance variable

			// for each coord for marker (i.e. x, y, z)
			for (int j = 0; j<3; j++) {
				
				// get spline from spline set which contains correct marker in simulation simulation
				spline = Extspline_set->getGCVSpline(3*k + j);
				
				// Obtain current time step in trc file -> note: spline.calcvalue requres vector input
				SimTK::Vector inputTime(1,TRC_TIME_ext[i]);
				
				// Caclulate the square of the eucledican distance (i.e (x-x1)^2 + (y1-y2)^2 + (z1-z2)^2 )
				////cout<<"\ninput time: "<<spline->;
				//cout<<"\ntrc time: "<<trc_data_flex[3*k + j][i];
				//cout<<"\nsimulation value: "<<spline->calcValue(inputTime);
				ed += (trc_data_ext[3*k + j][i] - spline->calcValue(inputTime))*(trc_data_ext[3*k + j][i] - spline->calcValue(inputTime));

			}
			
			// return eucliedian distance between measurements and model i.e ed = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
			// Store in opensim array structure
			//if (isNaN(ed)){
			//	cout<<"ERROR ed is negative";
			//	break;
			//}
			ed_array[i][k] = sqrt(ed);
		
		}
		//cout<<"\ned array ext : "<<ed_array[i]<<endl;
	}

	Array<Array<double>> ed_array_ext = ed_array;

	
	//ed_array_flex.append(ed_array_ext.gets
	//cout<<"\ned array: "<<ed_array_ext.getSize()<<endl;
	for (int i = 0 ; i<ed_array_ext.getSize(); i++){
		ed_array_flex[i].append(ed_array_ext[0].getSize(),&ed_array_ext[i][0]);
	}

	//cout<<"\nAppended array: "<<ed_array_flex<<endl;

	// Each column of ed_array consists of distance between markeri and model_markeri throughout time. 
	// I want to then get the rms of each column and store in rms_list (note: size 12 as 12 markers being included -> 6 from each simulaiton)
	Array<double> rms_list(0,12);

	// Loop through each marker column
	for (int j = 0; j<12; j++) {
		
		//cout<<"\nEUCLEDIAN ARRAY: ";	
		for (int i = 0; i<trc_data_flex[0].getSize(); i++) {
			//cout<<ed_array[i][j]<<"\t";

			// accumatively sum the distance squared (sq.err)
			rms_list[j]+=pow(ed_array_flex[i][j],(double)2);
		}
	}

	//cout<<"\n rms list: "<<rms_list<<endl;
	
	//cout<<"\n\nrms_list: "<<rms_list;
	// divide by length of trc file to get mean squared error and sum to return total mean squared error (rms_total)
	for (int i=0;i<rms_list.getSize();i++){
		rms_list[i] = rms_list[i]/(trc_data_flex[0].getSize()+trc_data_ext[0].getSize());
		rms_total += rms_list[i];
	}

	// Destruct pointers
	//spline->~GCVSpline();
	//spline_set->~GCVSplineSet();
	cout<<"\n\nrms total: "<<rms_total;
	cout<<"\nbestsoFar: "<<bestSoFar;

	double resample_dt = 0.1;
	//// Use an if statement to only store and print the results of an 
	////  optimization step if it is better than a previous result.
	if( rms_total < bestSoFar){
		// Save the states
		//Storage statesDegrees(manager.getStateStorage());
		//statesDegrees.resampleLinear(resample_dt);
		//statesDegrees.print(fd);
		bestSoFar = rms_total;

		cout << "\nobjective evaluation #: " << stepCount << "  PARAMS = " << PARAMS <<  " bestSoFar = " << rms_total << std::endl;

	}
	//////	//////////////////////////////
	////////////////////////////////////////
	//////
	
	return(rms_total);

}



Array<Array<double>> SimTools::DoIntegration(string& fd, Model& osimModel,Array<double>& ICs, double& initialTime, double& finalTime, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h)
{
	
	
	// Initialise and run simulation
	SimTK::State& s = osimModel.initSystem();
	osimModel.updMultibodySystem().updDefaultState();

	
	//m2h.deleteStorage();
	//m3h.deleteStorage();
	//m4h.deleteStorage();
	//m5h.deleteStorage();
	//m6h.deleteStorage();

	// Create the integrator, force reporter, and manager for the simulation.
	SimTK::RungeKuttaMersonIntegrator integrator(osimModel.getMultibodySystem());
	// default settings
	integrator.setAccuracy(0.00001);
	integrator.setMaximumStepSize(1);
	integrator.setMinimumStepSize(1e-8);
	
	// Create the force reporter
	// ForceReporter* reporter = new ForceReporter(&osimModel);
	// osimModel.addAnalysis(reporter);

	// Create the manager
	Manager manager(osimModel,  integrator);

	// Print out details of the model
	//osimModel.printDetailedInfo(s, std::cout);

	// Define non-zero (defaults are 0) states for the free joint
	// States are from initial ik step. States are passed via ICs
	OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
	//cout<<"\n\nCS = "<<CS;

	for (int i = 0; i<CS.getSize(); i++)//
	{
		//nb first 6 coords are prescribed (torso assumed to be fixed)
		CS[i].setValue(s,ICs[i]);
		CS[i].setSpeedValue(s,0);
		
	}

	// Print out the initial position and velocity states
	// s.getQ().dump("Initial q's"); // block positions
	// s.getU().dump("Initial u's"); // block velocities
	// std::cout << "Initial time: " << s.getTime() << std::endl;
	

	// Integrate from initial time to final time
	manager.setInitialTime(initialTime);
	manager.setFinalTime(finalTime);
	std::cout<<"\n\nIntegrating from "<<initialTime<<" to "<<finalTime<<std::endl;

	manager.integrate(s);

	cout<<"\n\nintegrating done\n\n";
	
	// return pk analysis results
	int no_markers = osimModel.getMarkerSet().getSize();

	// Declare 2D array with all pk analysis results in it. Can then directly compare with trc data
	Array<Array<double>> a(0,3*no_markers+1);
	double resample_dt = 0.01;

	// print results
	Storage statesDegrees(manager.getStateStorage());
	statesDegrees.resampleLinear(resample_dt);
	statesDegrees.print(fd);	


	
	//OpenSim::PointKinematics* m1h = dynamic_cast<PointKinematics*>(&as[0]);
	// Store marker 1 simulation results
	Storage* temp = m1h.getPositionStorage();
	temp->crop(initialTime,finalTime);

	Array<string> col_labels("",4);// = new Array<string>();
	col_labels[0] = "time";//,"m2","m3");		
	col_labels[1] = "m_x";
	col_labels[2] = "m_y";
	col_labels[3] = "m_z";
	temp->setColumnLabels(col_labels);
	
	temp->getTimeColumn(a[0]);
	temp->getDataColumn("m_x",a[1]);
	temp->getDataColumn("m_y",a[2]);
	temp->getDataColumn("m_z",a[3]);

	//cout<<"\n\na = "<<a[0];
	//cout<<"\n\na = "<<a[1];
	//cout<<"\n\na = "<<a[2];
	//cout<<"\n\na = "<<a[3];

	//OpenSim::PointKinematics* m2h = dynamic_cast<PointKinematics*>(&as[1]);
	// Store marker 2 simulation results
	temp = m2h.getPositionStorage();
	temp->crop(initialTime,finalTime);
	//temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);
	temp->getDataColumn("m_x",a[4]);
	temp->getDataColumn("m_y",a[5]);
	temp->getDataColumn("m_z",a[6]);	
	
	//OpenSim::PointKinematics* m3h = dynamic_cast<PointKinematics*>(&as[2]);
	// Store marker 2 simulation results
	temp = m3h.getPositionStorage();
	temp->crop(initialTime,finalTime);
	//temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);
	temp->getDataColumn("m_x",a[7]);
	temp->getDataColumn("m_y",a[8]);
	temp->getDataColumn("m_z",a[9]);	

	//OpenSim::PointKinematics* m4h = dynamic_cast<PointKinematics*>(&as[3]);
	// Store marker 2 simulation results
	temp = m4h.getPositionStorage();
	temp->crop(initialTime,finalTime);
	//temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);
	temp->getDataColumn("m_x",a[10]);
	temp->getDataColumn("m_y",a[11]);
	temp->getDataColumn("m_z",a[12]);	

	//OpenSim::PointKinematics* m5h = dynamic_cast<PointKinematics*>(&as[4]);
	// Store marker 2 simulation results
	temp = m5h.getPositionStorage();
	temp->crop(initialTime,finalTime);
	//temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);
	temp->getDataColumn("m_x",a[13]);
	temp->getDataColumn("m_y",a[14]);
	temp->getDataColumn("m_z",a[15]);	

	//OpenSim::PointKinematics* m6h = dynamic_cast<PointKinematics*>(&as[5]);
	// Store marker 2 simulation results
	temp = m6h.getPositionStorage();
	temp->crop(initialTime,finalTime);
	//temp->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version2/Marker2_positionv2.sto");
	temp->setColumnLabels(col_labels);
	temp->getDataColumn("m_x",a[16]);
	temp->getDataColumn("m_y",a[17]);
	temp->getDataColumn("m_z",a[18]);	

	return(a);
}

