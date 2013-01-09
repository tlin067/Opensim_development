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

	osimModel.invalidateSystem();

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
	Vec3 Tk1(0); // translational dof stiffness
	Vec3 Tdamping(0);//PARAMS[3]); // translational dof damping
	
	//// Define position in body 1 (p1), position in body 2 (p2), orientation in body 1 (o1), orientation in  body 2 (o2)
	//// These are position and orientation os joints wrt bodies. Can set to be constant throughout...
	Vec3 p1bushing(0),p2bushing(0),o1(0),o2(0);
	double transition = 2.0;

	double vert_mass = 0.018;
	double head_mass = 0.35;//PARAMS[4];
	
	o1 = Vec3(0,0,PARAMS[3]); // bushing offset

	OpenSim::BushingForce* t2t1 = new OpenSim::BushingForce("t2",p1bushing,o1,"t1",p2bushing,o2,Tk1,k1,Tdamping,Rdamping);



	// obtain pointer to model forces
	OpenSim::ForceSet& force_set = osimModel.updForceSet();

	// change limit force properties
	for (int i = 0; i<1; i++){
		
		//((OpenSim::CoordinateLimitForce*)&force_set.get(i))->setDamping(0);
		((OpenSim::CoordinateLimitForce*)&force_set.get(i))->setUpperLimit(ThetaStar);
		((OpenSim::CoordinateLimitForce*)&force_set.get(i))->setLowerLimit(-ThetaStar);
		((OpenSim::CoordinateLimitForce*)&force_set.get(i))->setUpperStiffness(k2);
		((OpenSim::CoordinateLimitForce*)&force_set.get(i))->setLowerStiffness(k2);
		//((OpenSim::CoordinateLimitForce*)&force_set.get(i))->setTransition(transition);

	}

	// change linear stiffness of bushing
	for (int i = 1; i<2; i++){

		((OpenSim::BushingForce*)&force_set.get(i))->set_rotational_stiffness(k1);
		((OpenSim::BushingForce*)&force_set.get(i))->set_orientation_body_1(o1); // chnage bushing offset
	}


	// set mass properties
	osimModel.updBodySet().get("sk").setMass(head_mass);
	osimModel.updBodySet().get("t2").setMass(vert_mass);


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/OpenSim new versions (from Version10)/OpenSim Output/OSIM_constrainedwFORCES_FLEX.osim");
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


	// DELETE POINTERS???
	

	return rms_total;	

}