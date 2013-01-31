#include "Model_setup.h"
#include "CoupledBushingForceEDIT.h"
#include "Force_Plugin\FunctionBasedBushingForce.h"
/* ################################################## */
// Intial set up functions
Model_setup::Model_setup(){}

void Model_setup::addExtForce(OpenSim::Model& osimModel, string& force_filename)


{
		//*** SET UP EXTERNAL FORCE IN MODEL ***//
		
		// Create storage object for force file
		OpenSim::Storage* force_storage = new Storage(force_filename);

		// smooths out all data....
		force_storage->smoothSpline(3,1); // order = 3, cutoff freq = 1Hz
		force_storage->resampleLinear(0.1); // resample to 10Hz
		

		// USE force file to prescribe forces to head_markers body
		// NOTE: There seemed to be problem invoking "getdatacolumn" and saving result directly to double. Done a work around where save result as array, ...
		// and then convert to double for input into GCVSpline
		Array<double> time = Array<double>();
		Array<double> fx = Array<double>();
		Array<double> fy = Array<double>();
		Array<double> fz = Array<double>();
		Array<double> tx = Array<double>();
		Array<double> ty = Array<double>();
		Array<double> tz = Array<double>();

		force_storage->getTimeColumn(time,-1);
		force_storage->getDataColumn("f1_vx",fx);
		force_storage->getDataColumn("f1_vy",fy);
		force_storage->getDataColumn("f1_vz",fz);
		force_storage->getDataColumn("t1_x",tx);
		force_storage->getDataColumn("t1_y",ty);
		force_storage->getDataColumn("t1_z",tz);

		vector<double> fx_v,fy_v,fz_v,tx_v,ty_v,tz_v;
		vector<double> time_v;

		for(int j = 0; j<force_storage->getSize(); j++)
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

		// Create function for all forces and torques
		OpenSim::Function* spline_fx = new OpenSim::GCVSpline(3,force_storage->getSize(),t,f1x,"fx",0);
		OpenSim::Function* spline_fy = new OpenSim::GCVSpline(3,force_storage->getSize(),t,f1y,"fy",0);
		OpenSim::Function* spline_fz = new OpenSim::GCVSpline(3,force_storage->getSize(),t,f1z,"fz",0);
		OpenSim::Function* spline_tx = new OpenSim::GCVSpline(3,force_storage->getSize(),t,t1x,"tx",0);
		OpenSim::Function* spline_ty = new OpenSim::GCVSpline(3,force_storage->getSize(),t,t1y,"ty",0);
		OpenSim::Function* spline_tz = new OpenSim::GCVSpline(3,force_storage->getSize(),t,t1z,"tz",0);
		/////////////////////////////////////

		OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
		BodySet& bodyset = osimModel.updBodySet();

		// Zero function -> define position where force acts
		OpenSim::Function* Zero_func = new OpenSim::Constant(0);

		// Add external prescribed force to head_markers body using spline functions created. 
		OpenSim::Body * head_markers = &bodyset.get("head_markers");
		PrescribedForce *force1 = new PrescribedForce(head_markers); 

		// Prescribe functions for specifying the point that force is applied to (function 1 = constant value of 0)
		force1->setPointFunctions(Zero_func,Zero_func,Zero_func);

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

void Model_setup::addInitialCoupledBushing(OpenSim::Model& osimModel)
{		
		OpenSim::Body& lower = osimModel.updBodySet().get("t2");
		OpenSim::Body& upper = osimModel.updBodySet().get("sk");
		Vec3 R_damp(10);
		Vec3 Tr_damp(10);// correct units? be multiplied by 180/pi?
		
		//CB = coupled bushing position (p) i and orientation (o) i
		Vec3 CBp1(0),CBo1(0),CBp2(0),CBo2(0);
		
		//CoupledBushingForceEDIT* F = new CoupledBushingForceEDIT(lower.getName(),CBp1,CBo1,upper.getName(),CBp2,CBo2,K,D);
		//F->setName("CoupledBushing1");
		//osimModel.addForce(F);

		OpenSim::FunctionBasedBushingForce* F = new OpenSim::FunctionBasedBushingForce(lower.getName(),CBp1,CBo1,upper.getName(),CBp2,CBo2);
		// set "diagonal terms"
		F->setName("CoupledBushing1");
		F->set_rotational_damping(R_damp);
		F->set_translational_damping(Tr_damp);
		osimModel.addForce(F);
}

void Model_setup::addPKAnalyses(OpenSim::Model& osimModel)
{
		// Set up point kinematic analyses -> to track marker positions for use in objective function

		//cout<<"\n\nINITAL ANALYSIS SET: "<<osimModel.getAnalysisSet()<<"\n\n";

		// Create list of point kinematics reporters which match marker positions from osim model
		PointKinematics* m1h = new PointKinematics(&osimModel);
		PointKinematics* m2h = new PointKinematics(&osimModel);
		PointKinematics* m3h = new PointKinematics(&osimModel);
		PointKinematics* m4h = new PointKinematics(&osimModel);
		PointKinematics* m5h = new PointKinematics(&osimModel);
		PointKinematics* m6h = new PointKinematics(&osimModel);

		// Set to return coordinates of point (specified wrt head_markers) but expressed in ground ref frame
		//m1h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m2h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m3h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m4h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m5h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m6h->setBody(&osimModel.updBodySet().get("head_markers"));

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
}