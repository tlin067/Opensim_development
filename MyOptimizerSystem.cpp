
#include "MyObjectiveFunc.h"
#include "my_classes.h"

extern double bestSoFar;
extern int stepCount;

//MyOptimizerSystem::MyOptimizerSystem(int& numParams, Model& osimModel, Storage& data_trc, double& ti, double& tf, Array<double>& ICs, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h, string& output_fd):
//numParams(numParams), OptimizerSystem(numParams), osimModel(osimModel), data_trc(data_trc), ti(ti), tf(tf), ICs(ICs), m1h(m1h), m2h(m2h), m3h(m3h), m4h(m4h), m5h(m5h), m6h(m6h), LIMIT_FORCE(LIMIT_FORCE), output_fd(output_fd)

MyOptimizerSystem::MyOptimizerSystem(int& numParams, OpenSim::Model& osimModel, double& TiFlex, double& TiExt, double& TfFlex, double& TfExt, Array<double>& ICs_flex, Array<double>& ICs_ext, Storage& force_storage_flex, Storage& force_storage_ext, OpenSim::PointKinematics& m1h, OpenSim::PointKinematics& m2h, OpenSim::PointKinematics& m3h, OpenSim::PointKinematics& m4h, OpenSim::PointKinematics& m5h, OpenSim::PointKinematics& m6h,Storage& data_trc_Flex,Storage& data_trc_Ext):
numParams(numParams), OptimizerSystem(numParams), osimModel(osimModel), TiFlex(TiFlex),	TiExt(TiExt), TfFlex(TfFlex), TfExt(TfExt),ICs_flex(ICs_flex), ICs_ext(ICs_ext), force_storage_flex(force_storage_flex), force_storage_ext(force_storage_ext), m1h(m1h), m2h(m2h), m3h(m3h), m4h(m4h), m5h(m5h), m6h(m6h), data_trc_Flex(data_trc_Flex), data_trc_Ext(data_trc_Ext)
{
	setNumEqualityConstraints( 0 );
	setNumInequalityConstraints( 0 );
}


int MyOptimizerSystem::objectiveFunc(  const Vector &newGuess, const bool new_coefficients, Real& f ) const 
{

	// Run simulation and save point kinematic reporter in data_sim storage object. Print simulation results to file.
	SimTools* stools = new SimTools();
	//Storage* temp;
	Array<Array<double>> a;
	//string fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version8/Version8_output/flex_output.sto";

	cout<<"\n";
	cout<<"OBJ FUN GUESS: "<<newGuess;
	cout<<"\n";

	Vector PARAMS(numParams);
	for (int i = 0; i<numParams; i++){
		PARAMS[i] = newGuess[i];
	}
	
	double t1 = clock();
	string output_fd = "";
	//a = stools->RunSimulation_LIMITSTOP(osimModel,PARAMS,ti,tf,ICs,false,fd,m1h,m2h,m3h,m4h,m5h,m6h);
	//f = stools->RunSimulation_wRMS(data_trc,osimModel,PARAMS,ti,tf,ICs,false,output_fd,m1h,m2h,m3h,m4h,m5h,m6h);
	f = stools->RunSimulation_FlexExt(osimModel,PARAMS,TiFlex,TiExt,TfFlex,TfExt,ICs_flex,ICs_ext,false,output_fd,force_storage_flex,force_storage_ext,m1h,m2h,m3h,m4h,m5h,m6h,data_trc_Flex,data_trc_Ext);
	double t2 = clock();
	cout<<"\nobjective function sim time: "<<((float)t2-(float)t1)/ CLOCKS_PER_SEC << " seconds";

	//double rms = stools->Calc_rms_VERSION2(data_trc,ti,tf,a);

	delete stools;

	//cout<<"\n\nRMS: "<< rms; 
	//double rms_scaled = (rms-0.02)/(0.22-0.02);
	//cout<<"\n\nRMS scaled: "<< rms_scaled; 
	//f = rms;//_scaled;

	stepCount++;

	return(0);
}

int MyOptimizerSystem::gradientFunc( const Vector &newGuess, bool new_coefficients, Vector &gradient ) const
{
   
	//Real out;
	//MyOptimizerSystem::f(0.0,out);
	//cout<<"\n\nOUT: "<<out;
	//cout<<"\n\nGRAD GUESS: "<<newGuess;
	//SinOmegaX test(newGuess,this);
	cout<<"\n\nGRAIDENT new guess: "<<newGuess;
	MyObjectiveFunc test(numParams,newGuess,this);
	

	Differentiator diff(test);
	
	//cout<<"\n\nmethod order: "<<diff.getMethodOrder(Differentiator::);
	gradient = diff.calcGradient(newGuess,Differentiator::ForwardDifference);//(newGuess,out,gradient);
	
	//cout<<"\n\n NEW GUESS after GRAD FUNCTION: "<<newGuess;
	
	
	cout<<"\n\nGRAIDENT: "<<gradient;

     return(0);

 }
