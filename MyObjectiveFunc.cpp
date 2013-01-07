#include "MyObjectiveFunc.h"

// ny = number parameters
// newGuess = current guess of optimizer system
// myOptSys = optimizer system to calculate objective function for
MyObjectiveFunc::MyObjectiveFunc(int ny, const Vector &newGuess, const MyOptimizerSystem *myOptSys):
Differentiator::GradientFunction(ny),time(0),newGuess(newGuess), myOptSys(myOptSys){};

Real MyObjectiveFunc::Call_obj(Vector GUESS) const{
	
	// Call obj function and save result to f
	Real f = NULL;
	myOptSys->objectiveFunc(GUESS,true,f);

	return f;
};

int MyObjectiveFunc::f(const Vector &newGuess, Real& fx) const {
    
	////Evaluate objective function of optimizer system and save scalar result to fx	
    //myOptSys->objectiveFunc(newGuess,true,fx);//
	fx = Call_obj(newGuess);	
	cout<<"\n\nX = "<<fx;

	//myOptSys->objectiveFunc(newGuess,true,fx);
	return 0; // success
};

// Already defined functions...
Real MyObjectiveFunc::getTime() const {return time;}

void MyObjectiveFunc::setTime(Real t) {time=t;}



//************************************//
//************************************//
SinOmegaX::SinOmegaX(const Vector &newGuess, const MyOptimizerSystem *myOptSys) : newGuess(newGuess), myOptSys(myOptSys){};

Real SinOmegaX::Call_obj() const{
	
	Real f;
	myOptSys->objectiveFunc(newGuess,true,f);

	return f;
};


int SinOmegaX::f(Real x, Real& fx) const {
    
	
	//fx = newGuess[0];//std::sin(w*x);
    fx = Call_obj();	
	
	cout<<"\n\nX = "<<fx;
	return 0; // success
};