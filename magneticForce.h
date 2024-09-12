#ifndef MAGNETICFORCE_H
#define MAGNETICFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class magneticForce
{
public:
	magneticForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~magneticForce();
	void computeFm(double m_time);
    
private:
	elasticPlate *plate;
    timeStepper *stepper;

    int currentIndex;

    Vector3d v1, v2, v3;
    int ind1, ind2, ind3;

    VectorXi arrayNum;
    VectorXd forceVec;

    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3; 

    double brRef;
    double ba1, ba2, ba3;
    double br1, br2, br3;

    double m11, m12, m13, m21, m22, m23, m31, m32, m33;

    VectorXd computeForce();
    VectorXd ListVec(double a1, double a2, double a3, double a4, double a5, 
    double a6, double a7, double a8, double a9);

};

#endif
