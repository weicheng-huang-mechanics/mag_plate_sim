#include "dampingForce.h"

dampingForce::dampingForce(elasticPlate &m_plate, timeStepper &m_stepper, double m_viscosity)
{
	plate = &m_plate;
    stepper = &m_stepper;

	viscosity = m_viscosity;

	TotalForceVec = VectorXd::Zero(plate->ndof);
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{
	
	TotalForceVec = VectorXd::Zero(plate->ndof);

	for (int i = 0; i < plate->nv; i++)
	{
		u = plate->getVelocity(i);

		f = - u * plate->massArray(3 * i) * viscosity;

		for (int k = 0; k < 3; k++)
		{
			ind = 3 * i + k;

			//TotalForceVec(ind) = TotalForceVec(ind) + f[k];
			
			stepper->addForce(ind, - f[k]);
		}
	}
	

	/*

	for (int i = 0; i < plate->triangularNum; i++)
	{
		int ind1 = plate->v_triangularElement[i].nv_1;
		int ind2 = plate->v_triangularElement[i].nv_2;
		int ind3 = plate->v_triangularElement[i].nv_3;

		Vector3d x1 = plate->getVertexOld(ind1);
		Vector3d x2 = plate->getVertexOld(ind2);
		Vector3d x3 = plate->getVertexOld(ind3);

		Vector3d v1 = plate->getVelocity(ind1);
		Vector3d v2 = plate->getVelocity(ind2);
		Vector3d v3 = plate->getVelocity(ind3);

		VectorXi arrayNum;

		arrayNum = VectorXi::Zero(9);

		arrayNum(0) = 3 * ind1 + 0;
		arrayNum(1) = 3 * ind1 + 1;
		arrayNum(2) = 3 * ind1 + 2;

		arrayNum(3) = 3 * ind2 + 0;
		arrayNum(4) = 3 * ind2 + 1;
		arrayNum(5) = 3 * ind2 + 2;

		arrayNum(6) = 3 * ind3 + 0;
		arrayNum(7) = 3 * ind3 + 1;
		arrayNum(8) = 3 * ind3 + 2;


		Vector3d v_c = (v1 + v2 + v3) / 3;

		Vector3d a1 = x3 - x2;
		Vector3d a2 = x1 - x3;
		Vector3d surfN = a1.cross(a2) / ( (a1.cross(a2)).norm() );

		Vector3d vN = v_c.dot(surfN) * surfN;

		Vector3d force = - 0.5 * viscosity * vN.norm() * vN * plate->v_triangularElement[i].area;

		for (int j = 0; j < 3; j++)
        {
            stepper->addForce(arrayNum(j), - force(j) / 3);
        }
        for (int j = 0; j < 3; j++)
        {
            stepper->addForce(arrayNum(j+3), - force(j) / 3);
        }
        for (int j = 0; j < 3; j++)
        {
            stepper->addForce(arrayNum(j+6), - force(j) / 3);
        }




	}

	*/


	/*
	for (int i=0;i<rod->ne;i++)
	{
		u = rod->getVelocity(i);

		f = - u * viscosity * rod->voronoiLen(i);
		
		for (int k = 0; k < 3; k++)
		{
			ind = 4*i + k;
			stepper->addForce(ind, - f[k]);
		}
	}
	*/
}

void dampingForce::computeJd()
{
	
	for (int i = 0; i < plate->nv; i++)
	{
		u = plate->getVelocity(i);

		jac(0) = - viscosity * plate->massArray(3 * i) / plate->dt;
		jac(1) = - viscosity * plate->massArray(3 * i) / plate->dt;
		jac(2) = - viscosity * plate->massArray(3 * i) / plate->dt;
		
		for (int j = 0; j < 3; j++)
		{
			ind = 3 * i + j;
			stepper->addJacobian(ind, ind, - jac(j));
		}
	}
	
}

void dampingForce::setFirstJacobian()
{
	for (int i = 0; i < plate->nv; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ind = 3 * i + j;
			stepper->addJacobian(ind, ind, 1);
		}
	}
}
