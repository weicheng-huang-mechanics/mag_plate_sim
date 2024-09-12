#include "magneticForce.h"
#include <iostream>

magneticForce::magneticForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;

  currentIndex = 0;
}

magneticForce::~magneticForce()
{
	;
}

void magneticForce::computeFm(double m_time)
{
  if (currentIndex < plate->timeSeries.size() - 1 && plate->timeSeries[currentIndex + 1] <= m_time)
  {
    currentIndex++;
  }

  Vector3d curt_ba = plate->ba_vec[currentIndex];

  ba1 = curt_ba(0);
  ba2 = curt_ba(1);
  ba3 = curt_ba(2);

	for (int i = 0; i < plate->triangularNum; i++)
	{
		ind1 = plate->v_triangularElement[i].nv_1;
		ind2 = plate->v_triangularElement[i].nv_2;
		ind3 = plate->v_triangularElement[i].nv_3;

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

		v1 = plate->getVertexOld(ind1);
		v2 = plate->getVertexOld(ind2);
		v3 = plate->getVertexOld(ind3);

   // cout << i << endl;

    //cout << plate->br_vec.size() << endl;

    Vector3d localBr = plate->br_vec[i];

    br1 = localBr(0);
    br2 = localBr(1);
    br3 = localBr(2);

    /*

    Vector3d vc = plate->v_triangularElement[i].v_c;

    //cout << vc.norm() << endl;

    if (vc.norm() > 0.1)
    {
      vc = vc / vc.norm();

      br1 = brRef * vc(0);
      br2 = brRef * vc(1);
      br3 = brRef * vc(2);
    }
    else
    {
      br1 = 0.0;
      br2 = 0.0;
      br3 = 0.0;
    }

    cout << i << " " << br1 << " " << br2 << " " << br3 << endl;

    */

		x1 = v1(0);
		y1 = v1(1);
		z1 = v1(2);
		x2 = v2(0);
		y2 = v2(1);
		z2 = v2(2);
		x3 = v3(0);
		y3 = v3(1);
		z3 = v3(2);

		Matrix3d aInv = plate->v_triangularElement[i].abarinv;

		m11 = aInv(0, 0);
		m12 = aInv(0, 1);
		m13 = aInv(0, 2);

		m21 = aInv(1, 0);
		m22 = aInv(1, 1);
		m23 = aInv(1, 2);

		m31 = aInv(2, 0);
		m32 = aInv(2, 1);
		m33 = aInv(2, 2);

		forceVec = - plate->v_triangularElement[i].area * plate->thickness * computeForce();

    //cout << plate->v_triangularElement[i].area << " " << plate->thickness << endl;

    //cout << "ba " << ba1 << " " << ba2 << " " << ba3 << endl;
    //cout << "br " << br1 << " " << br2 << " " << br3 << endl;
    //cout << "force " << computeForce().norm() << endl;

		for (int j = 0; j < 9; j++)
        {
            stepper->addForce(arrayNum(j), - forceVec(j));
        }

	}
}

VectorXd magneticForce::computeForce()
{
  VectorXd vecResult;

  vecResult = ListVec(ba1*(br1*(m21 - (m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(m22 - (m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(m23 - (m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba3*(br1*(-(m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(y2 - y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(y2 - y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(y2 - y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba2*(br1*(-(m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-z2 + z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-z2 + z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(y2 - y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z2 + z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-z2 + z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))),
   ba2*(br1*(m21 - (m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(m22 - (m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(m23 - (m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba3*(br1*(-(m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-x2 + x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-x2 + x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-x2 + x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba1*(br1*(-(m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(z2 - z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(z2 - z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x2 + x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z2 - z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(z2 - z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))),
   ba3*(br1*(m21 - (m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(m22 - (m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(m23 - (m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba2*(br1*(-(m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(x2 - x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(x2 - x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(x2 - x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba1*(br1*(-(m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-y2 + y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-y2 + y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x2 - x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y2 + y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-y2 + y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))),
   ba1*(br1*(-m11 - (m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(-m12 - (m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(-m13 - (m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba3*(br1*(-(m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-y1 + y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-y1 + y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-y1 + y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba2*(br1*(-(m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(z1 - z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(z1 - z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(-y1 + y3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(z1 - z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))),
   ba2*(br1*(-m11 - (m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(-m12 - (m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(-m13 - (m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba3*(br1*(-(m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(x1 - x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(x1 - x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(x1 - x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba1*(br1*(-(m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-z1 + z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-z1 + z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x1 - x3)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-z1 + z3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))),
   ba3*(br1*(-m11 - (m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(-m12 - (m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(-m13 - (m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba2*(br1*(-(m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-x1 + x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-x1 + x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-x1 + x3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba1*(br1*(-(m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(y1 - y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(y1 - y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x1 + x3)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(y1 - y3)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(y1 - y3))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))),
   ba1*(br1*(m11 - m21 - (m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(m12 - m22 - (m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(m13 - m23 - (m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
             (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba3*(br1*(-(m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(y1 - y2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(y1 - y2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(y1 - y2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba2*(br1*(-(m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-z1 + z2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-z1 + z2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(y1 - y2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-z1 + z2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))),
   ba2*(br1*(m11 - m21 - (m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(m12 - m22 - (m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(m13 - m23 - (m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
             (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
               2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba3*(br1*(-(m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-x1 + x2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-x1 + x2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
              (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-x1 + x2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba1*(br1*(-(m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(z1 - z2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(z1 - z2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(-x1 + x2)*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3) + 
                2*(z1 - z2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(z1 - z2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))),
   ba3*(br1*(m11 - m21 - (m31*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br2*(m12 - m22 - (m32*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5))) + 
       br3*(m13 - m23 - (m33*(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)*
             (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
               2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)))) + 
    ba2*(br1*(-(m31*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(x1 - x2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(x1 - x2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)*
              (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(x1 - x2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))) + 
    ba1*(br1*(-(m31*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m31*(-y1 + y2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br2*(-(m32*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m32*(-y1 + y2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2))) + 
       br3*(-(m33*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)*
              (2*(x1 - x2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + 
                2*(-y1 + y2)*(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3)))/
           (2.*pow(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
               pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
               pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2),1.5)) + 
          (m33*(-y1 + y2))/
           sqrt(pow(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3,2) + 
             pow(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3,2) + 
             pow(-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3,2)))));

  return vecResult;
}

VectorXd magneticForce::ListVec(double a1, double a2, double a3, double a4, double a5, 
	double a6, double a7, double a8, double a9)
{
  VectorXd vecResult;

  vecResult.setZero(9, 1);

  vecResult(0) = a1;
  vecResult(1) = a2;
  vecResult(2) = a3;
  vecResult(3) = a4;
  vecResult(4) = a5;
  vecResult(5) = a6;
  vecResult(6) = a7;
  vecResult(7) = a8;
  vecResult(8) = a9;

  return vecResult;
}