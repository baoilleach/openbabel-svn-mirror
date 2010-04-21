/**********************************************************************
qtpie.cpp - A OBChargeModel to handle QTPIE charges

Copyright (C) 2010 by Jiahao Chen <jiahao@mit.edu>
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifdef HAVE_EIGEN2

#include "qtpie.h"

using namespace std;

namespace OpenBabel
{
/////////////////////////////////////////////////////////////////
QTPIECharges theQTPIECharges("qtpie"); //Global instance

/////////////////////////////////////////////////////////////////


  Vector3d QTPIECharges::GetParameters(unsigned int Z, int Q)
  {
    //Returns a triple of numbers: electronegativity (in eV), hardness (in eV), and Gaussian exponent (in bohr^-2)

    Vector3d P;
    //For now, completely ignore the formal charge
    switch(Z) {
	case  1: P<<4.528*eV,	13.890*eV,	0.534337523756312; break;
	case  3: P<<3.006*eV,	 4.772*eV,	0.166838519142176; break;
	case  6: P<<5.343*eV,	10.126*eV,	0.206883838259186; break;
	case  7: P<<7.139*eV,	12.844*eV,	0.221439796025873; break;
	case  8: P<<8.741*eV,	13.364*eV,	0.223967308625516; break;
	case  9: P<<10.87*eV,	14.948*eV,	0.231257590182828; break;
	case 11: P<<2.843*eV,	 4.592*eV,	0.095892938712585; break;
	case 14: P<<4.168*eV,	 6.974*eV,	0.105219608142377; break;
	case 15: P<<5.463*eV,	 8.000*eV,	0.108476721661715; break;
	case 16: P<<6.084*eV,	10.660*eV,	0.115618357843499; break;
	case 17: P<<8.564*eV,	 9.892*eV,	0.113714050615107; break;
	case 19: P<<2.421*eV,	 3.840*eV,	0.060223294377778; break;
	case 35: P<<7.790*eV,	 8.850*eV,	0.070087547802259; break;
	case 37: P<<2.331*eV,	 3.692*eV,	0.041999054745368; break;
	case 53: P<<6.822*eV,	 7.524*eV,	0.068562697575073; break;
	case 55: P<<2.183*eV,	 3.422*eV,	0.030719481189777; break;
	default: P<<0., 1.e10, 1.e10; //Magic value
    }
    return P;
  }

  bool QTPIECharges::ComputeCharges(OBMol &mol)
  {

     ///////////////////////////////////////////////////////////////////////////////
    //Some OpenBabel bookkeeping that I copied from the Gasteiger scheme
    mol.SetPartialChargesPerceived();

    OBPairData *dp = new OBPairData;
    dp->SetAttribute("PartialCharges");
    dp->SetValue("QTPIE");
    dp->SetOrigin(perceived);
    mol.SetData(dp);


    ///////////////////////////////////////////////////////////////////////////////
    //Read in atomic information from OpenBabel molecule and parameterize

    //Read in total number of atoms
    int i, N = mol.NumAtoms();

    Hardness = MatrixXd::Zero(N+1, N+1);
    Voltage = VectorXd::Zero(N+1);
    Electronegativity = VectorXd::Zero(N);
    VectorXd BasisSet = VectorXd::Zero(N);

    Vector3d Parameters;

    FOR_ATOMS_OF_MOL(atom, mol)
    {
       	Parameters = GetParameters(atom->GetAtomicNum(), atom->GetFormalCharge());
	i = atom->GetIdx() - 1;

	if (Parameters[0] == 0.)
        {
		stringstream msg;
		msg << "Some QTPIE Parameters not found!" << endl
	            << "Parameters not found for atom no. " << i+1 << endl
		    << "Atom will be ignored in the charge computation.";
		obErrorLog.ThrowError(__FUNCTION__, msg.str(), obError);
	}

	Electronegativity[i] = Parameters[0];
	Hardness(i,i) = Parameters[1];
	BasisSet[i] = Parameters[2];
    }

    //Read in total charge of molecule
    double Z = (double)(mol.GetTotalCharge());

    if (Z != 0.) obErrorLog.ThrowError(__FUNCTION__, "Warning, total charge on molecule is not zero. QTPIE routine may give nonsense.", obWarning);
    

    ///////////////////////////////////////////////////////////////////////////////
    // Now populate integrals


    // Calculate integral pre-screening threshold
    double SmallestGaussianExponent = BasisSet.minCoeff();
    double CoulMaxDistance = 2 * sqrt(-log(CoulombThreshold) / SmallestGaussianExponent);

    int j;
    double R, Coulomb;
    FOR_ATOMS_OF_MOL(atom1, mol)
    {
	i = atom1->GetIdx() - 1;
        FOR_ATOMS_OF_MOL(atom2, mol)
	{
	    j = atom2->GetIdx() - 1;
	    if (i>j)
	    {
		//For some reason, this did _not_ produce the expected pairwise distances
		//(2,1) evaluated to the same distance as (2,0) in H2O: bug?
		// - jiahao@mit.edu 2010-04-20 r3691
		//R = atom1->GetDistance(atom2)*Angstrom;
		//
	        R = atom1->GetDistance(j+1)*Angstrom;
	        if (R<CoulMaxDistance)
		    Coulomb = CoulombInt(BasisSet[i], BasisSet[j], R);
	        else
		    Coulomb = 1./R;
	        Hardness(i,j) = Coulomb;
	        Hardness(j,i) = Coulomb;
	    }
	}
    }

    Hardness.block(N,0,1,N).setOnes();
    Hardness.block(0,N,N,1).setOnes();

    double OverlapMaxDistance = sqrt( log( (pi/cube(2*SmallestGaussianExponent))
			      / sqr(OverlapThreshold)) /SmallestGaussianExponent);

    //This matrix can be sparse, but I didn't get Eigen's SparseMatrix to
    //play well with this - jiahao@mit.edu 2010-04-20 eigen-2.0.12 r3691
    MatrixXd Overlap = MatrixXd::Zero(N,N);
    double OverlapVal;
    FOR_ATOMS_OF_MOL(atom1, mol)
    {
	i = atom1->GetIdx() - 1;
        FOR_ATOMS_OF_MOL(atom2, mol)
	{
	    j = atom2->GetIdx() - 1;
	    if (i>j)
	    { 
	        R = atom1->GetDistance(j+1)*Angstrom;
	        if (R<OverlapMaxDistance)
	        {
		    OverlapVal = OverlapInt(BasisSet[i], BasisSet[j], R);
		    Overlap(i,j) = OverlapVal;
		    Overlap(j,i) = OverlapVal;
	        }
	    }
	}
    }
    
    // Calculate normalization factors
    VectorXd OvNorm(N);
    for (i=0; i<N; i++) OvNorm[i] = 1.0 / (1.0 + Overlap.row(i).sum());

    // Calculate voltages
    double PotentialDiff, Norm, ThisVoltage, ThisOverlap;
    for (i=0; i<N; i++)
    {
	Norm = OvNorm[i];
	ThisVoltage = 0.;
	for (j=0; j<N; j++)
	{
		PotentialDiff = Electronegativity[i] - Electronegativity[j];
		ThisOverlap = Overlap(i,j);
		if (ThisOverlap > OverlapThreshold)
		       ThisVoltage -= PotentialDiff * Norm* ThisOverlap;
	}
    	
	Voltage[i] = ThisVoltage;
    }
    
    Voltage[N] = Z;

    ///////////////////////////////////////////////////////////////////////////////
    // Call linear algebra solver
    //
    bool status = solver(Hardness, Voltage, Charge);
    if (!status)
    {
	obErrorLog.ThrowError(__FUNCTION__, "Charges could not be computed.", obError);
	return false;
    }

    ChemicalPotential = Charge[N];

//For debugging purposes only
#if 0
    cout << "Ready to solve QTPIE problem" << endl << endl;
    cout << "Hardness matrix:" << endl << Hardness << endl << endl;
    cout << "Voltage vector:" << endl << Voltage << endl << endl;
    cout << "Overlap matrix:" << endl << Overlap << endl << endl;
    cout << "Overlap norm:" << endl << OvNorm << endl << endl;
    cout << "Charge vector:" << endl << Charge << endl << endl;
#endif

    //Now we are done calculating, pass all this back to OpenBabel molecule
    m_partialCharges.clear();
    m_partialCharges.reserve(mol.NumAtoms());
    m_formalCharges.clear();
    m_formalCharges.reserve(mol.NumAtoms());
    FOR_ATOMS_OF_MOL(atom, mol)
    {
        atom->SetPartialCharge(Charge[atom->GetIdx()-1]);
        m_partialCharges.push_back(Charge[atom->GetIdx() - 1]);
        m_formalCharges.push_back(atom->GetFormalCharge());
    }
    obErrorLog.ThrowError(__FUNCTION__, "QTPIE charges successfully assigned.", obInfo);
    return true;
  }


/// Calculates Coulomb integral
double QTPIECharges::CoulombInt(double a, double b, double R)
{
	double p = sqrt(a * b / (a + b));
	return erf(p * R) / R;
}

double QTPIECharges::OverlapInt(double a, double b, double R)
{
	double p = a + b;
	double q = a * b / p;
	return pow(4*q/p, 0.75) * exp(-q*R*R);
}

/// Here's a wrapper around the Eigen solver routine
bool QTPIECharges::solver(MatrixXd A, VectorXd b, VectorXd &x, const double NormThreshold)
{
	// using a LU factorization
	bool SolverOK = A.lu().solve(b, &x);
	//bool SolverOK = A.svd().solve(b, &x);
	
	VectorXd resid = A*x - b;
	double resnorm = resid.norm();
	if (isnan(resnorm) || resnorm > NormThreshold || !SolverOK)
	{
		stringstream msg;
		msg << "Warning, LU solver failed." << endl;
		if (!SolverOK) msg << "Solver returned error." << endl;
		if (isnan(resnorm)) msg << "NaNs were returned" << endl;
		if (resnorm > NormThreshold) msg << "Residual has norm " << resnorm
		       << " which exceeds the recommended threshold of " << NormThreshold
		       << endl;
		msg << "Proceeding with singular value decomposition.";

		obErrorLog.ThrowError(__FUNCTION__, msg.str(), obWarning);

		SolverOK = A.svd().solve(b, &x);
		resid = A*x - b;
		resnorm = resid.norm();
		
		if (isnan(resnorm) || !SolverOK)
		{
			obErrorLog.ThrowError(__FUNCTION__, "SVD solver returned an error. Charges may not be reliable!", obError);
			return false;
		}
	}

	stringstream msg_resid;
	msg_resid << "The residual of the solution has norm " << resnorm;
	obErrorLog.ThrowError(__FUNCTION__, msg_resid.str(), obInfo);

	if (resnorm > NormThreshold) {
		stringstream msg_reswarn;
		msg_reswarn << "Warning, the norm of the residual is " << resnorm
			    << "which exceeds the recommended threshold of " << NormThreshold;
		obErrorLog.ThrowError(__FUNCTION__, msg_reswarn.str(), obWarning);
	}
	return true;
}

}//namespace

#endif //HAVE_EIGEN2

//! \file qtpie.cpp
//! \brief Assign QTPIE partial charges.
