#include "GaussianCF.h"

using namespace std;
using namespace itpp;

GaussianCF::GaussianCF(double lengthscale, double var) : CovarianceFunction("Isotropic Gaussian")
{
	numberParameters = 2;
	setDefaultTransforms();

	range = lengthscale;
	variance = var;
}

GaussianCF::GaussianCF(vec parameters) : CovarianceFunction("Isotropic Gaussian")
{
	numberParameters = 2;
	assert(parameters.size() == getNumberParameters());
	range = parameters(0);
	variance = parameters(1);
	setDefaultTransforms();
}


GaussianCF::~GaussianCF()
{
}

inline double GaussianCF::computeElement(const vec& A, const vec& B) const
{
	return calcGaussian(A - B);
}

inline double GaussianCF::computeDiagonalElement(const vec& A) const
{
	return calcGaussianDiag();
}

inline double GaussianCF::calcGaussian(const vec& V) const
{
    return variance * exp( -0.5 * sum(sqr(V)) / sqr(range) );
}

inline double GaussianCF::calcGaussianDiag() const
{
	return variance;
}

void GaussianCF::setParameter(int parameterNumber, const double value)
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	switch(parameterNumber)
	{
		case 0 : range = value;
					break;
		case 1 : variance = value;
					break;
		default: assert(false);
					break;
	}
}

double GaussianCF::getParameter(int parameterNumber) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	switch(parameterNumber)
	{
		case 0 : return(range);
					break;
		case 1 : return(variance);
					break;
		default: assert(false);
					break;
	}
	cerr << "Warning: should not have reached here in GaussianCF::getParameter" << endl;
	return(0.0);
}

string GaussianCF::getParameterName(int parameterNumber) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	switch(parameterNumber)
	{
		case 0 : return("Range");
					break;
		case 1 : return("Variance");
					break;
		default: break;

	}
	return("Unknown parameter");
}

void GaussianCF::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	Transform* t = getTransform(parameterNumber);
	double gradientModifier = t->gradientTransform(getParameter(parameterNumber));

	switch(parameterNumber)
	{
		case 0 :
		{
			mat DM;
			DM.set_size(PD.rows(), PD.cols());
			computeSymmetric(PD, X);
			computeDistanceMatrix(DM, X);
			elem_mult_inplace(DM  * (gradientModifier / pow(range, 3.0)), PD);
			return;
			break;
		}

		case 1 :
		{
			computeSymmetric(PD, X);
			PD *= (gradientModifier / variance);
			return;
			break;
		}
	}
	cerr << "Warning: should not have reached here in GaussianCF::getParameterPartialDerivative" << endl;
}
