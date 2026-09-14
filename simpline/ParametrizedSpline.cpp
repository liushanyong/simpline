#include "Simpline.h"
#include "Constants.h"
#include <numeric>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

// taken from https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template<typename T>
std::vector<size_t> sortIndices(const std::vector<T>& unsortedVector)
{
	std::vector<size_t> sortedIndices(unsortedVector.size());
	std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
	std::sort(sortedIndices.begin(), sortedIndices.end(), [&unsortedVector](size_t i1, size_t i2){ return unsortedVector[i1] < unsortedVector[i2]; });
	return sortedIndices;
}

template<typename T>
simpline<T>::ParametrizedSpline::ParametrizedSpline():
		parameterValues(), points(), firstDerivatives(), secondDerivatives(), thirdDerivatives(), finiteDifferences(), closed(false)
{
}

template<typename T>
simpline<T>::ParametrizedSpline::ParametrizedSpline(const std::vector<T>& parameterValues, const std::vector<simpline<T>::Vector3>& points, const bool& closed):
		parameterValues(), points(), firstDerivatives(), secondDerivatives(), thirdDerivatives(), finiteDifferences(), closed(closed)
{
	if(points.size() < 2)
	{
		throw std::runtime_error("Point list must contain at least two items!");
	}

	if(parameterValues.size() != points.size())
	{
		throw std::runtime_error("Number of parameter values must be equal to number of points!");
	}

	std::vector<size_t> sortedIndices = sortIndices(parameterValues);
	std::vector<T> sortedParameterValues(parameterValues.size());
	std::vector<simpline<T>::Vector3> sortedPoints(points.size());
	for(size_t i = 0; i < sortedIndices.size(); i++)
	{
		sortedParameterValues[i] = parameterValues[sortedIndices[i]];
		sortedPoints[i] = points[sortedIndices[i]];

		if(i > 0 && sortedParameterValues[i - 1] == sortedParameterValues[i])
		{
			throw std::runtime_error("Multiple points cannot have the same parameter value.");
		}
	}

	if(this->closed && sortedPoints.size() > 2 && (sortedPoints.front() - sortedPoints.back()).norm() <= epsilon)
	{
		sortedPoints.pop_back();
		sortedParameterValues.pop_back();
	}

	if(this->closed && sortedPoints.size() < 3)
	{
		throw std::runtime_error("Closed spline requires at least three distinct points.");
	}

	this->parameterValues = sortedParameterValues;
	this->points = sortedPoints;

	if(!this->closed)
	{
		this->firstDerivatives.resize(this->points.size() - 1);
		this->secondDerivatives.resize(this->points.size());
		this->thirdDerivatives.resize(this->points.size() - 1);

		std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> A(3);
		for(size_t i = 0; i < 3; i++)
		{
			A[i] = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(this->points.size(), this->points.size());
			for(size_t j = 0; j < this->points.size() - 2; j++)
			{
				for(size_t k = 0; k < this->points.size(); k++)
				{
					T coefficient = 0.0;
					if((j + 1) == k)
					{
						coefficient = 2.0;
					}
					else if((j + 1) - k == -1)
					{
						coefficient = (this->parameterValues[j + 2] - this->parameterValues[j + 1]) / (this->parameterValues[j + 2] - this->parameterValues[j]);
					}
					else if((j + 1) - k == 1)
					{
						coefficient = (this->parameterValues[j + 1] - this->parameterValues[j]) / (this->parameterValues[j + 2] - this->parameterValues[j]);
					}

					A[i](j + 1, k) = coefficient;
				}
			}

			// second derivative at first point is known
			A[i](0, 0) = 1;
			for(size_t j = 1; j < this->points.size(); j++)
			{
				A[i](0, j) = 0;
			}

			// second derivative at last point is known
			A[i](this->points.size() - 1, this->points.size() - 1) = 1;
			for(size_t j = 0; j < this->points.size() - 1; j++)
			{
				A[i](this->points.size() - 1, j) = 0;
			}
		}

		std::vector<simpline<T>::VectorX> b(3);
		for(size_t i = 0; i < 3; i++)
		{
			b[i] = simpline<T>::VectorX(this->points.size());
			for(size_t j = 0; j < this->points.size() - 2; j++)
			{
				b[i](j + 1) = 6 * computeFiniteDifference(j, j + 2)[i];
			}

			// second derivative at first point is zero
			b[i](0) = 0;
			// second derivative at last point is zero
			b[i](this->points.size() - 1) = 0;
		}

		for(size_t i = 0; i < 3; i++)
		{
			simpline<T>::VectorX x(A[i].colPivHouseholderQr().solve(b[i]));
			for(size_t j = 0; j < this->points.size(); j++)
			{
				secondDerivatives[j](i) = x(j);
			}
		}

		for(size_t j = 0; j < this->points.size() - 1; j++)
		{
			firstDerivatives[j] = computeFiniteDifference(j, j + 1) -
							  ((this->parameterValues[j + 1] - this->parameterValues[j]) * secondDerivatives[j] / 3.0) -
							  ((this->parameterValues[j + 1] - this->parameterValues[j]) * secondDerivatives[j + 1] / 6.0);
			thirdDerivatives[j] = (secondDerivatives[j + 1] - secondDerivatives[j]) / (this->parameterValues[j + 1] - this->parameterValues[j]);
		}
		return;
	}

	const T closingInterval = (this->parameterValues.back() - this->parameterValues.front()) / static_cast<T>(this->parameterValues.size() - 1);
	if(closingInterval <= 0)
	{
		throw std::runtime_error("Failed to infer positive closing interval for closed spline.");
	}

	this->parameterValues.push_back(this->parameterValues.back() + closingInterval);
	this->points.push_back(this->points.front());

	const size_t uniquePointCount = this->points.size() - 1;
	this->firstDerivatives.resize(uniquePointCount);
	this->secondDerivatives.resize(uniquePointCount + 1);
	this->thirdDerivatives.resize(uniquePointCount);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(uniquePointCount, uniquePointCount);
	A.setZero();
	std::vector<simpline<T>::VectorX> b(3);
	for(size_t i = 0; i < 3; i++)
	{
		b[i] = simpline<T>::VectorX(uniquePointCount);
	}

	for(size_t i = 0; i < uniquePointCount; i++)
	{
		const size_t previousIndex = (i + uniquePointCount - 1) % uniquePointCount;
		const size_t nextIndex = (i + 1) % uniquePointCount;

		const T previousInterval = (i == 0) ?
				this->parameterValues[uniquePointCount] - this->parameterValues[uniquePointCount - 1] :
				this->parameterValues[i] - this->parameterValues[i - 1];
		const T nextInterval = this->parameterValues[i + 1] - this->parameterValues[i];

		if(previousInterval <= 0 || nextInterval <= 0)
		{
			throw std::runtime_error("Parameter values must be strictly increasing for closed spline.");
		}

		A(i, previousIndex) += previousInterval;
		A(i, i) = 2.0 * (previousInterval + nextInterval);
		A(i, nextIndex) += nextInterval;

		const simpline<T>::Vector3 rhs = 6.0 * (((this->points[i + 1] - this->points[i]) / nextInterval) -
													 ((this->points[i] - this->points[previousIndex]) / previousInterval));
		for(size_t d = 0; d < 3; d++)
		{
			b[d](i) = rhs(d);
		}
	}

	for(size_t i = 0; i < 3; i++)
	{
		simpline<T>::VectorX x(A.colPivHouseholderQr().solve(b[i]));
		for(size_t j = 0; j < uniquePointCount; j++)
		{
			secondDerivatives[j](i) = x(j);
		}
		secondDerivatives[uniquePointCount](i) = x(0);
	}

	for(size_t j = 0; j < uniquePointCount; j++)
	{
		const T interval = this->parameterValues[j + 1] - this->parameterValues[j];
		firstDerivatives[j] = ((this->points[j + 1] - this->points[j]) / interval) -
							  ((interval * secondDerivatives[j]) / 3.0) -
							  ((interval * secondDerivatives[j + 1]) / 6.0);
		thirdDerivatives[j] = (secondDerivatives[j + 1] - secondDerivatives[j]) / interval;
	}
}

template<typename T>
T simpline<T>::ParametrizedSpline::wrapParameterValue(const T& parameterValue) const
{
	if(!closed)
	{
		return parameterValue;
	}

	const T start = parameterValues.front();
	const T end = parameterValues.back();
	const T period = end - start;
	if(period <= 0)
	{
		return parameterValue;
	}

	T wrapped = std::fmod(parameterValue - start, period);
	if(wrapped < 0)
	{
		wrapped += period;
	}

	return start + wrapped;
}

template<typename T>
size_t simpline<T>::ParametrizedSpline::getSegmentIndex(const T& parameterValue) const
{
	size_t previousPointIndex = 0;
	while(previousPointIndex < parameterValues.size() - 1 && parameterValues[previousPointIndex + 1] <= parameterValue)
	{
		previousPointIndex++;
	}

	if(previousPointIndex >= firstDerivatives.size())
	{
		previousPointIndex = firstDerivatives.size() - 1;
	}

	return previousPointIndex;
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::computeFiniteDifference(const size_t& startPointIndex, const size_t& endPointIndex)
{
	std::pair<size_t, size_t> indexPair(startPointIndex, endPointIndex);
	if(finiteDifferences.find(indexPair) != finiteDifferences.end())
	{
		return finiteDifferences[indexPair];
	}

	simpline<T>::Vector3 finiteDifference;
	if(startPointIndex == endPointIndex)
	{
		finiteDifference = points[startPointIndex];
	}
	else
	{
		finiteDifference =
				(computeFiniteDifference(startPointIndex + 1, endPointIndex) - computeFiniteDifference(startPointIndex, endPointIndex - 1)) /
				(parameterValues[endPointIndex] - parameterValues[startPointIndex]);
	}

	finiteDifferences[indexPair] = finiteDifference;

	return finiteDifference;
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::getValue(const T& parameterValue) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get value from empty parametrized spline. Use non-default constructor to provide points.");
	}

	T evaluatedParameterValue = parameterValue;
	if(closed)
	{
		evaluatedParameterValue = wrapParameterValue(parameterValue);
	}
	else if(parameterValue < parameterValues[0] || parameterValue > parameterValues[parameterValues.size() - 1])
	{
		throw std::runtime_error("Value requested at parameterValue=" + std::to_string(parameterValue) + ". Parameter Value must be between " +
							 std::to_string(parameterValues[0]) + " and " + std::to_string(parameterValues[parameterValues.size() - 1]) + ".");
	}

	const size_t previousPointIndex = getSegmentIndex(evaluatedParameterValue);
	return points[previousPointIndex] +
		   firstDerivatives[previousPointIndex] * (evaluatedParameterValue - parameterValues[previousPointIndex]) +
		   (secondDerivatives[previousPointIndex] / 2.0) * std::pow((evaluatedParameterValue - parameterValues[previousPointIndex]), 2) +
		   (thirdDerivatives[previousPointIndex] / 6.0) * std::pow((evaluatedParameterValue - parameterValues[previousPointIndex]), 3);
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::getGradient(const T& parameterValue) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get gradient from empty parametrized spline. Use non-default constructor to provide points.");
	}

	T evaluatedParameterValue = parameterValue;
	if(closed)
	{
		evaluatedParameterValue = wrapParameterValue(parameterValue);
	}
	else if(parameterValue < parameterValues[0] || parameterValue > parameterValues[parameterValues.size() - 1])
	{
		throw std::runtime_error("Gradient requested at parameterValue=" + std::to_string(parameterValue) + ". Parameter Value must be between " +
							 std::to_string(parameterValues[0]) + " and " + std::to_string(parameterValues[parameterValues.size() - 1]) + ".");
	}

	const size_t previousPointIndex = getSegmentIndex(evaluatedParameterValue);
	return firstDerivatives[previousPointIndex] +
		   secondDerivatives[previousPointIndex] * (evaluatedParameterValue - parameterValues[previousPointIndex]) +
		   (thirdDerivatives[previousPointIndex] / 2.0) * std::pow((evaluatedParameterValue - parameterValues[previousPointIndex]), 2);
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::getSecondDerivative(const T& parameterValue) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get second derivative from empty parametrized spline. Use non-default constructor to provide points.");
	}

	T evaluatedParameterValue = parameterValue;
	if(closed)
	{
		evaluatedParameterValue = wrapParameterValue(parameterValue);
	}
	else if(parameterValue < parameterValues[0] || parameterValue > parameterValues[parameterValues.size() - 1])
	{
		throw std::runtime_error("Second derivative requested at parameterValue=" + std::to_string(parameterValue) + ". Parameter Value must be between " +
							 std::to_string(parameterValues[0]) + " and " + std::to_string(parameterValues[parameterValues.size() - 1]) + ".");
	}

	const size_t previousPointIndex = getSegmentIndex(evaluatedParameterValue);
	return secondDerivatives[previousPointIndex] +
		   thirdDerivatives[previousPointIndex] * (evaluatedParameterValue - parameterValues[previousPointIndex]);
}

template<typename T>
T simpline<T>::ParametrizedSpline::getLengthInInterval(const T& startParameterValue, const T& endParameterValue) const
{
	if(endParameterValue < startParameterValue)
	{
		throw std::runtime_error("End parameter value must be greater or equal to start parameter value.");
	}

	T length = 0.0;
	for(size_t i = 0; i < parameterValues.size() - 1; i++)
	{
		const T lowerBound = std::max(startParameterValue, parameterValues[i]);
		const T upperBound = std::min(endParameterValue, parameterValues[i + 1]);
		if(upperBound <= lowerBound)
		{
			continue;
		}

		const T intervalLength = upperBound - lowerBound;
		for(size_t j = 0; j < gaussianQuadratureAbcissa.size(); j++)
		{
			const T t = lowerBound + (((gaussianQuadratureAbcissa[j] + 1.0) / 2.0) * intervalLength);
			length += (intervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[j];
		}
	}

	return length;
}

template<typename T>
T simpline<T>::ParametrizedSpline::getLength() const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get length of empty parametrized spline. Use non-default constructor to provide points.");
	}

	return getLengthInInterval(parameterValues.front(), parameterValues.back());
}

template<typename T>
T simpline<T>::ParametrizedSpline::getLength(const T& startParameterValue, const T& endParameterValue) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get length of empty parametrized spline. Use non-default constructor to provide points.");
	}

	if(!closed)
	{
		if(startParameterValue > endParameterValue)
		{
			throw std::runtime_error("Length requested from parameterValue=" + std::to_string(startParameterValue) + " to parameterValue=" +
								 std::to_string(endParameterValue) + ". End parameter value must be greater or equal to start parameter value.");
		}

		if(startParameterValue < parameterValues[0] || endParameterValue > parameterValues[parameterValues.size() - 1])
		{
			throw std::runtime_error("Length requested from parameterValue=" + std::to_string(startParameterValue) + " to parameterValue=" +
								 std::to_string(endParameterValue) + ". Parameter value must be between " + std::to_string(parameterValues[0]) + " and " +
								 std::to_string(parameterValues[parameterValues.size() - 1]) + ".");
		}

		return getLengthInInterval(startParameterValue, endParameterValue);
	}

	const T baseStart = parameterValues.front();
	const T baseEnd = parameterValues.back();
	const T period = baseEnd - baseStart;
	if(period <= 0)
	{
		throw std::runtime_error("Closed spline period is invalid.");
	}

	const T totalLength = getLengthInInterval(baseStart, baseEnd);
	if(std::abs(endParameterValue - startParameterValue) <= epsilon)
	{
		return 0;
	}

	T deltaParameter = endParameterValue - startParameterValue;
	if(deltaParameter < 0)
	{
		const T missingCycles = std::floor((-deltaParameter) / period) + 1;
		deltaParameter += missingCycles * period;
	}

	const T fullCycles = std::floor(deltaParameter / period);
	const T residualParameter = deltaParameter - (fullCycles * period);

	T length = fullCycles * totalLength;
	if(residualParameter <= epsilon)
	{
		return length;
	}

	const T normalizedStart = wrapParameterValue(startParameterValue);
	const T residualEnd = normalizedStart + residualParameter;
	if(residualEnd <= baseEnd)
	{
		length += getLengthInInterval(normalizedStart, residualEnd);
	}
	else
	{
		length += getLengthInInterval(normalizedStart, baseEnd);
		length += getLengthInInterval(baseStart, residualEnd - period);
	}

	return length;
}

template<typename T>
T simpline<T>::ParametrizedSpline::getLengthFromStart(const T& parameterValue) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get length of empty parametrized spline. Use non-default constructor to provide points.");
	}

	if(!closed)
	{
		if(parameterValue < parameterValues[0] || parameterValue > parameterValues[parameterValues.size() - 1])
		{
			throw std::runtime_error("Length-from-start requested at parameterValue=" + std::to_string(parameterValue) + ". Parameter value must be between " +
							 std::to_string(parameterValues[0]) + " and " + std::to_string(parameterValues[parameterValues.size() - 1]) + ".");
		}
		return getLength(parameterValues[0], parameterValue);
	}

	const T period = parameterValues.back() - parameterValues.front();
	if(period <= 0)
	{
		throw std::runtime_error("Closed spline period is invalid.");
	}

	const T totalLength = getLengthInInterval(parameterValues.front(), parameterValues.back());
	T normalizedDelta = std::fmod(parameterValue - parameterValues.front(), period);
	if(normalizedDelta < 0)
	{
		normalizedDelta += period;
	}

	if(std::abs(normalizedDelta) <= epsilon && (parameterValue - parameterValues.front()) > 0)
	{
		return totalLength;
	}

	return getLengthInInterval(parameterValues.front(), parameterValues.front() + normalizedDelta);
}

template<typename T>
T simpline<T>::ParametrizedSpline::getParameterValueAtLengthFromStart(const T& lengthFromStart) const
{
	if(lengthFromStart < 0)
	{
		throw std::runtime_error("Length from start must be non-negative.");
	}

	const T totalLength = getLength();
	if(!closed && lengthFromStart > totalLength)
	{
		throw std::runtime_error("Length from start must be lower than or equal to spline length.");
	}

	if(totalLength <= epsilon)
	{
		return parameterValues.front();
	}

	T targetLength = lengthFromStart;
	if(closed)
	{
		if(std::abs(lengthFromStart - totalLength) <= epsilon)
		{
			return parameterValues.back();
		}
		targetLength = std::fmod(lengthFromStart, totalLength);
		if(targetLength < 0)
		{
			targetLength += totalLength;
		}
	}

	T lowerBound = parameterValues.front();
	T upperBound = parameterValues.back();
	while(upperBound - lowerBound > epsilon)
	{
		const T middleValue = (lowerBound + upperBound) / 2;
		const T middleLength = getLengthFromStart(middleValue);
		if(middleLength < targetLength)
		{
			lowerBound = middleValue;
		}
		else
		{
			upperBound = middleValue;
		}
	}

	return (lowerBound + upperBound) / 2;
}

template<typename T>
T simpline<T>::ParametrizedSpline::getClosestParameterValue(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get closest point from empty parametrized spline. Use non-default constructor to provide points.");
	}

	const T start = parameterValues.front();
	const T end = parameterValues.back();
	const size_t sampleCount = std::max(std::max(static_cast<size_t>(3), coarseSamples), (parameterValues.size() - 1) * static_cast<size_t>(20));
	const T step = (end - start) / static_cast<T>(sampleCount - 1);

	size_t bestSampleIndex = 0;
	T bestDistanceSquared = std::numeric_limits<T>::infinity();
	for(size_t i = 0; i < sampleCount; i++)
	{
		const T parameterSample = start + (step * static_cast<T>(i));
		const T distanceSquared = (getValue(parameterSample) - queryPoint).squaredNorm();
		if(distanceSquared < bestDistanceSquared)
		{
			bestDistanceSquared = distanceSquared;
			bestSampleIndex = i;
		}
	}

	const size_t leftSampleIndex = (bestSampleIndex == 0) ? bestSampleIndex : bestSampleIndex - 1;
	const size_t rightSampleIndex = std::min(bestSampleIndex + 1, sampleCount - 1);
	T leftBound = start + (step * static_cast<T>(leftSampleIndex));
	T rightBound = start + (step * static_cast<T>(rightSampleIndex));

	if(rightBound <= leftBound)
	{
		return closed ? wrapParameterValue(leftBound) : leftBound;
	}

	const T phi = (1.0 + std::sqrt(5.0)) / 2.0;
	for(size_t i = 0; i < 50; i++)
	{
		const T middleLeft = rightBound - ((rightBound - leftBound) / phi);
		const T middleRight = leftBound + ((rightBound - leftBound) / phi);
		const T distanceLeftSquared = (getValue(middleLeft) - queryPoint).squaredNorm();
		const T distanceRightSquared = (getValue(middleRight) - queryPoint).squaredNorm();
		if(distanceLeftSquared < distanceRightSquared)
		{
			rightBound = middleRight;
		}
		else
		{
			leftBound = middleLeft;
		}
	}

	const T closestParameterValue = (leftBound + rightBound) / 2;
	return closed ? wrapParameterValue(closestParameterValue) : closestParameterValue;
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::getClosestPoint(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples) const
{
	return getValue(getClosestParameterValue(queryPoint, coarseSamples));
}

template<typename T>
T simpline<T>::ParametrizedSpline::getClosestDistance(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples) const
{
	return (getClosestPoint(queryPoint, coarseSamples) - queryPoint).norm();
}

template<typename T>
std::vector<typename simpline<T>::Vector3> simpline<T>::ParametrizedSpline::resampleByCount(const size_t& sampleCount, const bool& includeEndPoint) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot resample empty parametrized spline. Use non-default constructor to provide points.");
	}

	if(sampleCount == 0)
	{
		throw std::runtime_error("Sample count must be above zero.");
	}

	if(!closed && sampleCount < 2)
	{
		throw std::runtime_error("At least two samples are required for non-closed spline resampling by count.");
	}

	std::vector<typename simpline<T>::Vector3> resampledPoints;
	resampledPoints.reserve(sampleCount);

	const T totalLength = getLength();
	if(sampleCount == 1)
	{
		resampledPoints.push_back(getValue(parameterValues.front()));
		return resampledPoints;
	}

	const T denominator = includeEndPoint ? static_cast<T>(sampleCount - 1) : static_cast<T>(sampleCount);
	for(size_t i = 0; i < sampleCount; i++)
	{
		T fraction = static_cast<T>(i) / denominator;
		if(includeEndPoint && i == sampleCount - 1)
		{
			fraction = 1.0;
		}
		const T parameterValue = getParameterValueAtLengthFromStart(fraction * totalLength);
		resampledPoints.push_back(getValue(parameterValue));
	}

	return resampledPoints;
}

template<typename T>
std::vector<typename simpline<T>::Vector3> simpline<T>::ParametrizedSpline::resampleByDistance(const T& distanceStep, const bool& includeEndPoint) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot resample empty parametrized spline. Use non-default constructor to provide points.");
	}

	if(distanceStep <= 0)
	{
		throw std::runtime_error("Distance step must be above zero.");
	}

	std::vector<typename simpline<T>::Vector3> resampledPoints;
	const T totalLength = getLength();
	if(totalLength <= epsilon)
	{
		resampledPoints.push_back(getValue(parameterValues.front()));
		return resampledPoints;
	}

	for(T lengthFromStart = 0.0; lengthFromStart < totalLength; lengthFromStart += distanceStep)
	{
		resampledPoints.push_back(getValue(getParameterValueAtLengthFromStart(lengthFromStart)));
	}

	if(includeEndPoint)
	{
		resampledPoints.push_back(getValue(getParameterValueAtLengthFromStart(totalLength)));
	}

	return resampledPoints;
}

template<typename T>
T simpline<T>::ParametrizedSpline::getStartParameterValue() const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get start parameter value from empty parametrized spline. Use non-default constructor to provide points.");
	}

	return parameterValues.front();
}

template<typename T>
T simpline<T>::ParametrizedSpline::getEndParameterValue() const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get end parameter value from empty parametrized spline. Use non-default constructor to provide points.");
	}

	return parameterValues.back();
}

template<typename T>
bool simpline<T>::ParametrizedSpline::isClosed() const
{
	return closed;
}

template class simpline<float>::ParametrizedSpline;

template class simpline<double>::ParametrizedSpline;
