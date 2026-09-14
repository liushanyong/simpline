#include "Simpline.h"
#include "Constants.h"
#include <cmath>
#include <stdexcept>

template<typename T>
simpline<T>::ConstantSpeedSpline::ConstantSpeedSpline():
		parametrizedSpline(), speed(0), duration(0), closed(false), initialized(false)
{
}

template<typename T>
simpline<T>::ConstantSpeedSpline::ConstantSpeedSpline(const std::vector<simpline<T>::Vector3>& points, const T& speed, const bool& closed):
		parametrizedSpline(), speed(speed), duration(0), closed(closed), initialized(false)
{
	if(points.size() < 2)
	{
		throw std::runtime_error("Point list must contain at least two items!");
	}

	if(speed <= 0)
	{
		throw std::runtime_error("Speed must be above 0!");
	}

	std::vector<simpline<T>::Vector3> splinePoints = points;
	if(this->closed && splinePoints.size() > 2 && (splinePoints.front() - splinePoints.back()).norm() <= epsilon)
	{
		splinePoints.pop_back();
	}

	if(this->closed && splinePoints.size() < 3)
	{
		throw std::runtime_error("Closed spline requires at least three distinct points.");
	}

	std::vector<T> parameterValues = { 0 };
	for(size_t i = 1; i < splinePoints.size(); i++)
	{
		parameterValues.push_back(parameterValues[i - 1] + (splinePoints[i] - splinePoints[i - 1]).norm());
	}

	parametrizedSpline = ParametrizedSpline(parameterValues, splinePoints, this->closed);
	duration = parametrizedSpline.getLength() / speed;
	initialized = true;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::computeParameterValue(const T& time) const
{
	if(time <= 0)
	{
		return parametrizedSpline.getStartParameterValue();
	}

	if(time >= duration)
	{
		return closed ? parametrizedSpline.getStartParameterValue() : parametrizedSpline.getEndParameterValue();
	}

	const T targetLength = speed * time;
	T lowerBound = parametrizedSpline.getStartParameterValue();
	T upperBound = parametrizedSpline.getEndParameterValue();

	while(upperBound - lowerBound > epsilon)
	{
		const T middleValue = (lowerBound + upperBound) / 2;
		const T middleLength = parametrizedSpline.getLengthFromStart(middleValue);
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
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getValue(const T& time) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get value from empty constant-speed spline. Use non-default constructor to provide points.");
	}

	if(time < 0.0 || time > duration)
	{
		throw std::runtime_error("Value requested at time=" + std::to_string(time) + ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}

	return parametrizedSpline.getValue(computeParameterValue(time));
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getGradient(const T& time) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get gradient from empty constant-speed spline. Use non-default constructor to provide points.");
	}

	if(time < 0.0 || time > duration)
	{
		throw std::runtime_error("Gradient requested at time=" + std::to_string(time) + ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}

	return parametrizedSpline.getGradient(computeParameterValue(time)).normalized() * speed;
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getSecondDerivative(const T& time) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get second derivative from empty constant-speed spline. Use non-default constructor to provide points.");
	}

	if(time < 0.0 || time > duration)
	{
		throw std::runtime_error("Second derivative requested at time=" + std::to_string(time) + ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}

	if(duration <= epsilon)
	{
		return simpline<T>::Vector3::Zero();
	}

	const T centralStep = std::min(std::max(duration * static_cast<T>(1e-3), static_cast<T>(epsilon)), duration / 2);
	if(time >= centralStep && (time + centralStep) <= duration)
	{
		return (getGradient(time + centralStep) - getGradient(time - centralStep)) / (2.0 * centralStep);
	}

	const T oneSidedStep = std::min(std::max(duration * static_cast<T>(1e-3), static_cast<T>(epsilon)), duration);
	if((time + oneSidedStep) <= duration)
	{
		return (getGradient(time + oneSidedStep) - getGradient(time)) / oneSidedStep;
	}

	return (getGradient(time) - getGradient(time - oneSidedStep)) / oneSidedStep;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getLength() const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get length of empty constant-speed spline. Use non-default constructor to provide points.");
	}

	return duration * speed;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getLength(const T& startTime, const T& endTime) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get length of empty constant-speed spline. Use non-default constructor to provide points.");
	}

	if(startTime > endTime)
	{
		throw std::runtime_error("Length requested from time=" + std::to_string(startTime) + " to time=" + std::to_string(endTime) +
							 ". End time must be greater or equal to start time.");
	}

	if(startTime < 0.0 || endTime > duration)
	{
		throw std::runtime_error("Length requested from time=" + std::to_string(startTime) + " to time=" + std::to_string(endTime) +
							 ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}

	return (endTime - startTime) * speed;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getLengthFromStart(const T& time) const
{
	return getLength(0.0, time);
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getClosestTime(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get closest point from empty constant-speed spline. Use non-default constructor to provide points.");
	}

	const T parameterValue = parametrizedSpline.getClosestParameterValue(queryPoint, coarseSamples);
	return parametrizedSpline.getLengthFromStart(parameterValue) / speed;
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getClosestPoint(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get closest point from empty constant-speed spline. Use non-default constructor to provide points.");
	}

	return parametrizedSpline.getClosestPoint(queryPoint, coarseSamples);
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getClosestDistance(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get closest distance from empty constant-speed spline. Use non-default constructor to provide points.");
	}

	return parametrizedSpline.getClosestDistance(queryPoint, coarseSamples);
}

template<typename T>
std::vector<typename simpline<T>::Vector3> simpline<T>::ConstantSpeedSpline::resampleByCount(const size_t& sampleCount, const bool& includeEndPoint) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot resample empty constant-speed spline. Use non-default constructor to provide points.");
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
	if(sampleCount == 1)
	{
		resampledPoints.push_back(getValue(0.0));
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
		resampledPoints.push_back(getValue(fraction * duration));
	}

	return resampledPoints;
}

template<typename T>
std::vector<typename simpline<T>::Vector3> simpline<T>::ConstantSpeedSpline::resampleByDistance(const T& distanceStep, const bool& includeEndPoint) const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot resample empty constant-speed spline. Use non-default constructor to provide points.");
	}

	if(distanceStep <= 0)
	{
		throw std::runtime_error("Distance step must be above zero.");
	}

	std::vector<typename simpline<T>::Vector3> resampledPoints;
	const T totalLength = getLength();
	if(totalLength <= epsilon)
	{
		resampledPoints.push_back(getValue(0.0));
		return resampledPoints;
	}

	for(T distance = 0.0; distance < totalLength; distance += distanceStep)
	{
		resampledPoints.push_back(getValue(distance / speed));
	}

	if(includeEndPoint)
	{
		resampledPoints.push_back(getValue(duration));
	}

	return resampledPoints;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getSpeed() const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get speed of empty constant-speed spline. Use non-default constructor to provide points.");
	}

	return speed;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getDuration() const
{
	if(!initialized)
	{
		throw std::runtime_error("Cannot get duration of empty constant-speed spline. Use non-default constructor to provide points.");
	}

	return duration;
}

template<typename T>
bool simpline<T>::ConstantSpeedSpline::isClosed() const
{
	return closed;
}

template class simpline<float>::ConstantSpeedSpline;

template class simpline<double>::ConstantSpeedSpline;
