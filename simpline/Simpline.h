#ifndef SIMPLINE_SIMPLINE_H
#define SIMPLINE_SIMPLINE_H

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <cstddef>

template<typename T>
struct simpline
{
	typedef typename Eigen::Matrix<T, 3, 1> Vector3;
	typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> VectorX;
	
	class ParametrizedSpline
	{
	public:
		ParametrizedSpline();
		
		ParametrizedSpline(const std::vector<T>& parameterValues, const std::vector<simpline<T>::Vector3>& points, const bool& closed = false);
		
		simpline<T>::Vector3 getValue(const T& parameterValue) const;
		
		simpline<T>::Vector3 getGradient(const T& parameterValue) const;

		simpline<T>::Vector3 getSecondDerivative(const T& parameterValue) const;
		
		T getLength() const;
		
		T getLength(const T& startParameterValue, const T& endParameterValue) const;

		T getLengthFromStart(const T& parameterValue) const;

		T getClosestParameterValue(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples = 200) const;

		simpline<T>::Vector3 getClosestPoint(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples = 200) const;

		T getClosestDistance(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples = 200) const;

		std::vector<simpline<T>::Vector3> resampleByCount(const size_t& sampleCount, const bool& includeEndPoint = true) const;

		std::vector<simpline<T>::Vector3> resampleByDistance(const T& distanceStep, const bool& includeEndPoint = true) const;

		T getStartParameterValue() const;

		T getEndParameterValue() const;

		bool isClosed() const;
	
	private:
		T wrapParameterValue(const T& parameterValue) const;
		T getLengthInInterval(const T& startParameterValue, const T& endParameterValue) const;
		T getParameterValueAtLengthFromStart(const T& lengthFromStart) const;
		size_t getSegmentIndex(const T& parameterValue) const;
		simpline<T>::Vector3 computeFiniteDifference(const size_t& startPointIndex, const size_t& endPointIndex);
		
		std::vector<T> parameterValues;
		std::vector<simpline<T>::Vector3> points;
		std::vector<simpline<T>::Vector3> firstDerivatives;
		std::vector<simpline<T>::Vector3> secondDerivatives;
		std::vector<simpline<T>::Vector3> thirdDerivatives;
		std::map<std::pair<size_t, size_t>, simpline<T>::Vector3> finiteDifferences;
		bool closed;
	};
	
	class ConstantSpeedSpline
	{
	public:
		ConstantSpeedSpline();
		
		ConstantSpeedSpline(const std::vector<simpline<T>::Vector3>& points, const T& speed, const bool& closed = false);
		
		simpline<T>::Vector3 getValue(const T& time) const;
		
		simpline<T>::Vector3 getGradient(const T& time) const;

		simpline<T>::Vector3 getSecondDerivative(const T& time) const;
		
		T getLength() const;
		
		T getLength(const T& startTime, const T& endTime) const;

		T getLengthFromStart(const T& time) const;

		T getClosestTime(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples = 200) const;

		simpline<T>::Vector3 getClosestPoint(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples = 200) const;

		T getClosestDistance(const simpline<T>::Vector3& queryPoint, const size_t& coarseSamples = 200) const;

		std::vector<simpline<T>::Vector3> resampleByCount(const size_t& sampleCount, const bool& includeEndPoint = true) const;

		std::vector<simpline<T>::Vector3> resampleByDistance(const T& distanceStep, const bool& includeEndPoint = true) const;
		
		T getSpeed() const;
		
		T getDuration() const;

		bool isClosed() const;
	
	private:
		T computeParameterValue(const T& time) const;
		
		ParametrizedSpline parametrizedSpline;
		T speed;
		T duration;
		bool closed;
		bool initialized;
	};
};

#endif
