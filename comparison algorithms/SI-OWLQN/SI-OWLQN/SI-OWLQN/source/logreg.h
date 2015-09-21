#pragma once

#include <deque>
#include <vector>
#include <cmath>
#include <iostream>

#include "OWLQN.h"

class LogisticRegressionProblem {
	std::deque<size_t> indices;
	std::deque<float> values;
	std::deque<size_t> instance_starts;
	std::deque<bool> labels;
	size_t numFeats;

public:
	LogisticRegressionProblem(size_t numFeats) : numFeats(numFeats) {
		instance_starts.push_back(0);
	}

	LogisticRegressionProblem(const char* mat, const char* labels);
	void AddInstance(const std::deque<size_t>& inds, const std::deque<float>& vals, bool label);
	void AddInstance(const std::vector<float>& vals, bool label);
	double ScoreOf(size_t i, const std::vector<double>& weights) const;

	bool LabelOf(size_t i) const {
		return labels[i];
	}

	void AddMultTo(size_t i, double mult, std::vector<double>& vec) const {
		if (labels[i]) mult *= -1;
		for (size_t j=instance_starts[i]; j < instance_starts[i+1]; j++) {
			size_t index = (indices.size() > 0) ? indices[j] : j - instance_starts[i];
			vec[index] += mult * values[j];
		}
	}

	size_t NumInstances() const {
		return labels.size();
	}

	size_t NumFeats() const {
		return numFeats;
	}
};

struct LogisticRegressionObjective : public DifferentiableFunction {
	const LogisticRegressionProblem& problem;
	const double l2weight;

	LogisticRegressionObjective(const LogisticRegressionProblem& p, double l2weight = 0) : problem(p), l2weight(l2weight) { }

	double Eval(const DblVec& input,  DblVec& gradient);

};
