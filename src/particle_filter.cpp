/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define EPS 0.00001

std::default_random_engine gen;

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  if (is_initialized) {
    return;
  }

  // Initializing the number of particles
  num_particles = 100;

// normal distribution of distribution x with std_x 
normal_distribution<double> dist_x(x, std[0]);

// normal distribution of distribution y with std_y
normal_distribution<double> dist_y(y, std[1]);

//normal distribution of distribution theta with std_theta
normal_distribution<double> angle_theta(theta, std[2]);

  // Generate particles with normal distribution with mean on GPS values.
for (int i = 0; i < num_particles; ++i) {
	//Using struct to make a particle structure and assign every information about each particles
	Particle particle;
	particle.id = i;
	particle.x = dist_x(gen);
	particle.y = dist_y(gen);
	particle.theta = angle_theta(gen);

	// assign weight=1 to each particle 
	particle.weight = 1.0;

	// add particle to ParticleFilter class =>  std::vector<Particle> particles;
	// with this method, every particle and vecotr particles can be generated.
	//add structure into a vector.
	particles.push_back(particle);
	}
//after initialized, is_initialized should be true. If not, paricle fitler will always become initialized and uselessful.
is_initialized = true;
//I wonder if the 'return' can be added.
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	//add noise to velocity and yaw rate., seems like vector:particles can't multiply directly with motion model.
	// the standard variance although extract from sigma_pos from ParticleFilter::init, the each function inclass is saperated , so 
	// if we want to use the std_sigma of each parameters, we should re-extract again.

	// Engine for later generation of particles


	//double std_x = std_pos[0];
	//double std_y = std_pos[1];
	//double std_theta = std_pos[2];

	// normal distribution of distribution x with zero mean and each std.
	normal_distribution<double> dist_x(0, std_pos[0]);
	// normal distribution of distribution y with std_y
	normal_distribution<double> dist_y(0, std_pos[1]);
	//normal distribution of distribution theta with std_theta
	normal_distribution<double> angle_theta(0, std_pos[2]);

	// it needs for loop

	for (int i = 0; i < num_particles; i++) {
		//double theta = particles[i].theta;
		if (fabs(yaw_rate) >= EPS) {
			particles[i].x = particles[i].x + (velocity / yaw_rate)*(sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			particles[i].y = particles[i].y + (velocity / yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			particles[i].theta = particles[i].theta + yaw_rate * delta_t;
		}
		else {// theta doesn't change
			particles[i].x = particles[i].x + velocity * delta_t *cos(particles[i].theta);
			particles[i].y = particles[i].y + velocity * delta_t *sin(particles[i].theta);
		}

		//add noise  to each particle in particles.
		particles[i].x = particles[i].x + dist_x(gen);
		particles[i].y = particles[i].y + dist_y(gen);
		particles[i].theta = particles[i].theta + angle_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	//observed measurement?
	int n_observation = observations.size();
	//size of predicted measurement? 
	int n_predictions = predicted.size();

	for (int i = 0; i < n_observation; i++) {
		//for each observation
		//initializing the min distance as really big number
		double min_dis = numeric_limits<double>::max();

		//initializing the found map that is not in map , this is made for return the nearset measurement around GT.
		int id_in_map = -100;
		//complexity is o(ij);
		for (int j = 0; j < n_predictions; j++) {
			//my code with helper function
			//double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

			//dari code with direct calculation
			double xDistance = observations[i].x - predicted[j].x;
			double yDistance = observations[i].y - predicted[j].y;

			double distance = xDistance * xDistance + yDistance * yDistance;

			// if distance is smaller than the distance, then save the id , then iterate all the predicted value
			//finally find the most nearest precited to GT value. 
			if (distance < min_dis) {
				min_dis = distance;
				id_in_map = predicted[j].id;
			}
		}
		//assign the observed measurement to this particular landmark.
		//for vehicle, it means, this observation is belong to this landmark.
		observations[i].id = id_in_map;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	//std_landmark provides uncertainty of measurement of landmark in x,y direction . 
	//Why use there name here . ?
	double stdLandmarkRange = std_landmark[0];
	double stdLandmarkBearing = std_landmark[1];

	for (int i = 0; i < num_particles; i++) {
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		//find landmarks in vehicle sensing range
		double sensor_range_2 = sensor_range * sensor_range;
		vector<LandmarkObs> inRangeLandmarks;

		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			float landmarkX = map_landmarks.landmark_list[j].x_f;
			float landmarkY = map_landmarks.landmark_list[j].y_f;
			int id = map_landmarks.landmark_list[j].id_i;
			double dX = x - landmarkX;
			double dY = y - landmarkY;

			//in this step, in range is constructed. After this step, we only calculate the landmarks in the range. 
			if (dX*dX + dY * dY <= sensor_range_2) {
				inRangeLandmarks.push_back(LandmarkObs{ id, landmarkX, landmarkY });
			}
		}

		// Transfrom observation coodinates from vehicle coordinate to map (global) coordinate.
		vector<LandmarkObs> mappedObservations;
		//Rotation
		for (int j = 0; j< observations.size(); j++) {
			double xx = x + cos(theta)*observations[j].x - sin(theta) * observations[j].y;
			double yy = y + sin(theta)*observations[j].x + cos(theta) * observations[j].y;
			//using struct defined in helperfunction.h LandmarkObs, to make a after transition and rotation transformed observation data.
			//The @param observations is a noise mixed sensor measurement data.
			mappedObservations.push_back(LandmarkObs{ observations[j].id, xx, yy });
		}

		//Observation association with landmark
		dataAssociation(inRangeLandmarks, mappedObservations);

		//reset the weight, I think this line can be deleted since when particles' weights were set before with 1.0
		particles[i].weight = 1.0;

		//calculate the weights
		for (int j = 0; j < mappedObservations.size(); j++) {
			double observationX = mappedObservations[j].x;
			double observationY = mappedObservations[j].y;

			int landmarkId = mappedObservations[j].id;

			double landmarkX, landmarkY;

			int k = 0;
			int nLandmarks = inRangeLandmarks.size();
			bool found = false;
			while (!found && k < nLandmarks) {
				if (inRangeLandmarks[k].id == landmarkId) {
					found = true;
					landmarkX = inRangeLandmarks[k].x;
					landmarkY = inRangeLandmarks[k].y;
				}
				k++;
			}

			//calculating weight
			double dX = observationX - landmarkX;
			double dY = observationY - landmarkY;

			//Since we assume the correlation between x direction and y direction is not exist, then rho in wiki is zero.
			//my weight update
			//double weight = (1 / 2 * M_PI*stdLandmarkRange*stdLandmarkBearing)*exp(-(1 / 2)*(dX*dX / (stdLandmarkRange* stdLandmarkRange) + (dY*dY) / (stdLandmarkBearing*stdLandmarkBearing)));

			//dari weight update
			double weight = (1 / (2 * M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp(-(dX*dX / (2 * stdLandmarkRange*stdLandmarkRange) + (dY*dY / (2 * stdLandmarkBearing*stdLandmarkBearing))));

			//if weight equal to zero. then multiply to the EPS. But I dont know why it have to multiply with EPS. 
			// just make weight become zero can not work?
			if (weight == 0) {
				particles[i].weight = particles[i].weight*EPS;
			}

			//if weight doesn't equal to zero, then weight should be multiply by i times. because it is multivariate define.
			else {
				particles[i].weight = particles[i].weight * weight;
			}
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // Get weights and max weight.
  vector<double> weights;
  double max_weight = numeric_limits<double>::min();
  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if ( particles[i].weight > max_weight ) {
      max_weight = particles[i].weight;
    }
  }

  // Creating distributions.
  uniform_real_distribution<float> dist_float(0.0, max_weight);
  uniform_int_distribution<int> dist_int(0, num_particles - 1);

  // Generating index.
  int index = dist_int(gen);

  double beta = 0.0;

  // the wheel
  vector<Particle> resampled_particles;
  for(int i = 0; i < num_particles; i++) {
    beta += dist_float(gen) * 2.0;
    while( beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampled_particles.push_back(particles[index]);
  }

  particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
