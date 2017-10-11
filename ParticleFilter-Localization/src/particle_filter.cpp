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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 20;
	particles.resize(num_particles);
	weights.resize(num_particles);

	default_random_engine gen;
	double sample_x, sample_y, sample_psi;
	//create a normal (Gaussian) distribution for x, y, psi.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);

	//Add random Gaussian noise to each particle
	for(int i=0; i<num_particles; ++i)
	{
		// Sample  and from these normal distrubtions 
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_psi = dist_psi(gen);

		Particle p;
		//p.id = i;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_psi;
		p.weight = 1;
		weights[i] = 1;

		particles[i] = p;		
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	double sample_x, sample_y, sample_psi;
	//create a normal (Gaussian) distribution for x, y, psi.
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_psi(0, std_pos[2]);

	for(int i=0; i<num_particles; i++)
	{
		Particle &p = particles[i];
		if(fabs(yaw_rate) < 0.00001)
		{
			//Ref: https://discussions.udacity.com/t/yaw-rate-theta-dot/243585/2
			p.x = p.x + velocity * delta_t * cos(p.theta);
			p.y = p.y + velocity * delta_t * sin(p.theta);
			//Theta will not change
		}
		else
		{
			p.x = (p.x + (velocity / yaw_rate) * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta)));
			p.y = (p.y + (velocity / yaw_rate) * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t)));
			p.theta = p.theta + yaw_rate * delta_t;
		}

		//Add random Gaussian noise to each particle
		p.x += dist_x(gen);
		p.y += dist_y(gen);
		p.theta += dist_psi(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Find the predicted measurement that is closest to each observed measurement and assign the 
	// observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	// implement this method and use it as a helper during the updateWeights phase.
	for(int i=0; i<observations.size(); i++)
	{
		int lm_id;
		double min_dist = std::numeric_limits<double>::max();
		for(int j=0; j<predicted.size(); j++)
		{
			double cur_dist = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
			if(cur_dist <= min_dist)
			{
				min_dist = cur_dist;
				lm_id = predicted[j].id;
			}
		}
		observations[i].id = lm_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	//   Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	//   NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	//Update particle weights
	for(int i=0; i<num_particles; i++)
	{
		Particle &p = particles[i];

		//Create valid map landmarks vector which are in sensor range of a given particle
		std::vector<LandmarkObs> map_land_marks_near_observed_land_marks;
		for(int j=0; j<map_landmarks.landmark_list.size(); j++)
		{
			double distance = dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, p.x, p.y);
			if(distance <= sensor_range)
			{
				//Consider land marks which are in sensor range
				LandmarkObs map_lm_ob;
				map_lm_ob.x = map_landmarks.landmark_list[j].x_f;
				map_lm_ob.y = map_landmarks.landmark_list[j].y_f;
				map_lm_ob.id = map_landmarks.landmark_list[j].id_i;
				map_land_marks_near_observed_land_marks.push_back(map_lm_ob);
			}
		}

		//Translate observed landmarks from car co-ordinate system to Map coordinate system
		//store it in global_observations vector
		std::vector<LandmarkObs> global_observations;
		for(int j=0; j<observations.size(); j++)
		{
			LandmarkObs trasnlated_ob;
			trasnlated_ob.x = p.x + (observations[j].x * cos(p.theta)) - (observations[j].y * sin(p.theta));
			trasnlated_ob.y = p.y + (observations[j].x * sin(p.theta)) + (observations[j].y * cos(p.theta)); 
			global_observations.push_back(trasnlated_ob);
		}

		//associate each obseravation with a landmark
		dataAssociation(map_land_marks_near_observed_land_marks, global_observations);

		//Calculate multi-variate gaussian probability & final weight of a particle
		double final_weight = 1;
		for(int j=0; j<global_observations.size(); j++)
		{
			LandmarkObs global_ob = global_observations[j];
			for(int k =0; k<map_land_marks_near_observed_land_marks.size(); k++)
			{
				if(map_land_marks_near_observed_land_marks[k].id == global_ob.id)
				{
					double mvg = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1])) * 
						exp(-1 * ((((global_ob.x - map_land_marks_near_observed_land_marks[k].x)*(global_ob.x - map_land_marks_near_observed_land_marks[k].x))/(2 * std_landmark[0] * std_landmark[0]))+
							(((global_ob.y - map_land_marks_near_observed_land_marks[k].y)*(global_ob.y - map_land_marks_near_observed_land_marks[k].y))/(2 * std_landmark[1] * std_landmark[1]))));
					//Multiply mvgs for all observations
					final_weight *= mvg;
				}
			}
		}

		//Update weights of particle
		p.weight = final_weight;
		weights[i] = final_weight;

		//Set assications, sense_x, sense_y for this particle
		std::vector<double> tmp_sense_x;
		std::vector<double> tmp_sense_y;
		std::vector<int> tmp_associations;
		for(int j=0; j<global_observations.size(); j++)
		{
			LandmarkObs tmp = global_observations[j];
			tmp_sense_x.push_back(tmp.x);
			tmp_sense_y.push_back(tmp.y);
			tmp_associations.push_back(tmp.id);
		}
		particles[i] = SetAssociations(particles[i], tmp_associations, tmp_sense_x, tmp_sense_y);
	}
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	// http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::vector<Particle> new_particles;
    new_particles.resize(num_particles);
    new_particles.clear();

    // Use discrete_distribution to obtain set of particles 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(weights.begin(), weights.end());
    for(int i=0; i<num_particles; ++i) 
    {
        new_particles.push_back(particles[d(gen)]);
    }
    particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	// particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
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
