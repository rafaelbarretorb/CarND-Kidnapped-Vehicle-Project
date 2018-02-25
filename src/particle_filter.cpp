/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random> // default_random_engine
#include <algorithm> // // std::min_element, std::max_element
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;

	std::default_random_engine gen;

	normal_distribution<double> N_x(x, std[0]);
	normal_distribution<double> N_y(y, std[1]);
	normal_distribution<double> N_theta(theta, std[2]);


	for (int i = 0; i < num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);
		weights.push_back(particle.weight);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	double x_f;
	double y_f;
	double theta_f;

	for (int i = 0; i < num_particles; i++)
	{
		// from Lesson 12: Motion Models
		if (yaw_rate == 0)
		{
			x_f = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			y_f = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			theta_f = particles[i].theta;
		}
		else
		{
			x_f = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			y_f = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			theta_f = particles[i].theta + yaw_rate*delta_t;
		}
		
		normal_distribution<double> N_x(x_f, std_pos[0]);
		normal_distribution<double> N_y(y_f, std_pos[1]);
		normal_distribution<double> N_theta(theta_f, std_pos[2]);

		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);
	}
}

void ParticleFilter::dataAssociation(Particle& particle, vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) 
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	vector<int> associations;
	vector<double> sense_x;
	vector<double> sense_y;

	vector<double> distance_list;
	double distance;

	for (int i = 0; i < predicted.size(); i++)
	{
		// Create a list of distances between a particular TOBS(predicted[i]) and all landmark associated with it. 
		for(int j = 0; j < observations.size(); j++)
		{
			distance = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
			distance_list.push_back(distance);
		}

		// Find the index of the smallest value in the list created above
		int minPos = 0;
	    for (int k = 0; k < distance_list.size(); k++)
	    {
	        if (distance_list[k] < distance_list[minPos]) // Found a smallest distance index
	            minPos = k;
    	}

    	// For each particle association (Landmarks close TOBS's)
    	associations.push_back(observations[minPos].id); // id
    	sense_x.push_back(observations[minPos].x);
    	sense_y.push_back(observations[minPos].y);

    	// Clear the distance list for the next TOBS
    	distance_list.clear();

    }
    // Set the three vectors in the struct Particle
    SetAssociations(particle, associations, sense_x, sense_y);

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const vector<LandmarkObs> &observations, const Map &map_landmarks) 
{
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


	double x_map;
	double y_map;
	double distance;

	vector<LandmarkObs> predicted;
	vector<LandmarkObs> landmark_select;

	bool cutOff_observations = false;

	for (int i = 0; i < num_particles; i++)
	{

		for (int z = 0; z < map_landmarks.landmark_list.size(); z++)
		{
			distance = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[z].x_f , map_landmarks.landmark_list[z].y_f);
			if (distance < sensor_range)
			{
				//
				LandmarkObs landmark;
				landmark.id = map_landmarks.landmark_list[z].id_i;
				landmark.x = map_landmarks.landmark_list[z].x_f;
				landmark.y = map_landmarks.landmark_list[z].y_f;

				landmark_select.push_back(landmark);
			}
		}

		if(observations.size() != landmark_select.size())
		{
			particles[i].weight = 0;
			weights[i] = 0;
		}
		else
		{

			// Cut off observations with length larger than sensor range???
			if (cutOff_observations)
			{
				for (int j = 0; j < observations.size(); j++)
				{
					distance = dist(particles[i].x, particles[i].y, observations[j].x , observations[j].x);
					if ( distance < sensor_range)
						predicted.push_back(observations[j]);
				}				
			}
			else
			{
				predicted = observations;
			}

			/*
				TRANSFORMATION
			*/
			for (int k = 0; k < predicted.size(); k++)
			{
				// Transform to map x coordinates
				// x_map = x_particle + x_obs*cos(theta) - y_obs*sin(theta)
				x_map = particles[i].x + predicted[k].x*cos(particles[i].theta) - predicted[k].y*sin(particles[i].theta);

				// Transform to map y coordinate
				// y_map = y_particle + x_obs*sin(theta) + y_obs*cos(theta)
				y_map = particles[i].y + predicted[k].x*sin(particles[i].theta) + predicted[k].y*cos(particles[i].theta);

				predicted[k].x = x_map;
				predicted[k].y = y_map;
				// predicted[j].id = observations[j].id
			}

			/*
				ASSOCIATION
			*/	

			// Call dataAssociation function
			dataAssociation(particles[i], predicted, landmark_select);

			/*
				UPDATE WEIGHTS

			Nearest-neighbors approach:

			Assign the closest landmark (id_i, x_f, y_f) in order in:

			particles[i].associations -> list of id'a of the closest landmarks for each TOBS
			particles[i].sense_x -> respective list of mu_x
			particles[i].sense_y -> respective list of mu_y

			*/
			// before the multiplication
			particles[i].weight = 1;

			for (int m = 0; m < predicted.size(); m++)// supondo que gravei ids das landmarks  das TOBS em particles[i].associations
			{
				double gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
		    	double exponent = (pow( (predicted[m].x - particles[i].sense_x[m]) , 2)/(2*pow(std_landmark[0], 2)) 
		    		+ pow( (predicted[m].y - particles[i].sense_y[m]) , 2)/(2*pow(std_landmark[1], 2)));		

				particles[i].weight *= gauss_norm*exp(-exponent);
			}

			// FINAL Weight
			weights[i] = particles[i].weight;

		}


		// Clear tmp vectors
		landmark_select.clear();
		predicted.clear();
	}
	// Normalize weights
	double sum_of_elems = accumulate(weights.begin(), weights.end(), 0.0);
	transform(weights.begin(), weights.end(), weights.begin(), std::bind1st(std::multiplies<double>(),(1/sum_of_elems)));
}

void ParticleFilter::resample() 
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;

	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resample_particles;

	for (int i = 0; i < num_particles; i++)
		resample_particles.push_back(particles[distribution(gen)]);

	particles = resample_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

   	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

    particle.associations = associations;
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