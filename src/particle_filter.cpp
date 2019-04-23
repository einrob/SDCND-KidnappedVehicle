/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 10; // TODO: Set the number of particles

  Particle tmp_particle;

  // Generate Normal Distribution with expected mean x,y,theata for Gaussian Noise
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta,  std[2]);


  for(int i=0; i < num_particles; i++)
  {
	  tmp_particle.id = i;
	  tmp_particle.weight = 1.0;
#ifdef ADD_NOISE
	  tmp_particle.x = dist_x(gen_);
	  tmp_particle.y = dist_y(gen_);
	  tmp_particle.theta = dist_theta(gen_);
#else
	  tmp_particle.x = x;
	  tmp_particle.y = y;
	  tmp_particle.theta = theta;

#endif

	  particles.push_back(tmp_particle);
  }

  is_initialized = true;
  std::cout << "IS INITIALIZED!" << std::endl;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

	std::default_random_engine gen;

	// Predict next position with noise

// I had a bug in the if statent below. I checked for yaw_rate > 0.000001 || yaw_rate < 0.000001 which lead to
// division by zero, because yaw_rate == ~0 was not catched of course
//	std::cout << "Particles before prediction: " << std::endl;
//	for (unsigned int i=0; i < particles.size(); i++)
//	{
//		std::cout << "Particle #" << i<< " X: " << particles[i].x << " Y: " <<  particles[i].y << std::endl;
//	}
//	std::cout << "Particles after prediction: " << std::endl;

	if (yaw_rate > 0.000001 || yaw_rate < -0.000001) // Yaw rate is greater zero
	{
		for (unsigned int i=0; i < particles.size(); i++)
		{
			particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t));
			particles[i].theta += yaw_rate*delta_t;

#ifdef ADD_NOISE
			std::normal_distribution<double> dist_x(0.0, std_pos[0]);
			std::normal_distribution<double> dist_y(0.0, std_pos[1]);
			std::normal_distribution<double> dist_theta(0.0,  std_pos[2]);

			particles[i].x +=  dist_x(gen);
			particles[i].y +=  dist_y(gen);
			particles[i].theta +=  dist_theta(gen);
			//std::cout << "Particle #" << i<< " X: " << particles[i].x << " Y: " <<  particles[i].y << std::endl;

#endif
		}
	}
	else // Yaw rate is zero, use constant yaw rate motion model
	{
		for (unsigned int i=0; i < particles.size(); i++)
		{
			particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			particles[i].y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			// theta unchanged

#ifdef ADD_NOISE
			std::normal_distribution<double> dist_x(0.0, std_pos[0]);
		    std::normal_distribution<double> dist_y(0.0, std_pos[1]);
			std::normal_distribution<double> dist_theta(0.0,  std_pos[2]);

			particles[i].x +=  dist_x(gen);
			particles[i].y +=  dist_y(gen);
			particles[i].theta +=  dist_theta(gen);
			//std::cout << "Particle #" << i<< " X: " << particles[i].x << " Y: " <<  particles[i].y << std::endl;

#endif
		}
	}


}

// Prediced measurements == praticle to map measurements in MAP coordinates
// observations == particle observations in MAP coordinates
void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */


	double closest_distance = 9999999.0;
	double obs_distance = 0.0;


	for(unsigned int obs_index = 0; obs_index < observations.size(); obs_index++)
	{
		 closest_distance = 9999999.0;

		for(unsigned int predicted_index = 0; predicted_index < predicted.size(); predicted_index++)
		{
			obs_distance = calcEuclideanDistance2D( predicted[predicted_index].x,
													predicted[predicted_index].y,
													observations[obs_index].x,
													observations[obs_index].y);

			if(obs_distance < closest_distance && obs_distance != 0)
			{
				observations[obs_index].id = predicted[predicted_index].id;
				observations[obs_index].map_index = predicted_index;
				closest_distance = obs_distance;
			}
		}
	}

//	std::cout << "Final observation association " << std::endl;
//	for(int obs = 0; obs < observations.size(); obs++)
//	{
//		std::cout << "Obs Index: " << obs << " Associated ID: " << observations[obs].id << std::endl;
//	}
}


double ParticleFilter::calcEuclideanDistance2D(double x1, double y1, double x2, double y2)
{
	return sqrt(std::pow(x1 - x2,2) + pow(y1 - y2,2));
}

vector<double> ParticleFilter::transformCarToMapCooridinate(double x_c, double y_c, double x_p, double y_p, double theta)
{
	vector<double> xy_map;

	xy_map.push_back(x_p + (cos(theta)*x_c) - sin(theta)*y_c);
	xy_map.push_back(y_p + (sin(theta)*x_c) + cos(theta)*y_c);

	return xy_map;
}


// xy obs
// x_mu, y_mu landmark

double ParticleFilter::calcMultivariateGaussian(double x, double y, float x_mu, float y_mu, double std_landmark[])
{
	double numerator;

	numerator = exp(-1.0 * ((pow((x - x_mu), 2)/(2.0 *  std_landmark[0])) + (pow((y - y_mu), 2)/(2.0 * std_landmark[1]))));

	return (numerator / (2 * M_PI * std_landmark[0] * std_landmark[1]));
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

	vector<LandmarkObs> particle_observed_landmark_positions; // Pseudo measurements to map landmarks in particles coordinate system
	vector<LandmarkObs> particle_predicted_landmark_positions;

	LandmarkObs tmp_landmark_obs;

	std::vector<double> tmp_point_xy;

	double measurement_dist = 0;
	double tmp_weight = 0.0;
	double weight_sum = 0.0;

	//std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx " << std::endl;



	for(unsigned int i = 0; i < particles.size(); i++)
	{


		//std::cout << "Particle " << i << std::endl;

		// Reste measurements
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();
		particles[i].associations.clear();

		particle_predicted_landmark_positions.clear();
		particle_observed_landmark_positions.clear();

		//std::cout << "Particle #" << i<< " X: " << particles[i].x << " Y: " <<  particles[i].y << std::endl;


		// Iterating through all map landmarks
		for(unsigned int map_landmark_index =0; map_landmark_index < map_landmarks.landmark_list.size(); map_landmark_index++)
		{
			measurement_dist = calcEuclideanDistance2D( (double)map_landmarks.landmark_list[map_landmark_index].x_f,
														(double)map_landmarks.landmark_list[map_landmark_index].y_f,
														particles[i].x,
														particles[i].y);

			if(measurement_dist <=  sensor_range)
			{
				tmp_landmark_obs.x = (double)map_landmarks.landmark_list[map_landmark_index].x_f; // WUUUUT why is landmark map not compatible with landmark obs!?!?
				tmp_landmark_obs.y = (double)map_landmarks.landmark_list[map_landmark_index].y_f;

				// Add measurement noise
				std::normal_distribution<double> dist_x(tmp_landmark_obs.x, std_landmark[0]);
				std::normal_distribution<double> dist_y(tmp_landmark_obs.y, std_landmark[1]);

				tmp_landmark_obs.x =  dist_x(gen_);
				tmp_landmark_obs.y =  dist_y(gen_);

				tmp_landmark_obs.id = map_landmarks.landmark_list[map_landmark_index].id_i;

				particle_predicted_landmark_positions.push_back(tmp_landmark_obs);


//				std::cout << "Pred Landmark sense X:  " << tmp_landmark_obs.x << std::endl;
//				std::cout << "Pred Landmark sense Y:  " << tmp_landmark_obs.y << std::endl;
//				std::cout << "Predicted landmark distance: " << measurement_dist << std::endl;

			}

		}

		std::cout << "# of predicted landmark observations: " << particle_predicted_landmark_positions.size() << std::endl;

		if(particle_predicted_landmark_positions.size() > 0)
		{
	//		std::cout << "PX: " << particles[i].x << std::endl;
	//		std::cout << "PY: " << particles[i].y << std::endl;
	//		std::cout << "PTheta: " << particles[i].theta << std::endl;
	//		std::cout << std::endl;

			// transform all current observations of particle into map coordinate system

			for(unsigned int obs =0; obs < observations.size(); obs++)
			{

				tmp_point_xy = transformCarToMapCooridinate(observations[obs].x,
															observations[obs].y,
															particles[i].x,
															particles[i].y,
															particles[i].theta);


				particles[i].sense_x.push_back(tmp_point_xy[0]);
				particles[i].sense_y.push_back(tmp_point_xy[1]);

				tmp_landmark_obs.x = tmp_point_xy[0];  // Why the do the particles have single vectors for each observation!?
				tmp_landmark_obs.y = tmp_point_xy[1];


				particle_observed_landmark_positions.push_back(tmp_landmark_obs);

	//			std::cout << "Obs Landmark sense X:  " << tmp_landmark_obs.x << std::endl;
	//			std::cout << "Obs Landmark sense Y:  " << tmp_landmark_obs.y << std::endl;

			}

			// Now find associated landmark to each observation the current particle

			dataAssociation(particle_predicted_landmark_positions, particle_observed_landmark_positions);


			// calculate likelihood for each observation of a particle
			// take the product of all likelihoods

			particles[i].weight = 1.0;
			tmp_weight = 0.0;


			for(unsigned int obs = 0; obs < observations.size(); obs++)
			{
				//std::cout << "SenseX: " << particles[i].sense_x[obs] << " SenseY: " << particles[i].sense_y[obs] << " map_id: " << particle_observed_landmark_positions[obs].id << " Map index: "<< particle_observed_landmark_positions[obs].map_index<< std::endl;
				particles[i].associations.push_back(particle_observed_landmark_positions[obs].id);



				tmp_weight = calcMultivariateGaussian( 	particles[i].sense_x[obs],
														particles[i].sense_y[obs] ,
														particle_predicted_landmark_positions[particle_observed_landmark_positions[obs].map_index].x,
														particle_predicted_landmark_positions[particle_observed_landmark_positions[obs].map_index].y,
														std_landmark);


//
//				std::cout << 	"obsX: " << particles[i].sense_x[obs] <<
//								"obsY: " << particles[i].sense_y[obs] << std::endl <<
//								"LX: " 	<<	particle_predicted_landmark_positions[particle_observed_landmark_positions[obs].map_index].x <<
//								"LY: "  << particle_predicted_landmark_positions[particle_observed_landmark_positions[obs].map_index].y << std::endl;
//
				if(tmp_weight > 0)
				{
					particles[i].weight *= tmp_weight;
					//std::cout << "Obs"<<obs<<" weight: " << tmp_weight << std::endl << std::endl;
				}

			}
		}
		else
		{
			particles[i].weight = 0.0000001;
		}


		weight_sum += particles[i].weight;

//		std::cout << "Particle weight: " << particles[i].weight << std::endl;
//		std::cout << "-------------------------------------------------------------------" << std::endl;
	}

	  for (unsigned int i = 0; i < particles.size(); i++)
	  {
		particles[i].weight /= weight_sum;
	  }


}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

	std::vector<Particle> new_particles;
	int index = std::rand() % particles.size();

	double max_weight =   -1.0;
	double beta = 0.0;
	float random;
	std::vector<double> weights;

	//std::cout << "-------------------------------------------------------------" << std::endl;

	for(unsigned int i=0; i<particles.size(); i++)
	{
		if(max_weight<particles[i].weight)
			max_weight = particles[i].weight;

		//std::cout << particles[i].weight << std::endl;
		weights.push_back(particles[i].weight);
	}
//
//	std::cout << "Max  Weight: " << max_weight << std::endl;
//	std::cout << "-------------------------------------------------------------" << std::endl;
//
//	std::cout << std::endl;


	//std::cout << "MAX WEIGHT: " << max_weight << std::endl;

	for(unsigned int i=0; i<weights.size();i++)
	{

		random = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		//std::cout << random << std::endl;
		beta +=random * 2.0 * max_weight;

		while(beta > weights[index])
		{
		   beta -= weights[index];
		   index = (index + 1) % weights.size();
		}

		new_particles.push_back(particles[index]);

	}

	particles = new_particles;

/*
    for i in range(N):
        w.append(p[i].measurement_prob(Z))

    p3 = []
    index = int(random.random() * N)
    beta = 0.0
    mw = max(w)
    for i in range(N):
        beta += random.random() * 2.0 * mw
        while beta > w[index]:
            beta -= w[index]
            index = (index + 1) % N
        p3.append(p[index])
    p = p3
*/




}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}


