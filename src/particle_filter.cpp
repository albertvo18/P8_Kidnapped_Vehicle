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
#define NUM_OF_PARTICLES 100

using namespace std;
static default_random_engine random_gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  if (is_initialized) {
    return;
  }

  // Initializing the number of particles
  int num_of_particles = NUM_OF_PARTICLES;

  //  standard deviations
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  // Creating normal distributions for sensor noise
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // Generate particles with normal distribution 
  // with mean on GPS values.
  for (int i = 0; i < NUM_OF_PARTICLES; i++) {

    Particle particle;
    particle.id = i;
    particle.x = dist_x(random_gen);
    particle.y = dist_y(random_gen);
    particle.theta = dist_theta(random_gen);
    particle.weight = 1.0;

    particles.push_back(particle);
	}

  // The filter is now initialized.
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  const double vel_X_dt = velocity * delta_t;
  const double yaw_rate_X_dt = yaw_rate * delta_t;
  const double vel_div_yaw_rate = velocity/yaw_rate;





  //  standard deviations
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

  // Creating normal distributions
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  // Calculate new state.
  for (int i = 0; i < NUM_OF_PARTICLES; i++) {

  	double theta = particles[i].theta;

    // yaw rate not changing
    if ( fabs(yaw_rate) < EPS ) { 

      particles[i].x += vel_X_dt * cos( theta );
      particles[i].y += vel_X_dt * sin( theta );

      // yaw continue to be the same.
    } else {


      particles[i].x += vel_div_yaw_rate * ( sin( theta + yaw_rate_X_dt ) - sin( theta ) );
      particles[i].y += vel_div_yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
      particles[i].theta += yaw_rate_X_dt;


    }

    // Adding noise.
    particles[i].x += dist_x(random_gen);
    particles[i].y += dist_y(random_gen);
    particles[i].theta += dist_theta(random_gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

  unsigned int num_observations = observations.size();
  unsigned int num_predictions = predicted.size();

  for (unsigned int i = 0; i < num_observations; i++) { 

    // Initialize min distance
    double min_distance = numeric_limits<double>::max();

    // Initialize the found map 
    int map_id = -1;

    // iterate through predictions
    for (unsigned j = 0; j < num_predictions; j++ ) { 
      double x_distance = observations[i].x - predicted[j].x;
      double y_distance = observations[i].y - predicted[j].y;

      double distance = x_distance * x_distance + y_distance * y_distance;

      // If the "distance" is less than min, stored the id and update min.
      if ( distance < min_distance ) {
        min_distance = distance;
        map_id = predicted[j].id;
      }
    }

    // Update the observation identifier.
    observations[i].id = map_id;
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
  double stdLandmarkRange = std_landmark[0];
  double stdLandmarkBearing = std_landmark[1];

  for (int i = 0; i < NUM_OF_PARTICLES; i++) {

    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;

    // Find landmarks in particle's range.
    double sensor_range_2 = sensor_range * sensor_range;
    vector<LandmarkObs> inRangeLandmarks;
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      float landmark_X = map_landmarks.landmark_list[j].x_f;
      float landmark_Y = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      double dX = x - landmark_X;
      double dY = y - landmark_Y;
      if ( dX*dX + dY*dY <= sensor_range_2 ) {
        inRangeLandmarks.push_back(LandmarkObs{ id, landmark_X, landmark_Y });
      }
    }

    // Transform observation coordinates.
    vector<LandmarkObs> mappedObservations;
    for(unsigned int j = 0; j < observations.size(); j++) {
      double xx = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double yy = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      mappedObservations.push_back(LandmarkObs{ observations[j].id, xx, yy });
    }

    // Observation association to landmark.
    dataAssociation(inRangeLandmarks, mappedObservations);

    // Reseting weight.
    particles[i].weight = 1.0;
    // Calculate weights.
    for(unsigned int j = 0; j < mappedObservations.size(); j++) {
      double observation_X = mappedObservations[j].x;
      double observation_Y = mappedObservations[j].y;

      int landmarkId = mappedObservations[j].id;

      double landmark_X, landmark_Y;
      unsigned int k = 0;
      unsigned int nLandmarks = inRangeLandmarks.size();
      bool found = false;
      while( !found && k < nLandmarks ) {
        if ( inRangeLandmarks[k].id == landmarkId) {
          found = true;
          landmark_X = inRangeLandmarks[k].x;
          landmark_Y = inRangeLandmarks[k].y;
        }
        k++;
      }

      // Calculating weight.
      double dX = observation_X - landmark_X;
      double dY = observation_Y - landmark_Y;

      double weight = ( 1/(2*M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp( -( dX*dX/(2*stdLandmarkRange*stdLandmarkRange) + (dY*dY/(2*stdLandmarkBearing*stdLandmarkBearing)) ) );
      if (weight == 0) {
        particles[i].weight *= EPS;
      } else {
        particles[i].weight *= weight;
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
  for(int i = 0; i < NUM_OF_PARTICLES; i++) {
    weights.push_back(particles[i].weight);
    if ( particles[i].weight > max_weight ) {
      max_weight = particles[i].weight;
    }
  }

  // Creating distributions.
  uniform_real_distribution<double> distDouble(0.0, max_weight);
  uniform_int_distribution<int> distInt(0, NUM_OF_PARTICLES - 1);

  // Generating index.
  int index = distInt(random_gen);

  double beta = 0.0;

  // the wheel
  vector<Particle> resampledParticles;
  for(int i = 0; i < NUM_OF_PARTICLES; i++) {
    beta += distDouble(random_gen) * 2.0;
    while( beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % NUM_OF_PARTICLES;
    }
    resampledParticles.push_back(particles[index]);
  }

  particles = resampledParticles;
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
    s = s.substr(0, s.length()-1);  
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  
    return s;
}
