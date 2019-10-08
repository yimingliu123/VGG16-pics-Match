/**
 *
 *  This code is part of the publication:
 *
 *      "Contours, corners and T-junctions detection algorithm" by Antoni Buades,
 *      Rafael Grompone von Gioi, and Julia Navarro. Image Processing On Line, 2018.
 *
 *  Antoni Buades <toni.buades@uib.es> Universitat de les Illes Balears, Spain
 *  Rafael Grompone von Gioi <grompone@gmail.com> CMLA, ENS Cachan, France
 *  Julia Navarro <julia.navarro@uib.es> Universitat de les Illes Balears, Spain
 *
 *
 **/

//
//  libDetectionIpol.h
//  Sources1.3
//
//  Created by antoni buades on 04/09/15.
//  Copyright (c) 2015 Antoni Buades. All rights reserved.
//

#ifndef __libDetectionIpol__
#define __libDetectionIpol__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <limits.h>
#include <float.h>

#include "libBasic.h"


using namespace std;


#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )


typedef struct point_descriptor {
    int po;
    int px;
    int py;
    int ncorbes;
    float den;

} str_point_descriptor ;


typedef struct dijkstra_point_descriptor
{
    int po;
    int px;
    int py;
    int index;

} str_dijkstra_point_descriptor ;



//! Updated code
void compute_directional_convolutions(float *input, float **dimages, int nDir, float fSigma, float anglePrecision, int width, int height);

void compute_corners(float **dimages, float *localMax, vector<str_point_descriptor> &cCorners, int nDir, float fSigma, float minRespCorner, float fCorThreshold,int width, int height);

void lateral_inhibition_each_scale(float **g, float **gh, float sigma, int n, int X, int Y, float fAngleFiltPrec, float fThreshold, float anglePrec, int flagKeepValue);

void good_continuation_filter(float **input, float **output, float Sigma, float fAngleFiltPrec, int nDir, float anglePrec, int width, int height);

void grouping_gestalt_binomial(float **dimages, float **curvesID, float *localMax, float fThreshold, std::vector < std::vector <str_point_descriptor> > &llista, int width, int height, int nDir, char *filename, float ntestsM, float lnfaThr, int connectivity);

void compute_corners(float **dimages, float *localMax, vector<str_point_descriptor> &cCorners, int nDir, float fSigma, float minRespCorner, float fCorThreshold, float radForZone, int width, int height);

void differentiate_corners_junctions(vector<str_point_descriptor> &cCorners, float **curvesID, float fSigma, int nDir,  float rMult, int width, int height);



void clear_corners(vector<str_point_descriptor> &cCorners, float *localMax, float **curvesID, float fSigma, int nDir, float anglePrecision, float rMult, float radForZone, float fCorValidationDegree, int width, int height);


#endif

