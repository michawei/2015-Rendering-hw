#include "stdafx.h"
#include "cameras/realistic.h"
#include "sampler.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
    // YOUR CODE HERE -- build and store datastructures representing the given lens
    // Hither is the minimum distance where the object are seen in the camera windows
    Hither = hither;
    // Yon is the maximum distance
    Yon = yon;
    distance = 0.f;
    // Read Dat file
    if ( !ParsingDat(specfile) ) {
        cerr << "Cannot open dat file : " << specfile << endl;
        exit(1);
    }
    distance -= filmdistance;
    float diagonal = sqrtf(f->xResolution * f->xResolution + f->yResolution * f->yResolution);
    float ratio = filmdiag / diagonal;
    // and film placement.
    placeFilm = Translate(Vector(ratio * 0.5 * f->xResolution, -ratio * 0.5 * f->yResolution, distance)) *
    Scale(-ratio, ratio, 1.f);
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    // YOUR CODE HERE -- make that ray!
    // use sample->imageX and sample->imageY to get raster-space coordinates
    // of the sample point on the film.
    Point cameraPos;
    placeFilm(Point(sample.imageX, sample.imageY, 0.f), &cameraPos);
    
    int nowLen = lens.size() - 1;
    // use sample->lensU and sample->lensV to get a sample position on the lens
    // ppt Another Method
    float r = lens[nowLen].aperture * sqrtf(sample.lensU);
    float theta = 2 * M_PI * sample.lensV;
    float z = sqrtf(lens[nowLen].radius * lens[nowLen].radius - r * r);
    float lensZ;
    
    if (lens[nowLen].radius < 0.f)
        lensZ = lens[nowLen].z - lens[nowLen].radius - z;
    else
        lensZ = lens[nowLen].z - lens[nowLen].radius + z;
    
    Point hit = Point(r * cosf(theta), r * sinf(theta), lensZ);
    Vector rayT = hit - cameraPos;
    
    // Ray Tracing
    for (int i = nowLen; i >= 0; --i){
        if ( lens[i].zero ) {
            rayT = Normalize(rayT);
            float t = (lens[i].z - hit.z) / rayT.z;
            hit += t * rayT;
            if ( hit.x * hit.x + hit.y * hit.y > lens[i].aperture * lens[i].aperture )
                return 0.f;
        }
        else {
            if (!RaySphereIntersection(hit, rayT, lens[i]))
                return 0.f;
        }
    }
    
    ray->o = hit;
    ray->d = Normalize(rayT);
    // Hither is the minimum distance where the object are seen in the camera windows
    ray->mint = Hither;
    // Yon is the maximum distance
    ray->maxt = (Yon - Hither) / ray->d.z;
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    
    // weight
    float w = Normalize(hit - cameraPos).z;
    w = w * (w / fabsf(distance));
    w = w * w * (lens[0].aperture * lens[0].aperture * M_PI);
    
    return w;
}

bool RealisticCamera::ParsingDat(const string file){
    FILE *fp = fopen(file.c_str(), "r");
    char *line = NULL;
    size_t length = 0;
    vector<float> ns;
    ns.push_back(1.f);
    
    if (fp == NULL)
        return false;
    
    float radius, axpos, n, aperture;
    while (getline(&line, &length, fp) != -1) {
        //printf("%s", line);
        if (line[0]=='#')
            continue;
        sscanf(line, "%f%f%f%f", &radius, &axpos, &n, &aperture);
        Len len = {(n == 0), radius, distance, 0.f, aperture * 0.5f};
        distance -= axpos;
        lens.push_back(len);
        ns.push_back((n == 0) ? 1.f : n);
    }
    
    for ( int i=1 ; i!=ns.size() ; ++i )
        lens[i-1].n_ratio = ns[i] / ns[i-1];

    fclose(fp);
    
    return true;
}

bool RealisticCamera::RaySphereIntersection(Point& O, Vector& D, Len len) const {
    Point C = Point(0.f, 0.f, len.z - len.radius);
    Vector O_C = O - C;
    Vector a = Normalize(D);
    float b = Dot(O_C, a);
    float c = O_C.LengthSquared() - len.radius * len.radius;
    float det = b * b - c;
    float t;
    Vector N;
    if ( det < 0.f )
        return false;
    else {
        det = sqrtf(det);
        if ( len.radius < 0.f )
            t = - b - det;
        else
            t = - b + det;
    }
    
    O = O + t * a;
    if ( O.x * O.x + O.y * O.y > len.aperture * len.aperture )
        return false;
    
    if ( len.radius < 0.f )
        N = Normalize(O - C);
    else
        N = Normalize(C - O);
    
    // Heckber's Method
    float c1 = Dot(-a, N);
    float c2 = 1.f - len.n_ratio * len.n_ratio * (1.f - c1 * c1);
    if (c2 <= 0.f)
        return false;
    else
        c2 = sqrtf(c2);
    
    D = len.n_ratio * a + ( len.n_ratio * c1 - c2 ) * N;
    
    return true;
}

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}
