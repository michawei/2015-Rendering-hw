
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

struct Len {
    bool zero;
    float radius;
    float z;
    float n_ratio;
    float aperture;
};

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
	// RealisticCamera Private Methods
    bool ParsingDat(const string datFile);
    bool RaySphereIntersection(Point& O, Vector& D, Len len) const;
    Transform placeFilm;
    vector<Len> lens;
    float distance;
    float Hither;
    float Yon;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H