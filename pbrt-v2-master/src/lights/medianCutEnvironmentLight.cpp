
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// lights/medianCutEnvironmentLight.cpp*
#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"
#include <math.h>

// medianCutEnvironmentLight Utility Classes
struct InfiniteAreaCube {
    // MedianCutEnvironmentLight Public Methods
    InfiniteAreaCube(const MedianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};



// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete distribution;
    delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
    // delete[] texels;
    // Initialize sampling PDFs for infinite area light
        
    float solidAngle = ((2.f * M_PI) / float(width - 1)) * (M_PI / float(height - 1));
        
    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta;
            texels[u+v*width] *= (solidAngle * sinTheta);
        }
    }
    
    /* Here */
    int totalSize = width * height;
    float *table = new float[totalSize];
    for (int i = 0 ; i < height ; ++i){
        for ( int j = 0 ; j < width ; ++j ){
                
            int self = i * width + j;
            int left = self - 1;
            int up = self - width;
                
            if ( i == 0 && j == 0 )
                table[self] = img[self] * solidAngle;
            else if ( i == 0 )
                table[self] = table[left] + img[self] * solidAngle;
            else if ( j == 0 )
                table[self] = table[up] + img[self] * solidAngle;
            else
                table[self] = table[left] + table[up] - table[up - 1] + img[self] * solidAngle;
        }
    }
    
    //int leftcorner[2] = {0, 0};
    int W = width - 1;
    int H = height - 1;
    bool isHorizontalCut = (H > W);
    float halfenergy = 0.5 * table[totalSize - 1];
    vector<Region> regions;
    Region reg = {{0, 0}, W, H, {0.5f, 0.5f}, isHorizontalCut, halfenergy};
    regions.push_back(reg);
    
    // Cut
    // int cutTimes = Log2(float(ns));
    int cutTimes = 6;
    float x_ratio = 1.f / float(width - 1);
    float y_ratio = 1.f / float(height - 1);
    for ( int time = 0 ; time < cutTimes ; time++ ){
        vector<Region> next_regions;
        for ( int i = 0 ; i < regions.size() ; i++ ){
            int move_index;
            int cut_x, cut_y;
            float E_left_up = 0.0f, E_left = 0.0f, E_up = 0.0f;
            
            if ( regions[i].leftCorner[0] == 0 && regions[i].leftCorner[1] != 0 )
                E_up = table[ regions[i].W + (regions[i].leftCorner[1] - 1) * width ];
            else if ( regions[i].leftCorner[0] != 0 && regions[i].leftCorner[1] == 0 )
                E_left = table[ regions[i].leftCorner[0] - 1 + regions[i].H * width ];
            else if ( regions[i].leftCorner[0] != 0 && regions[i].leftCorner[1] != 0 ){
                E_left_up   = table[ regions[i].leftCorner[0] - 1 + ( regions[i].leftCorner[1] - 1 ) * width ];
                E_left      = table[ regions[i].leftCorner[0] - 1 + ( regions[i].leftCorner[1] + regions[i].H ) * width ];
                E_up        = table[ regions[i].leftCorner[0] + regions[i].W + ( regions[i].leftCorner[1] - 1 ) * width ];
            }
            
            // Half Cut, Compute Energy
            if ( regions[i].horizontalCut ){    // --- Cut
                move_index = width;
                cut_x = regions[i].leftCorner[0] + regions[i].W;
                cut_y = regions[i].leftCorner[1] + 0.5 * regions[i].H;
                if ( regions[i].leftCorner[0] != 0 )
                    E_left = table[ regions[i].leftCorner[0] - 1 + cut_y * width ];
            }
            else{                               // ||| Cut
                move_index = 1;
                cut_x = regions[i].leftCorner[0] + 0.5 * regions[i].W;
                cut_y = regions[i].leftCorner[1] + regions[i].H;
                if ( regions[i].leftCorner[1] != 0 )
                    E_up = table[ cut_x + ( regions[i].leftCorner[1] - 1 ) * width ];
            }
            
            // move left or right, move up or down
            int cut_point = cut_x + cut_y * width;
            float cutEnergy = table[cut_point] - E_left - E_up + E_left_up;
            move_index *= ( cutEnergy > regions[i].halfEnergy ) ? -1 : 1;
            
            // Choose right cut
            while(true){
                int next_index = cut_point + move_index;
                if ( regions[i].horizontalCut && regions[i].leftCorner[0] != 0 )
                    E_left = table[ next_index - (regions[i].W + 1) ];
                else if ( !regions[i].horizontalCut && regions[i].leftCorner[1] != 0 )
                    E_up = table[ next_index - (regions[i].H + 1) * width ];
                
                float nextEnergy = table[next_index] - E_left - E_up + E_left_up;
                // This Point Cut or Not
                if ( (cutEnergy - regions[i].halfEnergy) * ( nextEnergy - regions[i].halfEnergy ) < 0 )
                    break;
                else{
                    cut_point = next_index;
                    cutEnergy = nextEnergy;
                }
            }
            
            // New cut regions
            cut_x = cut_point % width;
            cut_y = (cut_point - cut_x) / width;
            float halfenergy_1 = cutEnergy * 0.5;
            float halfenergy_2 = regions[i].halfEnergy - halfenergy_1;
            int leftcorner_1[2] = {regions[i].leftCorner[0], regions[i].leftCorner[1]};
            int leftcorner_2[2] = {0, 0};
            int W_1 = regions[i].W;
            int W_2 = regions[i].W;
            int H_1 = regions[i].H;
            int H_2 = regions[i].H;
            
            if ( regions[i].horizontalCut ){    // --- Cut
                H_1 = cut_y - regions[i].leftCorner[1];
                H_2 = regions[i].H - H_1 - 1;
                leftcorner_2[0] = regions[i].leftCorner[0];
                leftcorner_2[1] = cut_y + 1;
            }
            else{                               // ||| Cut
                W_1 = cut_x - regions[i].leftCorner[0];
                W_2 = regions[i].W - W_1 - 1;
                leftcorner_2[0] = cut_x + 1;
                leftcorner_2[1] = regions[i].leftCorner[1];
            }
            
            float center_1[2] = {(float(leftcorner_1[0]) + 0.5 * W_1) * x_ratio, (float(leftcorner_1[1]) + 0.5 * H_1) * y_ratio };
            float center_2[2] = {(float(leftcorner_2[0]) + 0.5 * W_2) * x_ratio, (float(leftcorner_2[1]) + 0.5 * H_2) * y_ratio };
            Region reg_1 = {{leftcorner_1[0], leftcorner_1[1]}, W_1, H_1, {center_1[0], center_1[1]}, (H_1 > W_1 * sinf(center_1[1] * M_PI)), halfenergy_1};
            Region reg_2 = {{leftcorner_2[0], leftcorner_2[1]}, W_2, H_2, {center_2[0], center_2[1]}, (H_2 > W_2 * sinf(center_2[1] * M_PI)), halfenergy_2};
            next_regions.push_back(reg_1);
            next_regions.push_back(reg_2);
        }
        regions = next_regions;
    }
    
    int numRegions = regions.size();
    PDF = 1.f / float(numRegions);
    for ( int i = 0 ; i < numRegions ; i++ ){
        RGBSpectrum tmp_spectrum = RGBSpectrum(0.f);
        for ( int y = 1 ; y <= regions[i].H ; y++ ){
            for ( int x = 1 ; x <= regions[i].W ; x++){
                tmp_spectrum += texels[ (regions[i].leftCorner[0] + x) + ( regions[i].leftCorner[1] + y ) * width ];
            }
        }
        regionSpectrum regS = {{regions[i].center[0], regions[i].center[1]}, tmp_spectrum};
        regionSpectrums.push_back(regS);
    }
        
    //delete[] texels;
    //delete[] table;
    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);
    delete[] img;
}

Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _InfiniteAreaLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _InfiniteAreaLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _InfiniteAreaLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _InfiniteAreaLight_ to SH from cube map sampling
        SHProjectCube(InfiniteAreaCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Find $(u,v)$ sample coordinates in infinite light texture
    // float uv[2], mapPdf;
    // distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    // if (mapPdf == 0.f) return 0.f;
    
    regionSpectrum regS = regionSpectrums[int(ls.uComponent * regionSpectrums.size())];
    
    // Convert infinite light sample point to direction
    float theta = regS.pos[1] * M_PI, phi = regS.pos[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));


    /* HERE */
    // Compute PDF for sampled infinite light direction
    // *pdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    // if (sintheta == 0.f) *pdf = 0.f;
    *pdf = PDF;
    Spectrum Ls = Spectrum(regS.spectrum, SPECTRUM_ILLUMINANT);
    // Return radiance value for infinite light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
    // Spectrum Ls = Spectrum(radianceMap->Lookup(uv[0], uv[1]),
    //                       SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
    PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Compute direction for infinite light sample ray

    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return Spectrum(0.f);

    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute _InfiniteAreaLight_ ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


