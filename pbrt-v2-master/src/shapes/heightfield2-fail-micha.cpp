
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


// shapes/heightfield.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"
#include <iostream>
#include <cmath>

using namespace std;

// Heightfield Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o, bool ro, int x, int y, const float *zs): Shape(o2w, w2o, ro) {
    
    nx = x;
    ny = y;
    nxMinus1 = nx-1.f;
    nyMinus1 = ny-1.f;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));
    
    numTriangles = 2*(nx-1)*(ny-1);
    numVertices = nx*ny;
    widthX = 1.f / (nx-1);
    widthY = 1.f / (ny-1);
    quadMaxZ = new float[(nx-1)*(ny-1)];
    quadMinZ = new float[(nx-1)*(ny-1)];
    vertexNormal = new Vector[nx * ny];
    P = new Point[nx * ny];
    uvs = new float[2 * nx * ny];
    triangleVertexIndex = new int[3 * numTriangles];
    
    computeNormal();
}


Heightfield2::~Heightfield2() {
    delete [] z;
    delete [] quadMaxZ;
    delete [] quadMinZ;
    delete [] vertexNormal;
    delete [] uvs;
    delete [] triangleVertexIndex;
}

void Heightfield2::computeNormal(){
    Vector *surfaceNormal = new Vector [numTriangles];
    
    // Vertex position
    int position = 0;
    int x;
    int y;
    for (y = 0; y < ny; ++y){
        for (x = 0; x < nx; ++x){
            P[position].x = (float)x * widthX;
            P[position].y = (float)y * widthY;
            P[position].z = z[position];
            uvs[2*position] = P[position].x;
            uvs[2*position + 1] = P[position].y;
            
            vertexNormal[position] = Vector(0, 0, 0);
            ++position;
        }
    }
    
//    position = 0;
//    for (y = 0; y < ny-1; ++y){
//        for (x = 0; x < nx-1; ++x){
//            quadMaxZ[position] = 0;
//            quadMinZ[position] = 0;
//            ++position;
//        }
//    }
    
    // Surface normal
    int *vp = triangleVertexIndex;
    int index[4];
    Vector point[4];
    int nowTriangleIndex = 0;
    int nowQuadIndex;
    for (y = 0; y < ny-1; ++y){
        for (x = 0; x < nx-1; ++x){
    #define VERT(x,y) ((x)+(y)*nx)
            *vp++ = index[0] = VERT(x, y);
            *vp++ = index[1] = VERT(x+1, y);
            *vp++ = index[2] = VERT(x+1, y+1);
            
            point[0] = Vector(P[index[0]]);
            point[1] = Vector(P[index[1]]);
            point[2] = Vector(P[index[2]]);
            
            surfaceNormal[nowTriangleIndex] = Normalize( Cross((point[0]-point[2]), (point[1]-point[2])) );
            
            vertexNormal[index[0]] += surfaceNormal[nowTriangleIndex];
            vertexNormal[index[1]] += surfaceNormal[nowTriangleIndex];
            vertexNormal[index[2]] += surfaceNormal[nowTriangleIndex];
            
            ++nowTriangleIndex;
            
            *vp++ = index[0] = VERT(x, y);
            *vp++ = index[1] = VERT(x+1, y+1);
            *vp++ = index[2] = VERT(x, y+1);
            
            point[0] = Vector(P[index[0]]);
            point[1] = Vector(P[index[1]]);
            point[2] = Vector(P[index[2]]);
            
            surfaceNormal[nowTriangleIndex] = Normalize( Cross((point[0]-point[2]), (point[1]-point[2])) );
            
            vertexNormal[index[0]] += surfaceNormal[nowTriangleIndex];
            vertexNormal[index[1]] += surfaceNormal[nowTriangleIndex];
            vertexNormal[index[2]] += surfaceNormal[nowTriangleIndex];
            
            ++nowTriangleIndex;
            
            index[0] = VERT(x, y);
            index[1] = VERT(x+1, y);
            index[2] = VERT(x+1, y+1);
            index[3] = VERT(x, y+1);
            
            point[0] = Vector(P[index[0]]);
            point[1] = Vector(P[index[1]]);
            point[2] = Vector(P[index[2]]);
            point[3] = Vector(P[index[3]]);
            
            nowQuadIndex = nowTriangleIndex/2 - 1;
            quadMaxZ[nowQuadIndex] = max(max(point[0].z, point[3].z), max(point[1].z, point[2].z) );
            quadMinZ[nowQuadIndex] = min(min(point[0].z, point[3].z), min(point[1].z, point[2].z) );
    #undef VERT
        }
    }
    
    position = 0;
    for (y = 0; y < ny; ++y){
        for (x = 0; x < nx; ++x){
            vertexNormal[position] = Normalize(vertexNormal[position]);
            ++position;
        }
    }
    
    delete [] surfaceNormal;
    return;
}

void Heightfield2::computeMaxMinZ(const Ray &ray, float rayEpsilon, int traceX, int traceY, float &minZ, float &maxZ) const{
    float x1 = traceX * widthX;
    float y1 = traceY * widthY;
    float x2 = (traceX + 1) * widthX;
    float y2 = (traceY + 1) * widthY;
    
    float raydx = ray.d.x;
    float raydy = ray.d.y;
    if ( raydx == 0 )
        raydx = rayEpsilon;
    if ( raydy == 0 )
        raydy = rayEpsilon;
    
    /* try to get middle */
    float tx1 = (x1 - ray.o.x) / raydx;
    float ty1 = (y1 - ray.o.y) / raydy;
    float tx2 = (x2 - ray.o.x) / raydx;
    float ty2 = (y2 - ray.o.y) / raydy;
    
    float max1 = max(tx1, ty1);
    float min1 = min(tx1, ty1);
    float max2 = max(tx2, ty2);
    float min2 = min(tx2, ty2);
    
    //cout << max1 << " / " << min1 << " / " << max2 << " / " << min2 << endl;
    
    float z1 = ray.o.z + min(max1, max2) * ray.d.z;
    float z2 = ray.o.z + max(min1, min2) * ray.d.z;
    maxZ = max(z1, z2);
    minZ = min(z1, z2);
    
    return;
}

bool Heightfield2::triangleIntersection(/*Output*/ DifferentialGeometry *dg, float *tHit, /*Input*/ const Ray &ray,  int *index) const {
    
    const Point &p1 = P[index[0]];
    const Point &p2 = P[index[1]];
    const Point &p3 = P[index[2]];
    
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;
    // Compute first barycentric coordinate
    Vector d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;
    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;
    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
    
    // Phong shading, Phong interpolation
    float tmp_uvs[3][2];
    tmp_uvs[0][0] = uvs[2 * index[0]];
    tmp_uvs[0][1] = uvs[2 * index[0] + 1];
    tmp_uvs[1][0] = uvs[2 * index[1]];
    tmp_uvs[1][1] = uvs[2 * index[1] + 1];
    tmp_uvs[2][0] = uvs[2 * index[2]];
    tmp_uvs[2][1] = uvs[2 * index[2] + 1];
    
    float b0 = 1.f - b1 - b2;
    float tu = b0 * tmp_uvs[0][0] + b1 * tmp_uvs[1][0] + b2 * tmp_uvs[2][0];
    float tv = b0 * tmp_uvs[0][1] + b1 * tmp_uvs[1][1] + b2 * tmp_uvs[2][1];
    
    //float tu = ray(ray.mint).x;
    //float tv = ray(ray.mint).y;
    
    Vector phongNormal = b0 * vertexNormal[index[0]] + b1 * vertexNormal[index[1]] + b2 * vertexNormal[index[2]];
    Vector tmpTangent(1, 0, 0);
    Vector tmpSurface(1, 0, 0);
    //phongNormal = Normalize(phongNormal);
    tmpSurface = Normalize( Cross(phongNormal, tmpTangent) );
    tmpTangent = Normalize( Cross(tmpSurface, phongNormal) );
    
    Vector dpdu = tmpTangent;
    Vector dpdv = tmpSurface;
    Normal dndu = Normal(1, 0, 0);
    Normal dndv = Normal(0, 0, 1);
    
    *dg = DifferentialGeometry( (*ObjectToWorld)(ray(t)),  (*ObjectToWorld)(dpdu), (*ObjectToWorld)(dpdv), (*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv), tu, tv, this);
    *tHit = t;
    
    return true;
}

bool Heightfield2::triangleIntersectionP(/*Input*/ const Ray &ray, int *index ) const {
    const Point &p1 = P[index[0]];
    const Point &p2 = P[index[1]];
    const Point &p3 = P[index[2]];
    
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;
    // Compute first barycentric coordinate
    Vector d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;
    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;
    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
    
    return true;
}

BBox Heightfield2::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool Heightfield2::CanIntersect() const {
    return true;
}

int sign(float a){
    if ( a == 0. )
        return 0;
    return (a > 0) ? 1: -1;
}

bool Heightfield2::Intersect(const Ray &ray, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const {
    
    //cout << "EPS" << endl;
    //cout << *rayEpsilon << endl;
    // Transform ray to object space
    Ray thisRay;
    (*WorldToObject)(ray, &thisRay);
    
    // Bounding box
    BBox BBoxHeightField2 = this->ObjectBound();
    float tBegin, tEnd;
    bool BBoxIntersect = BBoxHeightField2.IntersectP(thisRay, &tBegin, &tEnd);
    bool triangleIntersect = false;
    
    // Intesect or not
    if (!BBoxIntersect)
        return false;
    else{
        
        // First Hit
        Point firstHit = thisRay.o + tBegin * thisRay.d;
        int firstIndexX = firstHit.x * nxMinus1;
        int firstIndexY = firstHit.y * nyMinus1;
        
        /* check */
        if (firstIndexX < 0)
            firstIndexX = 0;
        else if (firstIndexX >= nxMinus1)
            firstIndexX = nxMinus1-1;   // nx-2
        
        if (firstIndexY < 0)
            firstIndexY = 0;
        else if (firstIndexY >= nyMinus1)
            firstIndexY = nyMinus1-1;   // ny-2
        
        // Final Hit
        //tEnd += *rayEpsilon;
        Point finalHit = thisRay.o + tEnd * thisRay.d;
        int finalIndexX = finalHit.x * nxMinus1;
        int finalIndexY = finalHit.y * nyMinus1;
        
//        cout << "---- final -----" << endl;
//        cout << finalIndexX << endl;
//        cout << finalIndexY << endl;
        
        /* check */
        if (finalIndexX < 0)
            finalIndexX = 0;
        else if (finalIndexX >= nxMinus1)
            finalIndexX = nxMinus1-1;   // nx-2
        
        if (finalIndexY < 0)
            finalIndexY = 0;
        else if (finalIndexY >= nyMinus1)
            finalIndexY = nyMinus1-1;   // ny-2
        
        float x1 = firstIndexX * widthX;
        float y1 = firstIndexY * widthY;
        float x2 = (firstIndexX + 1) * widthX;
        float y2 = (firstIndexY + 1) * widthY;
        
        /* vector does not effect */
        if ( thisRay.d.x == 0 )
            thisRay.d.x = *rayEpsilon;
        if ( thisRay.d.y == 0)
            thisRay.d.y = *rayEpsilon;
        
        float tx1 = ( x1 - thisRay.o.x ) / thisRay.d.x;
        float ty1 = ( y1 - thisRay.o.y ) / thisRay.d.y;
        float tx2 = ( x2 - thisRay.o.x ) / thisRay.d.x;
        float ty2 = ( y2 - thisRay.o.y ) / thisRay.d.y;
        float maxX = max(tx1, tx2);
        float maxY = max(ty1, ty2);
        
        /* begin trace ray */
        int traceX = firstIndexX;
        int traceY = firstIndexY;
        int index[3];
        
        float thisRaydx = thisRay.d.x;
        float thisRaydy = thisRay.d.y;
        float tDeltaX = widthX / thisRaydx;
        float tDeltaY = widthY / thisRaydy;
        
        while (1){
            
            if ( traceX < 0 || traceX >= nxMinus1 || traceY < 0 || traceY >= nyMinus1 )
                break;
            
            int position = traceY * nxMinus1 + traceX;
            int triangleIndex1 = 2 * position;
            int triangleIndex2 = 2 * position + 1;
            float maxZ, minZ;
            computeMaxMinZ(thisRay, *rayEpsilon, traceX, traceY, minZ, maxZ);
            
            /* not overflow or underflow */
            if ( !(minZ > (quadMaxZ[position]/* + *rayEpsilon*/)) && !(maxZ < (quadMinZ[position]/* - *rayEpsilon*/)) ){
            
                index[0] = triangleVertexIndex[triangleIndex1 * 3];
                index[1] = triangleVertexIndex[triangleIndex1 * 3 + 1];
                index[2] = triangleVertexIndex[triangleIndex1 * 3 + 2];
                triangleIntersect = triangleIntersection(dg, tHit, thisRay, index);
                
                if ( triangleIntersect )
                    break;
                
                index[0] = triangleVertexIndex[triangleIndex2 * 3];
                index[1] = triangleVertexIndex[triangleIndex2 * 3 + 1];
                index[2] = triangleVertexIndex[triangleIndex2 * 3 + 2];
                triangleIntersect = triangleIntersection(dg, tHit, thisRay, index);
                
                if ( triangleIntersect )
                    break;
                
                if ( traceX == finalIndexX && traceY == finalIndexY )
                    break;
            }
            
            if ( maxX < maxY ){
                maxX += abs(tDeltaX);
                traceX += sign(ray.d.x);
            }
            else {
                maxY += abs(tDeltaY);
                traceY += sign(ray.d.y);
            }
        }
    }
    
    return triangleIntersect;
}


bool Heightfield2::IntersectP(const Ray &ray) const {
    // Transform ray to object space
    Ray thisRay;
    (*WorldToObject)(ray, &thisRay);
    
    // Bounding box
    BBox BBoxHeightField2 = this->ObjectBound();
    float tBegin, tEnd;
    bool BBoxIntersect = BBoxHeightField2.IntersectP(thisRay, &tBegin, &tEnd);
    bool triangleIntersect = false;
    
    // Intesect or not
    if (!BBoxIntersect)
        return false;
    else{
        
        // First Hit
        Point firstHit = thisRay.o + tBegin * thisRay.d;
        int firstIndexX = firstHit.x * nxMinus1;
        int firstIndexY = firstHit.y * nyMinus1;
        
        /* check */
        if (firstIndexX < 0)
            firstIndexX = 0;
        else if (firstIndexX >= nxMinus1)
            firstIndexX = nxMinus1-1;   // nx-2
        
        if (firstIndexY < 0)
            firstIndexY = 0;
        else if (firstIndexY >= nyMinus1)
            firstIndexY = nyMinus1-1;   // ny-2
        
        
        // Final Hit
        //float EPS = 1e-3f * thisRay.mint;
        //tEnd += EPS;
        Point finalHit = thisRay.o + tEnd * thisRay.d;
        int finalIndexX = finalHit.x * nxMinus1;
        int finalIndexY = finalHit.y * nyMinus1;
        
        /* check */
        if (finalIndexX < 0)
            finalIndexX = 0;
        else if (finalIndexX >= nxMinus1)
            finalIndexX = nxMinus1-1;   // nx-2
        
        if (finalIndexY < 0)
            finalIndexY = 0;
        else if (finalIndexY >= nyMinus1)
            finalIndexY = nyMinus1-1;   // ny-2
        
        float x1 = firstIndexX * widthX;
        float y1 = firstIndexY * widthY;
        float x2 = (firstIndexX + 1) * widthX;
        float y2 = (firstIndexY + 1) * widthY;
        
        float tx1 = ( x1 - thisRay.o.x ) / thisRay.d.x;
        float ty1 = ( y1 - thisRay.o.y ) / thisRay.d.y;
        float tx2 = ( x2 - thisRay.o.x ) / thisRay.d.x;
        float ty2 = ( y2 - thisRay.o.y ) / thisRay.d.y;
        float maxX = max(tx1, tx2);
        float maxY = max(ty1, ty2);
        
        /* begin trace ray */
        int traceX = firstIndexX;
        int traceY = firstIndexY;
        int index[3];
        
        float thisRaydx = thisRay.d.x;
        float thisRaydy = thisRay.d.y;
        float tDeltaX = 1.f / (thisRaydx * nxMinus1);
        float tDeltaY = 1.f / (thisRaydy * nyMinus1);
        
        while (1){
            
            if ( traceX < 0 || traceX >= nxMinus1 || traceY < 0 || traceY >= nyMinus1 )
                break;
            
            int position = traceY * nxMinus1 + traceX;
            int triangleIndex1 = 2 * position;
            int triangleIndex2 = 2 * position + 1;
            float maxZ, minZ;
            computeMaxMinZ(thisRay, 0/*EPS*/, traceX, traceY, minZ, maxZ);
            
            /* not overflow or underflow */
            if ( !(minZ > (quadMaxZ[position]/* + EPS*/)) && !(maxZ < (quadMinZ[position]/* - EPS*/)) ){
                
                index[0] = triangleVertexIndex[triangleIndex1 * 3];
                index[1] = triangleVertexIndex[triangleIndex1 * 3 + 1];
                index[2] = triangleVertexIndex[triangleIndex1 * 3 + 2];
                triangleIntersect = triangleIntersectionP(thisRay, index);
                
                if ( triangleIntersect )
                    break;
                
                index[0] = triangleVertexIndex[triangleIndex2 * 3];
                index[1] = triangleVertexIndex[triangleIndex2 * 3 + 1];
                index[2] = triangleVertexIndex[triangleIndex2 * 3 + 2];
                triangleIntersect = triangleIntersectionP(thisRay, index);
                
                if ( triangleIntersect )
                    break;
                
                if ( traceX == finalIndexX && traceY == finalIndexY )
                    break;
            }
            
            if ( maxX < maxY ){
                maxX += abs(tDeltaX);
                traceX += sign(ray.d.x);
            }
            else {
                maxY += abs(tDeltaY);
                traceY += sign(ray.d.y);
            }
        }
    }
    
    return triangleIntersect;
}

void Heightfield2::GetShadingGeometry(const Transform &obj2world,
    const DifferentialGeometry &dg, DifferentialGeometry *dgShading) const{
    
    // Give Phong Normal
    * dgShading = dg;
    return ;
}


void Heightfield2::Refine(vector<Reference<Shape> > &refined) const {
    /* Default Not Use */
    /*int ntris = 2*(nx-1)*(ny-1);
    refined.reserve(ntris);
    int *verts = new int[3*ntris];
    Point *P = new Point[nx*ny];
    float *uvs = new float[2*nx*ny];
    int nverts = nx*ny;
    int x, y;
    // Compute heightfield vertex positions
    int pos = 0;
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
            P[pos].z = z[pos];
            ++pos;
        }
    }

    // Fill in heightfield vertex offset array
    int *vp = verts;
    for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y);
            *vp++ = VERT(x+1, y+1);
    
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y+1);
            *vp++ = VERT(x, y+1);
        }
#undef VERT
    }
    ParamSet paramSet;
    paramSet.AddInt("indices", verts, 3*ntris);
    paramSet.AddFloat("uv", uvs, 2 * nverts);
    paramSet.AddPoint("P", P, nverts);
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
    delete[] P;
    delete[] uvs;
    delete[] verts;*/
}


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}


