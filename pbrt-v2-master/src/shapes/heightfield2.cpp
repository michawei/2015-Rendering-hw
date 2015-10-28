
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
using namespace std;

#define DEFAULT 1
#define PHONG_INTERPOLATION 0

// Heightfield Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o, bool ro, int x, int y, const float *zs): Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    nxMinus1 = nx-1;
    nyMinus1 = ny-1;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));
    vertexNormal = new Normal[nx * ny];
    widthX = 1.f / (nx-1);
    widthY = 1.f / (ny-1);
    computeNormal();
}


Heightfield2::~Heightfield2() {
    delete [] z;
    delete [] vertexNormal;
}

void Heightfield2::computeNormal(){
    
    // Neighbor
    int neighbor1[12] = {-1, 0, -1, 1, 0, 1, 1, 0, 1, -1, 0, -1};
    int neighbor2[12] = {1, 0, 1, -1, 0, -1, -1, 0, -1, 1, 0, 1};
    int *v;
    int x, y;
    for (y = 0; y < ny; ++y){
        for (x = 0; x < nx; ++x){
            if (((x == nxMinus1) && (y > 0 && y <= nyMinus1)) || (((y == nyMinus1) && (x > 0 && x < nxMinus1))))
                v = neighbor2;
            else
                v = neighbor1;
            vector<Point> points;
            Normal result = Normal(0, 0, 0);
            // centor normal
            Point p = Point(x * widthX, y * widthY, z[x + nx * y]);
            
            int px, py, pz;
            for (int k = 0; k < 12; k += 2){
                px = x + v[k];
                py = y + v[k+1];
                pz = px + nx * py;
                if ((px >= 0) && (px <= (nxMinus1)) && (py >= 0) && (py <= (nyMinus1))){
                    Point point = Point(px * widthX, py * widthY, z[pz]);
                    points.push_back(point);
                }
            }
            vector<Point>::iterator it;
            Vector v1, v2;
            for (it = points.begin(); (it+1)!= points.end(); it++){
                v1 = *it - p;
                v2 = *(it+1) - p;
                result += Normalize(Normal(Cross(v1, v2)));
            }
            if (points.size() == 6) {
                it = points.begin();
                v2 = *it - p;
                v1 = *(it+5) - p;
                result += Normalize(Normal(Cross(v1, v2)));
            }
            vertexNormal[x + nx * y] = Normalize(result);
        }
    }
    
    // Surface normal - 2 triangle normal
    /*
     
    //Normal *surfaceNormal = new Normal [numTriangles];
    int *vp = triangleVertexIndex;
    int index[4];
    Vector point[4];
    int nowTriangleIndex = 0;
    //int nowQuadIndex;
    for (y = 0; y < ny-1; ++y){
        for (x = 0; x < nx-1; ++x){
    #define VERT(x,y) ((x)+(y)*nx)
            *vp++ = index[0] = VERT(x, y);
            *vp++ = index[1] = VERT(x+1, y);
            *vp++ = index[2] = VERT(x+1, y+1);
            
            point[0] = Vector(P[index[0]]);
            point[1] = Vector(P[index[1]]);
            point[2] = Vector(P[index[2]]);
            
            surfaceNormal[nowTriangleIndex] = Normalize(Normal(Cross((point[0]-point[2]), (point[1]-point[2]))));
            
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
            
            surfaceNormal[nowTriangleIndex] = Normalize(Normal(Cross((point[0]-point[2]), (point[1]-point[2]))));
            
            vertexNormal[index[0]] += surfaceNormal[nowTriangleIndex];
            vertexNormal[index[1]] += surfaceNormal[nowTriangleIndex];
            vertexNormal[index[2]] += surfaceNormal[nowTriangleIndex];
            
            ++nowTriangleIndex;
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
    */
    //delete [] surfaceNormal;
    return;
}


bool Heightfield2::triangleIntersection(/*Output*/ DifferentialGeometry *dg, float *tHit, float *rayEpsilon, /*Input*/ const Ray &ray,  Point *triangle) const {
    
    const Point &p1 = (*ObjectToWorld)(triangle[0]);
    const Point &p2 = (*ObjectToWorld)(triangle[1]);
    const Point &p3 = (*ObjectToWorld)(triangle[2]);
    
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
    
    
    // Compute triangle partial derivatives
    Vector dpdu, dpdv;
    
    // Compute deltas for triangle partial derivatives
    float du1 = triangle[0].x - triangle[2].x;
    float du2 = triangle[1].x - triangle[2].x;
    float dv1 = triangle[0].y - triangle[2].y;
    float dv2 = triangle[1].y - triangle[2].y;
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }
    
    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*triangle[0].x + b1*triangle[1].x + b2*triangle[2].x;
    float tv = b0*triangle[0].y + b1*triangle[1].y + b2*triangle[2].y;
    
    // Fill in _DifferentialGeometry_ from triangle hit
    *dg = DifferentialGeometry(ray(t), dpdu, dpdv,
                               Normal(0,0,0), Normal(0,0,0),
                               tu, tv, this);

    *tHit = t;
    *rayEpsilon = 1e-3f * *tHit;
    return true;
}

bool Heightfield2::triangleIntersectionP(/*Input*/ const Ray &ray, Point *triangle) const {
    const Point &p1 = (*ObjectToWorld)(triangle[0]);
    const Point &p2 = (*ObjectToWorld)(triangle[1]);
    const Point &p3 = (*ObjectToWorld)(triangle[2]);
    
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
#if DEFAULT
    return false;
#else
    return true;
#endif
}

int sign(float a){
    return (a >= 0) ? 1: -1;
}

bool Heightfield2::Intersect(const Ray &ray, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const {
    
    //cout << *rayEpsilon << endl;
    // Transform ray to object space
    Ray thisRay;
    (*WorldToObject)(ray, &thisRay);
    
    // Bounding box
    BBox BBoxHeightField2 = this->ObjectBound();
    float tMin;
    if (BBoxHeightField2.Inside(thisRay(thisRay.mint)))
        tMin = thisRay.mint;
    else if (!BBoxHeightField2.IntersectP(thisRay, &tMin))
        return false;
    
    // First Hit
    Point firstHit = thisRay(tMin);
    
    /* check */
    int indexX, indexY;
    float tNextX, tNextY;
    float tDeltaX, tDeltaY;
    int signX, signY;
    int outX, outY;
    
    // DDA set
    indexX = Clamp((int)(firstHit.x * nxMinus1), (int)0, (int)(nxMinus1-1));
    signX = sign(thisRay.d.x);
    if ( thisRay.d.x >= 0 ) {
        tNextX = tMin + ( (indexX + 1) * widthX - firstHit.x ) / thisRay.d.x;
        tDeltaX = widthX / thisRay.d.x;
        outX = nxMinus1;
    }
    else {
        tNextX = tMin + ( (indexX) * widthX - firstHit.x ) / thisRay.d.x;
        tDeltaX = (-1) * widthX / thisRay.d.x;
        outX = -1;
    }
    
    indexY = Clamp((int)(firstHit.y * nyMinus1), (int)0, (int)(nyMinus1-1));
    signY = sign(thisRay.d.y);
    if ( thisRay.d.y >= 0 ) {
        tNextY = tMin + ( (indexY + 1) * widthY - firstHit.y ) / thisRay.d.y;
        tDeltaY = widthY / thisRay.d.y;
        outY = nyMinus1;
    }
    else {
        tNextY = tMin + ( (indexY) * widthY - firstHit.y ) / thisRay.d.y;
        tDeltaY = (-1) * widthY / thisRay.d.y;
        outY = -1;
    }
    
    // run on DDA line
    Point triangle[3];
    while (1){
        
        triangle[0] = Point(indexX * widthX, indexY * widthY, z[indexX + indexY * nx]);
        triangle[1] = Point((indexX + 1) * widthX, indexY * widthY, z[(indexX + 1) + indexY * nx]);
        triangle[2] = Point(indexX * widthX, (indexY + 1) * widthY, z[indexX + (indexY + 1) * nx]);
        
        if ( triangleIntersection(dg, tHit, rayEpsilon, ray, triangle) )
            return true;
        
        triangle[0] = Point((indexX + 1) * widthX, (indexY + 1) * widthY, z[(indexX + 1) + (indexY + 1) * nx]);
        
        if ( triangleIntersection(dg, tHit, rayEpsilon, ray, triangle) )
            return true;
        
        if ( tNextX < tNextY ){
            if ( thisRay.maxt < tNextX )
                break;
            indexX += signX;
            if ( indexX == outX )
                break;
            tNextX += tDeltaX;
        }
        else{
            if ( thisRay.maxt < tNextY )
                break;
            indexY += signY;
            if ( indexY == outY )
                break;
            tNextY += tDeltaY;
        }
    }
    return false;
}


bool Heightfield2::IntersectP(const Ray &ray) const {
    // Transform ray to object space
    Ray thisRay;
    (*WorldToObject)(ray, &thisRay);
    
    // Bounding box
    BBox BBoxHeightField2 = this->ObjectBound();
    float tMin;
    if (BBoxHeightField2.Inside(thisRay(thisRay.mint)))
        tMin = thisRay.mint;
    else if (!BBoxHeightField2.IntersectP(thisRay, &tMin))
        return false;
    
    // First Hit
    Point firstHit = thisRay(tMin);
    
    /* check */
    int indexX, indexY;
    float tNextX, tNextY;
    float tDeltaX, tDeltaY;
    int signX, signY;
    int outX, outY;
    
    // DDA set
    indexX = Clamp((int)(firstHit.x * nxMinus1), (int)0, (int)(nxMinus1-1));
    signX = sign(thisRay.d.x);
    if ( thisRay.d.x >= 0 ) {
        tNextX = tMin + ( (indexX+1) * widthX - firstHit.x ) / thisRay.d.x;
        tDeltaX = widthX / thisRay.d.x;
        outX = nxMinus1;
    }
    else {
        tNextX = tMin + ( (indexX) * widthX - firstHit.x ) / thisRay.d.x;
        tDeltaX = (-1) * widthX / thisRay.d.x;
        outX = -1;
    }
    
    indexY = Clamp((int)(firstHit.y * nyMinus1), (int)0, (int)(nyMinus1-1));
    signY = sign(thisRay.d.y);
    if ( thisRay.d.y >= 0 ) {
        tNextY = tMin + ( (indexY+1) * widthY - firstHit.y ) / thisRay.d.y;
        tDeltaY = widthY / thisRay.d.y;
        outY = nyMinus1;
    }
    else {
        tNextY = tMin + ( (indexY) * widthY - firstHit.y ) / thisRay.d.y;
        tDeltaY = (-1) * widthY / thisRay.d.y;
        outY = -1;
    }
    
    // run on DDA line
    Point triangle[3];
    while (1){
        
        triangle[0] = Point(indexX * widthX, indexY * widthY, z[indexX + indexY * nx]);
        triangle[1] = Point((indexX + 1) * widthX, indexY * widthY, z[(indexX + 1) + indexY * nx]);
        triangle[2] = Point(indexX * widthX, (indexY + 1) * widthY, z[indexX + (indexY + 1) * nx]);
        
        if ( triangleIntersectionP(ray, triangle) )
            return true;
        
        triangle[0] = Point((indexX + 1) * widthX, (indexY + 1) * widthY, z[(indexX + 1) + (indexY + 1)*nx]);
        
        if ( triangleIntersectionP(ray, triangle) )
            return true;
        
        if ( tNextX < tNextY ){
            if ( thisRay.maxt < tNextX )
                break;
            indexX += signX;
            if ( indexX == outX )
                break;
            tNextX += tDeltaX;
        }
        else{
            if ( thisRay.maxt < tNextY )
                break;
            indexY += signY;
            if ( indexY == outY )
                break;
            tNextY += tDeltaY;
        }
    }
    
    return false;
}

void Heightfield2::GetShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry *dgShading) const{
    
#if PHONG_INTERPOLATION
    // Phong shading
    // find hit triangle
    // Hit point index
    int x = Clamp((int)(dg.u * nxMinus1), (int)0, (int)nxMinus1);
    int y = Clamp((int)(dg.v * nyMinus1), (int)0, (int)nyMinus1);
    
    Point pointXAdd1 = Point((x + 1) * widthX, y * widthY, z[y * nx + x + 1]);
    Point pointYAdd1 = Point(x * widthX, (y + 1) * widthY, z[(y + 1) * nx + x]);
    Point pointDiagonal = Point(dg.u, pointYAdd1.y - ( dg.u - pointYAdd1.x ), 0);
    Point pointOriginal;
    Normal n[3];
    
    // triangle[0] is (x, y) or (x+1, y+1), pick point0 and their normal
    if ( dg.v <= pointDiagonal.y ) {
        pointOriginal = Point(x * widthX, y * widthY, z[y * nx + x]);
        n[0] = Normal(vertexNormal[y * nx + x]);
        n[1] = Normal(vertexNormal[y * nx + (x + 1)]);
        n[2] = Normal(vertexNormal[(y + 1) * nx + x]);
    }
    else {
        pointOriginal = Point((x + 1) * widthX, (y + 1) * widthY, z[(y + 1) * nx + x + 1]);
        n[0] = Normal(vertexNormal[(y + 1) * nx + (x + 1)]);
        n[1] = Normal(vertexNormal[y * nx + (x + 1)]);
        n[2] = Normal(vertexNormal[(y + 1) * nx + x]);
    }

    // Phong Interpolation begin
    // set z = 0 for use "Length()" function
    pointXAdd1.z = 0;
    pointYAdd1.z = 0;
    
    Normal n1 = Normalize((pointDiagonal - pointXAdd1).Length() * n[2] + (pointDiagonal - pointYAdd1).Length() * n[1]);
    Normal n2;

    if ( dg.v < pointDiagonal.y )
        n2 = Normalize(fabs(dg.u - pointOriginal.x) * n[1] + fabs(pointXAdd1.x - dg.u) * n[0]);
    else
        n2 = Normalize(fabs(dg.u - pointYAdd1.x) * n[0] + fabs(pointOriginal.x - dg.u) * n[2]);
    
    Normal hitNormal = (*ObjectToWorld)(Normalize(fabs(dg.v - pointOriginal.y) * n1 + fabs(pointDiagonal.y - dg.v) * n2));
				
    // triangle mesh dndu, dndv
    Normal dndu, dndv;
    float du1 = pointOriginal.x - pointYAdd1.x;
    float du2 = pointXAdd1.x    - pointYAdd1.x;
    float dv1 = pointOriginal.y - pointYAdd1.y;
    float dv2 = pointXAdd1.y    - pointYAdd1.y;
    Normal dn1 = n[0] - n[2];
    Normal dn2 = n[1] - n[2];
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        dndu = dndv = Normal(0, 0, 0);
    }
    else {
        float invdet = 1.f / determinant;
        dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
        dndv = (-du2 * dn1 + du1 * dn2) * invdet;
    }

    Vector ss = Normalize(dg.dpdu);
    Vector ts = Cross(ss, hitNormal);
    if (ts.LengthSquared() > 0.f) {
        ts = Normalize(ts);
        ss = Cross(ts, hitNormal);
    }
    else
        CoordinateSystem((Vector)hitNormal, &ss, &ts);
    
    *dgShading = DifferentialGeometry(dg.p, ss, ts, (*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv), dg.u, dg.v, dg.shape);
    dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
    dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
    dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
#else
    *dgShading = dg;
#endif
    return;
}


void Heightfield2::Refine(vector<Reference<Shape> > &refined) const {
    /* Default Not Use */
    int ntris = 2*(nx-1)*(ny-1);
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
    delete[] verts;
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


