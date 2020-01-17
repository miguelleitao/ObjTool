/*
 *	intersection.c
 *
 *	Part of the ObjTool command line toolset.
 */


#define _GNU_SOURCE
 	
#include "ObjTool.h"

float Distance( Vert v1, Vert v2 ) {
    double total = 0.;
    for( int i=0 ; i<3 ; i++ ) {
        double d = v1[i]-v2[i];
        total += d*d;
    }
    return sqrt(total);
}

double findVerticalIntersection_master(ObjFile *obj, Vert coords) {
    // master version
    // under development
    int i;
    int zCoord = 3-findIntersection;
    printf("Finding vertical intersection %f %f %f\n", coords[0], coords[1], coords[2]);
    for( i=0 ; i<obj->stats.faces ; i++ ) {
        printf("  testing face %d with %d nodes\n", i, obj->faces[i].nodes  );
        float xMin =  1e8;
        float xMax = -1e8;
        float yMin =  1e8;
        float yMax = -1e8;
        for( int v=0 ; v<obj->faces[i].nodes ; v++ ) {
            printf("    testing node %d\n", v);
            int vidx = obj->faces[i].Node[v][0];
            printf("    testing node %d, vidx=%d\n", v, vidx);
            float vx = obj->verts[vidx][0];
            float vy = obj->verts[vidx][findIntersection];
            if ( vx<xMin ) xMin = vx;
            if ( vx>xMax ) xMax = vx;
            if ( vy<yMin ) yMin = vy;
            if ( vy>yMax ) yMax = vy;
        }
        if ( coords[0]<xMin || coords[0]>xMax ) continue;
        if ( coords[1]<yMin || coords[1]>yMax ) continue;
        // face intersected.
        
        printf("    face intersected !! \n");
        Vert medium = { 0., 0., 0. };
        double totDist = 0.;
        for( int v=0 ; v<obj->faces[i].nodes ; v++ ) {
            int vidx = obj->faces[i].Node[v][0];
            coords[zCoord] = obj->verts[vidx][zCoord];
            double dist = Distance(coords,obj->verts[vidx]);
            if ( dist<1e8 ) return obj->verts[vidx][zCoord];
            double invDist = 1. / dist;
            for( int i=0 ; i<3 ; i++ ) 
                medium[i] += (obj->verts[vidx][i]) * invDist;
            totDist += invDist;
        }
        for( i=0 ; i<3 ; i++ ) 
            medium[i] /= totDist;
        return medium[zCoord];
    }
    // No intersetion found.
    return 0.;
}

double findVerticalIntersection(ObjFile *obj, Vert coords) {
    // simplified version
    
    int zCoord = 3-findIntersection;
    int vidx = getVertexIdx(obj, 0, findIntersection, coords[0], coords[findIntersection]);
    if ( vidx<0 ) return 0.;
    
    return obj->verts[vidx][zCoord];
    
}
