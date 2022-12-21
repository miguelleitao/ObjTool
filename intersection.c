/*
 *	intersection.c
 *
 *	Part of the ObjTool command line toolset.
 */


#define _GNU_SOURCE
 	
#include <stdbool.h>
#include "ObjTool.h"

float Distance( Vert v1, Vert v2 ) {
    double total = 0.;
    for( int i=0 ; i<3 ; i++ ) {
        double d = v1[i]-v2[i];
        total += d*d;
    }
    return sqrt(total);
}

double signed_tetra_volume( Vert a, Vert b, Vert c, Vert d) {
    Vert t1, t2, t3, tc;
    vec3_subed(t1, b, a);
    vec3_subed(t2, c, a);
    vec3_subed(t3, d, a);
    vec3_cross(tc, t1, t2);
    return vec3_dot(tc, t3) / 6.;
}

int sign(double n) {
    if ( n>0. ) return  1;
    if ( n<0. ) return -1;
    return 0;
}

int tetra_volume_sign(Vert a, Vert b, Vert c, Vert d) {
    return sign( signed_tetra_volume( a, b, c, d) );
}
    
int intersect_line_triangle(Vert q1, Vert q2, Vert p1, Vert p2, Vert p3) {
    int s1 = tetra_volume_sign(q1,p1,p2,p3);
    int s2 = tetra_volume_sign(q2,p1,p2,p3);
    
    if ( s1 == s2 ) return false;
    
    int s3 = tetra_volume_sign(q1,q2,p1,p2);
    int s4 = tetra_volume_sign(q1,q2,p2,p3);
    int s5 = tetra_volume_sign(q1,q2,p3,p1);
    
    if ( s3!=s4 ) return false;
    if ( s4!=s5 ) return false;
    // Intersection
    // Intersection point (P) is:
    //    n = np.cross(p2-p1,p3-p1)
    //    t = np.dot(p1-q1,n) / np.dot(q2-q1,n)
    //    P = q1 + t * (q2-q1)
    return true;
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
        // This is an aproximation.
        // Intersetion test was performed only with the axis aligned bounding rectangle.
        // Works well if face is an axis aligned rectangle.
        
        
        
        
        // Face intersected !!
        printf("    face intersected !! \n");
        
        // Now zelect Z by the weighted average of Zi from all vertexes.
        // Zi are weigthed by the inverse horizontal distance.
        
        Vert medium = { 0., 0., 0. };	// Medium point.
        double totDist = 0.;
        for( int v=0 ; v<obj->faces[i].nodes ; v++ ) {
            int vidx = obj->faces[i].Node[v][0];
            coords[zCoord] = obj->verts[vidx][zCoord];
            double dist = Distance(coords, obj->verts[vidx]);
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

int getVertexIdx(ObjFile *obj, int xCoord, int yCoord, double x, double y) {
    int minI = -1;
    double minD = 1e8;
    
    if ( Verbose>2 )
        fprintf(stderr,"Looking for vertex %f %f ...", x,y);
    for( int i=0 ; i<obj->stats.verts ; i++ ) {
        double distX = obj->verts[i][xCoord] - x;
        double distY = obj->verts[i][yCoord] - y;
        double dist = distX*distX + distY*distY;
        if ( dist<minD ) {
            minD = dist;
            minI = i;
        }
    }
    if ( Verbose>2 ) {
        if (minI>=0) fprintf(stderr,"found vertex %d at dist %f\n", minI, minD);
        else fprintf(stderr,"erro\n");
    }
    return minI;
}

double findVerticalIntersection(ObjFile *obj, Vert coords) {
    // simplified version
    
    int zCoord = 3-findIntersection;
    int vidx = getVertexIdx(obj, 0, findIntersection, coords[0], coords[findIntersection]);
    if ( vidx<0 ) return 0.;
    printf("findVerticalIntersection retornando %f %f %f, z=%d, vidx=%d\n", obj->verts[vidx][0], obj->verts[vidx][1], obj->verts[vidx][2],
           zCoord, vidx);
    return obj->verts[vidx][zCoord];
    
}


