/*
 *  ObjTool.c
 *
 *  Command line toolset to work with Wavefront 3D Obj files.
 *
 *  Uses OBJ file format as defined in http://paulbourke.net/dataformats/obj/
 *  Not all OBJ file format specs are implemented.
 */

#include "ObjTool.h"

// Globals
Vert VMin = { -1e30, -1e30, -1e30 };
Vert VMax = {  1e30,  1e30,  1e30 };
Vert Translate = { 0., 0., 0. };
Vert Scale = { 1., 1., 1. };
Vert Rotate = { 0., 0., 0. };
Vert Vertical = { 0., 0., 0. };
char *Material = NULL;
char *OutputFile = NULL;
char *ShadowOutputFile = NULL;
short int SolidCut = 0;
short int Relative = 0;
short int Verbose = 0;
short int info = 0;
short int findIntersection = 0;
short int Negate = 0;
short int Explode = 0;
short int InvertNormals = 0;
char *SelectGroup[20];
char *SelectObject[20];
int SelGroups = 0;
int SelObjects = 0;
short int genTextCoords = 0;
char upAxis = 'y';  // Used for shadow calculation.
                    // Possible values: {y,z}

#define SinF(x) ( sinf(x) )
#define CosF(x) ( cosf(x) )

#define AngNormalize(x) ( (x)>PI ? (x-2*PI) : ( (x)<-PI ? (x+2*PI) : (x) ) )

int JoinObjFiles(int nObjs, ObjFile ObjSet[], ObjFile *obj);
void SetUseCounters(ObjFile *obj);
void SetIndexs(ObjFile *obj);
int getVertexIdx(ObjFile *obj, int xCoord, int yCoord, double x, double y);

char* strrstr(char *Str, const char *SubStr)
{
    int i=0;
    int find_len = strlen(SubStr);

    for( i=strlen(Str)-find_len; i>=0; i--) {
        if( ! strncmp(Str+i, SubStr, find_len) )
                return Str+i;
    }
    return NULL;    
}

int VertIndex(ObjFile *obj, int i, int total) {
    if ( Relative )
    return -total+obj->order.verts[i-1] ;
    else
    return 1+obj->order.verts[i-1];
}

int TextIndex(ObjFile *obj, int i, int total) {
    if ( Relative )
    return -total+obj->order.texts[i-1];
    else
    return 1+obj->order.texts[i-1];
}

int NormIndex(ObjFile *obj, int i, int total) {
    if ( Relative )
    return -total+obj->order.norms[i-1];
    else
    return 1+obj->order.norms[i-1];
}

void PrintObjStats(ObjStats *os) {
    printf("# Obj Stats\n");
    printf("  Vertices: %d\n", os->verts);
    if ( os->verts > 0 ) {
        printf("    Min: %.3f %.3f %.3f\n", os->vmin[0], os->vmin[1], os->vmin[2]);
        printf("    Max: %.3f %.3f %.3f\n", os->vmax[0], os->vmax[1], os->vmax[2]);
    }
    printf("  Normals: %d\n", os->norms);
    printf("  Texture Coords: %d\n", os->texts);
    if ( os->texts > 0 ) {
        printf("    Min: %.3f %.3f\n", os->tmin[0], os->tmin[1]);
        printf("    Max: %.3f %.3f\n", os->tmax[0], os->tmax[1]);
    }
    printf("  Faces: %d\n", os->faces);
    printf("  Objects: %d\n", os->objs);
    printf("  ObjGroups: %d\n", os->grps);
    printf("  UseMtl: %d\n", os->mats);
    printf("  MtlLib:   %d\n", os->libs);
    printf("  ShadGrp:  %d\n", os->shds);
    printf("#\n");
}

void PrintLongObjStats(ObjFile *obj) {
        PrintObjStats(&(obj->stats));
    int i;
    for( i=0 ; i<obj->stats.mats ; i++ )
            printf("use mat[%d]: %d : %s\n", i+1, obj->mats[i].line, obj->mats[i].name);
    for( i=0 ; i<obj->stats.objs ; i++ )
            printf("objs[%d]: %d : %s\n", i+1, obj->objs[i].line, obj->objs[i].name);
    for( i=0 ; i<obj->stats.grps ; i++ )
            printf("grps[%d]: %d : %s\n", i+1, obj->grps[i].line, obj->grps[i].name);           
}

void PrintFullObjStats(ObjFile *obj) {
    PrintObjStats(&(obj->stats));
    int i;
    for( i=0 ; i<obj->stats.verts ; i++ )
        printf("  vert[%d]: %d\n", i+1, obj->counts.verts[i]);
    for( i=0 ; i<obj->stats.texts ; i++ )
        printf("  text[%d]: %d\n", i+1, obj->counts.texts[i]);
    for( i=0 ; i<obj->stats.norms ; i++ )
        printf("  norm[%d]: %d\n", i+1, obj->counts.norms[i]);
    for( i=0 ; i<obj->stats.mats ; i++ )
        printf("  use mat[%d]: %d - %s\n", i+1, obj->mats[i].line, obj->mats[i].name);
}



int VerticalCrossSegment(Vert *p, Vert v1, Vert v2) {
    // Test if vector UP from p intersect segment v1,v2
    // p,v1,v2 are 2d points

    double dx = v2[0] - v1[0];
    if ( fabs(dx)<0.00001 )     // Segment [v1,v2] is vertical
        return 0;
    int res = 0;
    double dy = v2[1] - v1[1];
    double m = dy / dx;     
    double y0 = v1[1]-m*(v1[0]-(*p)[0]);
    if ( y0>(*p)[1] ) {
        if ( dy>0 ) {
            if ( y0-v1[1]> .1 && y0-v2[1]<-.1 ) res = 1;
        }
        else {
            if ( y0-v1[1]<-.1 && y0-v2[1]> .1 ) res = 1;
        }
    }
    //printf( "        v1 %f %f, v2 %f %f, res:%d\n", v1[0],v1[1],v2[0],v2[1],res );
    return res; 
}

int PointInPolygon2D(Vert *p, Vert *verts, FaceNode *node, int n) {
    // Test if vertex p is outside 2D polygon verts[node[]]
    int cross = 0; 
    if ( n<3 ) return cross;    // polygon must be closed   
    int i;
    for( i=1 ; i<n ; i++ )
        cross += VerticalCrossSegment(p,verts[node[i-1][0]],verts[node[i][0]]);
    cross += VerticalCrossSegment(p,verts[node[0][0]],verts[node[n-1][0]]);
//  printf("    PointInPoly n:%d vert %d, point p %f %f: cross: %d\n", n, v, (*p)[0], (*p)[1], cross);
    return cross % 2;
}

int VertexInPolygon2D(int v, Vert *verts, FaceNode *node, int n) {
    // Test if vertex verts[v] is outside 2D polygon verts[node[]]
    Vert *p = verts+v;
    return PointInPolygon2D(p, verts, node, n);
}

/*
int VertexInPolygon2D(int v, Vert *verts, FaceNode *node, int n) {
    // Test if vertex verts[v] is outside 2D polygon verts[node[]]
    Vert *p = verts+v;
    int cross = 0; 
    if ( n<3 ) return cross;    // polygon must be closed   
    int i;
    for( i=1 ; i<n ; i++ )
        cross += VerticalCrossSegment(p,verts[node[i-1][0]],verts[node[i][0]]);
    cross += VerticalCrossSegment(p,verts[node[0][0]],verts[node[n-1][0]]);
//  printf("    PointInPoly n:%d vert %d, point p %f %f: cross: %d\n", n, v, (*p)[0], (*p)[1], cross);
    return cross % 2;
}
*/

double SegmentsIntersection(Vec2 p1, Vec2 p2, Vec2 q1, Vec2 q2) {
    // Test if 2D line segment [p1,p2] intersects 2D line segment [q1,q2]
    vec2 pmin, pmax, qmin, qmax;
    int i;

    //printf("interseting (%f,%f) (%f,%f) - (%f,%f) (%f,%f)\n",p1[0],p1[1],p2[0],p2[1],q1[0],q1[1],q2[0],q2[1]);
    for( i=0 ; i<2 ; i++ ) {
        if ( p1[i]<p2[i] ) {
            pmin[i] = p1[i];
            pmax[i] = p2[i];
        } else {
            pmin[i] = p2[i];
            pmax[i] = p1[i];
        } 
        if ( q1[i]<q2[i] ) {
            qmin[i] = q1[i];
            qmax[i] = q2[i];
        } else {
            qmin[i] = q2[i];
            qmax[i] = q1[i];
        }
        if ( pmax[i]<qmin[i] ) return -1.;  // Bounding boxes non-overlapping
        if ( qmax[i]<pmin[i] ) return -1.;
    }
    //printf("bb interseted\n");       
    vec2 p, q;
    vec2_subed(p, p2, p1);
    vec2_subed(q, q2, q1);
    
    vec2 d1;
    vec2_subed(d1, q1, p1);
    double u0 = vec2_cross(d1, p);
    double pq = vec2_cross(p, q);
    
    //printf("     u0=%f, pq=%f\n",u0,pq);
    if ( abs(pq)<1e-10 ) {          //   parallel segments
        if ( abs(u0)<1e-10 )        //   collinear segments
            return 0;
        return -1.;                 //   non-intersecting
    }
    
    double t = vec2_cross(d1, q)/pq;
    double u = u0/pq;
    if ( t>=0 && t<=1. && u>=0 && u<=1. )            // intersecting
        return t;
    return -1;                      // lines intersecting but outside segments.
    
}

double FindClosestIntersection(float x1, float y1,
                               float x2, float y2, 
                               ObjFile *obj, 
                               int *nfi, int *nni) {
    int ycoord = 1; // must be the same as the one defined at CreateShadowObj
    
    if ( upAxis=='y' )  ycoord = 2;
    
    // finds edge from obj that crosses segment (x1,y1)-(x2,y2) nearest to (x1,y1)
    vec2 p1 = { x1, y1 };
    vec2 p2 = { x2, y2 };
    double dist = MAX_FLOAT;
    int fi, ni;
    //double dy = y2-y1;
    //double dx = x2-x1;
    for( fi=0 ; fi<obj->stats.faces ; fi++ ) {
        for( ni=0 ; ni<obj->faces[fi].nodes ; ni++ ) {
            vec2 q1, q2;
            int vi = obj->faces[fi].Node[ni][0]-1;
            q1[0] = obj->verts[vi][0];
            q1[1] = obj->verts[vi][ycoord];
            int next = ni+1;
            if ( next>=obj->faces[fi].nodes ) next = 0;
            vi = obj->faces[fi].Node[next][0]-1;
            q2[0] = obj->verts[vi][0];
            q2[1] = obj->verts[vi][ycoord];
            double idist = SegmentsIntersection(p1,p2,q1,q2);
            if ( idist>0 && idist<dist ) {
                printf(" intersected, d=%f\n",idist);
                dist = idist;
                *nfi = fi;
                *nni = ni;
            }
        }
    }
    if ( dist>MAX_FLOAT/2. ) return -1.;    // No intersection found
    return dist;
}

ObjFile *CreateShadowObj_v1(ObjFile *obj) {
  //  printf("Creating shadow\n");
    int i, si;  // face indexs
/*
    double max[2];
    double min[2];
    for( i=0 ; i<2 ; i++ ) {
    max[i] = -MAX_FLOAT;
    min[i] =  MAX_FLOAT;
    }
*/
    ObjFile *shadow = Malloc(1,sizeof(ObjFile));

    shadow->stats.verts = obj->stats.verts;
    shadow->verts = Malloc(obj->stats.verts, sizeof(Vert));
    memcpy(shadow->verts, obj->verts, obj->stats.verts * sizeof(Vert));

    shadow->faces = Malloc(obj->stats.faces, sizeof(Face));
    shadow->stats.faces = 0;
  
    shadow->stats.norms = 1;    // Only one norm, facing up.

    shadow->stats.texts = obj->stats.texts;

    shadow->norms = Malloc(1, sizeof(Norm));
    shadow->norms[0][0] = 0.;
    shadow->norms[0][1] = 0.;
    shadow->norms[0][2] = 1.;   // Up vector

    shadow->texts = obj->texts;
    shadow->grps = obj->grps ;
    shadow->mats = obj->mats;
    shadow->objs = obj->objs;
    shadow->libs = obj->libs;
    shadow->shds = obj->shds;

    shadow->counts.verts = obj->counts.verts;
    shadow->counts.norms = obj->counts.norms;
    shadow->counts.texts = obj->counts.texts;

    shadow->order.verts = obj->order.verts;
    shadow->order.texts = obj->order.texts;
    shadow->order.norms = obj->order.norms ;

    // Find z min
    float zmin = MAX_FLOAT;
    for( i=0 ; i<shadow->stats.verts ; i++ )
    if ( shadow->verts[i][2]<zmin )
        zmin = shadow->verts[i][2];
    // smachdown zz coords to zmin
    for( i=0 ; i<shadow->stats.verts ; i++ )
    shadow->verts[i][2] = zmin;

    for( i=0 ; i<obj->stats.faces ; i++ ) {
        int n;
        short int inside = 0;
//printf("Testing face %d, %d nodes\n", i, obj->faces[i].nodes);
        for( n=0 ; n<obj->faces[i].nodes ; n++ ) {
        //printf(" Testing node %d\n",n);
        int v = obj->faces[i].Node[n][0];
        if ( v<1 || v>obj->stats.verts ) {
            fprintf(stderr, "Vertice %d not found on face %d, node %d/%d.\n", v, i, n, obj->faces[i].nodes);
            break;
        }
        // Test if vertex is outside polygon already draw
        int pp;
        for( pp=0 ; pp<shadow->stats.faces ; pp++ ) {
            if ( VertexInPolygon2D( v-1, shadow->verts, shadow->faces[pp].Node, shadow->faces[pp].nodes ) ) {
                inside += 1;
                break;
            }
        }
        }
        if ( inside<=n ) {  // Add face
        //printf("Adding shadow face\n");
        // Reusing (pointing to) nodes from original obj
        si = shadow->stats.faces;
        memcpy( &(shadow->faces[si]), &(obj->faces[i]), sizeof(Face) );
/*
        shadow->faces[si].nodes = obj->faces[i].nodes; */
        for( n=0 ; n<shadow->faces[si].nodes ; n++ )
            shadow->faces[si].Node[n][2] = 1;
    
        shadow->stats.faces += 1;
        }
        else printf("### FACE NOT ADDED\n");
    }       
    return shadow;
}

ObjFile *CreateShadowObj(ObjFile *obj) {
    // New version using convex hull
    // Only works in y or z direction.
    // Result is always a convex poligon.
    // ToDo: Support for non-convex shadow poligons.
    //       This will require working with edges instead of vertexes. 
    
    int ycoord = 1;
    int zcoord = 2;
    
    if ( upAxis=='y' )  {
        ycoord = 2;
        zcoord = 1;
    }
    if (Verbose>2) fprintf(stderr,"Creating shadow\n");
  
    int i, si;  // face indexs
    ObjFile *shadow = Malloc(1, sizeof(ObjFile));

    shadow->stats.verts = 0;
    shadow->verts = Malloc(obj->stats.verts, sizeof(Vert));

    shadow->stats.faces = 0;

    // Add one normal
    shadow->stats.norms = 1;    // Only one norm, facing up.
    shadow->norms = Malloc(shadow->stats.norms, sizeof(Norm));
    shadow->norms[0][0] = 0.;
    shadow->norms[0][ycoord] = 0.;
    shadow->norms[0][zcoord] = 1.;  // Up vector

    // Add two texture coords
    shadow->stats.texts = 2;
    shadow->texts = Malloc(shadow->stats.texts, sizeof(Text));
    shadow->texts[0][0] = 0.01;         // Dark texture point
    shadow->texts[0][1] = 0.01;
    shadow->texts[1][0] = 0.91;         // Light texture point
    shadow->texts[1][1] = 0.91;
        
    // Add one MtlLib
    shadow->stats.libs = 1;
    shadow->libs = Malloc(shadow->stats.libs, sizeof(ObjGroup));
    shadow->libs[0].line = 0;       // before all faces
    shadow->libs[0].name = "shadow.mtl";
    
    // Add one Material
    shadow->stats.mats = 1;
    shadow->mats = Malloc(shadow->stats.mats, sizeof(ObjGroup));
    shadow->mats[0].line = 0;       // before all faces
    shadow->mats[0].name = "shadow";
    
    shadow->grps = NULL;
    shadow->objs = NULL;
    shadow->shds = NULL;

    shadow->counts.verts = NULL;
    shadow->counts.norms = NULL;
    shadow->counts.texts = NULL;

    shadow->order.verts = NULL;
    shadow->order.texts = NULL;
    shadow->order.norms = NULL;
        
   // Find z min
   // Shadow is flat poligon at Zmin
   float zmin = MAX_FLOAT;
   for( i=0 ; i<obj->stats.verts ; i++ )
       if ( obj->verts[i][zcoord]<zmin )
           zmin = obj->verts[i][zcoord];
   if ( Verbose>2 ) printf("  ShadowZ: %f\n", zmin);

   // Find leftest vertex
   float xmin = MAX_FLOAT;
   float ymin = MAX_FLOAT;
   for( i=0 ; i<obj->stats.verts ; i++ )
       if ( obj->verts[i][0]<xmin ) {
           xmin = obj->verts[i][0];
           ymin = obj->verts[i][ycoord];
       }
   if ( Verbose>2 ) printf("  Shadow Min(x,y): %f, %f\n", xmin, ymin);
   
   // Find surrounding convex polygon
   float xmax, ymax;
   float sum_x = 0.;
   float sum_y = 0.;
   // Star at leftest point.
   float xp = xmin; // xx coord at current vertex (cursor).
   float yp = ymin; // yy coord at current vertex (cursor).
   // Looking from left to right
   while (1) {
       double mmax = -MAX_FLOAT;   // max slope
       ymax = xmax = -MAX_FLOAT;
       // Find next vertex, with greatest slope
       for( i=0 ; i<obj->stats.verts ; i++ ) {
           float x = obj->verts[i][0];
           float y = obj->verts[i][ycoord];
           if ( x<xp ) continue;
               if ( x==xp ) {      // Vertical
                   if ( y<=yp ) continue;
                   mmax = MAX_FLOAT;
                   xmax = x;
                   ymax = y;
                   break;
               }
           double m = (y-yp)/(x-xp);   // slope
           if ( m>mmax ) {
               mmax = m;
               xmax = x;
               ymax = y;
           }
       }
       if ( Verbose>3 ) printf("    mmax: %f, x,y: %f, %f\n", mmax, xmax, ymax);
       if ( mmax<=-MAX_FLOAT/2. ) {
           if ( Verbose>3 ) printf("    No next point. Rightest point reched\n");
           break;
       }
       // new border vertex found at (xmax,ymax)
       xp = xmax;
       yp = ymax;
       sum_x += xmax;
       sum_y += ymax;
       si = shadow->stats.verts;
       shadow->verts[si][0] = xp;
       shadow->verts[si][ycoord] = yp;
       shadow->stats.verts = si+1;
   }
   // Looking from right to left
   if ( Verbose>2 ) printf("  Shadow Max(x,y): %f, %f\n", xp, yp);
   while (1) {
       double mmax = -MAX_FLOAT;
       ymax = xmax = -MAX_FLOAT;
       // Find next vertex, with largest slope
       for( i=0 ; i<obj->stats.verts ; i++ ) {
           float x = obj->verts[i][0];
           float y = obj->verts[i][ycoord];
           if ( x>xp ) continue;
           if ( x==xp ) {      // Vertical
               if ( y>=yp ) continue;
               mmax = MAX_FLOAT;
               xmax = x;
               ymax = y;
               break;
           }
           double m = (y-yp)/(x-xp);
           if ( m>mmax ) {
               mmax = m;
               xmax = x;
               ymax = y;
           }
       }
       if ( Verbose>3 ) printf("    mmax: %f, x,y: %f, %f\n", mmax, xmax, ymax);
       if ( mmax<=-MAX_FLOAT/2. ) break;
       // new border vertex found
       xp = xmax;
       yp = ymax;
       sum_x += xmax;
       sum_y += ymax;
       si = shadow->stats.verts;
       shadow->verts[si][0] = xp;
       shadow->verts[si][ycoord] = yp;
       shadow->stats.verts = si+1;
   }
   
   if ( shadow->stats.verts<3 ) {
       fprintf(stderr, "Invalid shadow\n");
       return shadow;
   }
   
   // Add central vertex
   si = shadow->stats.verts;
   float xmed = sum_x / si;
   float ymed = sum_y / si;
   shadow->verts[si][0] = xmed;
   shadow->verts[si][ycoord] = ymed;
   shadow->stats.verts = si+1;
   printf("Shadow has %d vertexes\n", shadow->stats.verts);
   
   // smachdown zz coords to zmin
   for( i=0 ; i<shadow->stats.verts ; i++ )
       shadow->verts[i][zcoord] = zmin;
   
   // alloc counters
   shadow->counts.verts = Malloc(shadow->stats.verts, sizeof(int));
   shadow->counts.norms = Malloc(shadow->stats.norms, sizeof(int));
   shadow->counts.texts = Malloc(shadow->stats.texts, sizeof(int));
   shadow->order.verts  = Malloc(shadow->stats.verts, sizeof(int));
   shadow->order.norms  = Malloc(shadow->stats.norms, sizeof(int));
   shadow->order.texts  = Malloc(shadow->stats.texts, sizeof(int));
   
   // Add faces
   shadow->stats.faces = shadow->stats.verts - 1;
   shadow->faces = Malloc(shadow->stats.faces, sizeof(Face));   
   for( i=0 ; i<shadow->stats.faces ; i++ ) {
       shadow->faces[i].nodes = 3;   // shadow is build using triangles.
       shadow->faces[i].Node = Malloc(3, sizeof(FaceNode));
       shadow->faces[i].Node[0][0] = i+1;
       shadow->faces[i].Node[0][1] = 2;         // Light texture point
       shadow->faces[i].Node[0][2] = 1;     // Up normal vector
       shadow->faces[i].Node[1][0] = i+2;
       shadow->faces[i].Node[1][1] = 2;         // Light texture point
       shadow->faces[i].Node[1][2] = 1;     // Up normal vector
       // 3rd face node is always central vertex
       shadow->faces[i].Node[2][0] = shadow->stats.verts;
       shadow->faces[i].Node[2][1] = 1;         // Dark texture point
       shadow->faces[i].Node[2][2] = 1;     // Up normal vector
   }
   shadow->faces[i-1].Node[1][0] = 1;       // Last polygon uses first vertex 
   if ( Verbose>2 ) printf("Shadow done\n");
   return shadow;
}

ObjFile *CreateShadowObj_Projection(ObjFile *obj) {
    // future version using vertical projection.
    // Under development...
    int ycoord = 2;
    int zcoord = 1;
    if (Verbose>2) fprintf(stderr, "Creating shadow\n");
    
    int i, si;  // face indexs
    ObjFile *shadow = Malloc(1, sizeof(ObjFile));

    shadow->stats.verts = 0;
    shadow->verts = Malloc(obj->stats.verts+1, sizeof(Vert));

    shadow->stats.faces = 0;

    // Add one normal
    shadow->stats.norms = 1;    // Only one norm, facing up.
    shadow->norms = Malloc(1, sizeof(Norm));
    shadow->norms[0][0] = 0.;
    shadow->norms[0][ycoord] = 0.;
    shadow->norms[0][zcoord] = 1.;  // Up vector
    
    // Add two texture coords
    shadow->stats.texts = 2;
    shadow->texts = Malloc(shadow->stats.texts, sizeof(Text));
    shadow->texts[0][0] = 0.01;
    shadow->texts[0][1] = 0.01;
    shadow->texts[1][0] = 0.91;
    shadow->texts[1][1] = 0.91;
    
    // Add one MtlLib
    shadow->stats.libs = 1;
    shadow->libs = Malloc(shadow->stats.libs, sizeof(ObjGroup));
    shadow->libs[0].line = 0;       // before all faces
    shadow->libs[0].name = "shadow.mtl";
    
    // Add one Material
    shadow->stats.mats = 1;
    shadow->mats = Malloc(shadow->stats.mats, sizeof(ObjGroup));
    shadow->mats[0].line = 0;       // before all faces
    shadow->mats[0].name = "shadow";

    shadow->grps = NULL;
    shadow->objs = NULL;
    shadow->shds = NULL;

    shadow->counts.verts = NULL;
    shadow->counts.norms = NULL;
    shadow->counts.texts = NULL;

    shadow->order.verts = NULL;
    shadow->order.texts = NULL;
    shadow->order.norms = NULL;
        
    // Find z min
    float zmin = MAX_FLOAT;
    for( i=0 ; i<obj->stats.verts ; i++ )
    if ( obj->verts[i][zcoord]<zmin )
        zmin = obj->verts[i][zcoord];

    /* Flags may be use to reject vertexes
    unsigned char flags[obj->stats.verts];
    // Reset all flags;
    for( i=0 ; i<obj->stats.verts ; i++ )
        flags[i] = 0;
    */

    // Find leftest vertex
    float xmin = MAX_FLOAT;
    float ymin = MAX_FLOAT;
    int   vmin = -1;
    int   fmin = -1;
    //int   nmin = -1;

    int fi, ni;
    for( fi=0 ; fi<obj->stats.faces ; fi++ )
    for( ni=0 ; ni<obj->faces[fi].nodes ; ni++ ) {
    int i = obj->faces[fi].Node[ni][0]-1;
    if ( obj->verts[i][0]<xmin ) {
        xmin = obj->verts[i][0];
        ymin = obj->verts[i][ycoord];
        vmin = i;
        fmin = fi;
        //nmin = ni;
    }
    }
   
    // Find surrounding convex polygon
    // cursor variables. Register current point.
    float xp = xmin;
    float yp = ymin;
    double ap = PI/2.;
    int vp = vmin;
    int fp = fmin;
    //int np = nmin;

    float sum_x = xp;
    float sum_y = yp;
   
    {    // Register first vertex
        int si = shadow->stats.verts;
        shadow->verts[si][0] = xp;
        shadow->verts[si][ycoord] = yp;
        shadow->stats.verts = si+1;
    }
    //if ( vp>=0 ) flags[vp] = 1;
    printf("starting vertex %d, %f, %f\n", vp, xp, yp);
   
    // Looking for boundary
    while (1) {
        double amax = -MAX_FLOAT;
        double ymax = -MAX_FLOAT;
        double xmax = -MAX_FLOAT;
        int    vmax = -1;
        int    fmax = -1;
        int    nmax = -1;
    
        // Find next edge, with greatest slope
        // Edge must contain current vertex
        for( int fi=0 ; fi<obj->stats.faces ; fi++ ) {
            //printf("   testing face %d, nodes %d\n", fi, obj->faces[fi].nodes);
            int ni, vidx=-1;
            for( ni=0 ; ni<obj->faces[fi].nodes ; ni++ ) {
                vidx = obj->faces[fi].Node[ni][0]-1;
                assert(vidx>=0);
                if ( vp==vidx ) break;   // vertex_idx matches
                if ( xp==obj->verts[vidx][0] && yp==obj->verts[vidx][ycoord] )
                    break;      // vertex position matches
            }
            if ( ni>=obj->faces[fi].nodes ) continue; // No matches in this face. Next face.
            printf(" # face match %d, node %d vidx %d, ap=%f\n", fi, ni, vidx,ap);
            // Match. Face[fi] contains current vertex in node ni.
            // Since faces are convex, analyze only edges to adjacent vertexes.
            int ei1, ei2;
            // previous node
            ei1 = ni-1;
            if ( ei1<0 ) ei1 = obj->faces[fi].nodes-1;
            // next node
            ei2 = ni+1;
            if ( ei2>=obj->faces[fi].nodes ) ei2 = 0;
            printf("ei1 ei2 %d %d\n", ei1, ei2); 
            // previous node
            vidx = obj->faces[fi].Node[ei1][0]-1;
            printf("    previous node, vidx=%d\n", vidx);
            assert(vidx<obj->stats.verts);
            assert(vidx>=0);
            if ( vp!=vidx ) {   // Avoid duplicated vertexes
                float x = obj->verts[vidx][0];
                float y = obj->verts[vidx][ycoord];
                if ( x!=xp || y!=yp ) {
                    double ang = atan2(y-yp,x-xp);
                    double dang = AngNormalize(ang-ap);
                    printf(" analysing edge to %f %f, ang %f (%.2f), dang %f (%.2f)\n", x, y, ang, ang*180./PI, dang, dang*180./PI);
                    if ( dang<PI && dang>amax ) {
                        amax = dang;
                        xmax = x;
                        ymax = y;
                        vmax = vidx;
                        fmax = fi;
                        nmax = ni;
                    }       
                }
            }
            // next node
            vidx = obj->faces[fi].Node[ei2][0]-1;
            assert(vidx>=0);
            if ( vp!=vidx ) {   // Avoid duplicated vertexes
                float x = obj->verts[vidx][0];
                float y = obj->verts[vidx][ycoord];
                //printf(" analysing edge to %f %f\n",x,y);
                if ( x!=xp || y!=yp ) {
                    double ang = atan2(y-yp,x-xp);
                    double dang = AngNormalize(ang-ap);
                    printf(" analysing edge to %f %f, ang %f (%.2f), dang %f (%.2f)\n", x, y, ang, ang*180./PI, dang, dang*180./PI);
                    if ( dang<PI && dang>amax ) {
                        amax = dang;
                        xmax = x;
                        ymax = y;
                        vmax = vidx;
                        fmax = fi;
                        nmax = ni;
                    }       
                }
            }
        } // for all faces

        if ( amax<-MAX_FLOAT/2. ) break;  // no edges to follow
        if ( vmax==vmin ) break;    // round complete
   
        // best edge got
        double dx = xmax-xp;
        double dy = ymax-yp;
        ap = amax = atan2(dy,dx);
        // assert(vp>=0);
                
        printf(" loop exited. amax=%f, vmax=%d, ap= %f (%.2f) len %f\n",amax,vmax,ap,ap*180/PI,sqrt(dx*dx+dy*dy));
        
              //////
        // if (vmax==34682) break;    //DEBUG
        //////////////////////////////////////////////////////////////
        
        double tdist = -1;
        do {    
            // Follow this edge.
            // Next vertex will be vmax(xmax,ymax) only if edge is not intersected
            // Must test intersection with all edges in every face.
            int nfi, nni;
            printf("  testing intersection with xp,yp - xmax,ymax (%f %f) (%f %f)\n",xp,yp,xmax,ymax);
            tdist = FindClosestIntersection(xp,yp,xmax,ymax,obj,&nfi,&nni);
            //printf("    intersection dist: %f\n", tdist);
            if ( tdist>0. && tdist<1.) {   // Intersection found
                assert(nfi>=0);
                assert(nni>=0);
                printf("    intersection dist: %f\n", tdist);
                //double len = sqrt(dx*dx+dy*dy);
                xp += dx*tdist;
                yp += dy*tdist;
                vp = -1;    // Unexistant/new vertex
                fp = nfi;
                //np = nni;
                
                sum_x += xp;
                sum_y = yp;
    
                {    // Register vertex            
                    printf("## got shadow new vertex (%f %f)\n", xp, yp);
                    int si = shadow->stats.verts;
                    shadow->verts[si][0] = xp;
                    shadow->verts[si][ycoord] = yp;
                    shadow->stats.verts = si+1;
                    if ( shadow->stats.verts>obj->stats.verts ) {
                        fprintf(stderr, "Too many vertexes in shadow\n");
                        return shadow;
                    }
                }
                printf("nfi %d nni %d\n", nfi, nni);
                fmax = fp;
                nmax = nni+1;
                if ( nmax>=obj->faces[fmax].nodes ) nmax = 0;
                vmax = obj->faces[fmax].Node[nmax][0]-1;
                xmax = obj->verts[vmax][0];
                ymax = obj->verts[vmax][ycoord];
                dx = xmax-xp;
                dy = ymax-yp;
                amax = atan2(dy,dx);
                double dang = AngNormalize(amax-ap);
                printf("fmax %d nmax %d vmax %d xmax %f ymax %f amax %f ap %f dang %f\n", fmax, nmax, vmax, xmax, ymax, amax,ap,dang);
                if ( dang<0. ) {
                    dang += PI;
                    ap = AngNormalize(ap-PI);
                    nmax = nni;
                    vmax = obj->faces[fmax].Node[nmax][0]-1;
                    xmax = obj->verts[vmax][0];
                    ymax = obj->verts[vmax][ycoord];
                    dx = xmax-xp;
                    dy = ymax-yp;
                    amax = atan2(dy, dx);
                    printf(" fmax %d nmax %d vmax %d xmax %f ymax %f amax %f ap %f dang %f\n", fmax, nmax, vmax, xmax, ymax, amax,ap,dang);
                }
    
            }
        } while ( tdist>0. && tdist<1.);
        
        {  // no nore intersections found
            xp = xmax;
            yp = ymax;
            vp = vmax;
            fp = fmax;
            //np = nmax;
            ap = amax;
        }
        // Register new vertex
        printf("## got shadow MAIN vertex %d (%f %f), ap=%f\n",vp,xp,yp,ap);
        sum_x += xp;
        sum_y += yp;
        int si = shadow->stats.verts;
        shadow->verts[si][0] = xp;
        shadow->verts[si][ycoord] = yp;
        shadow->stats.verts = si+1;
        if ( shadow->stats.verts>obj->stats.verts ) {
            fprintf(stderr,"Too many vertexes in shadow\n");
            return shadow;
        }
    }  // while(1) 
  
    printf("Shadow add central vertex. normal %f %f %f\n", shadow->norms[0][0], shadow->norms[0][1], shadow->norms[0][2]);
 
    if ( shadow->stats.verts<3 ) {
       fprintf(stderr,"Invalid shadow\n");
       return shadow;
    }
   
    // Add central vertex
    float xmed = sum_x / shadow->stats.verts;
    float ymed = sum_y / shadow->stats.verts;
    si = shadow->stats.verts;
    shadow->verts[si][0] = xmed;
    shadow->verts[si][ycoord] = ymed;
    shadow->stats.verts = si+1;
    printf("Shadow has %d vertexes\n", shadow->stats.verts);
 
    // smachdown zz coords to zmin
    for( i=0 ; i<shadow->stats.verts ; i++ )
        shadow->verts[i][zcoord] = zmin;
 
    // alloc counters
    shadow->counts.verts = Malloc(shadow->stats.verts, sizeof(int));
    shadow->counts.norms = Malloc(shadow->stats.norms, sizeof(int));
    shadow->order.verts = Malloc(shadow->stats.verts, sizeof(int));
    shadow->order.norms = Malloc(shadow->stats.norms, sizeof(int));
    shadow->counts.texts = Malloc(shadow->stats.texts, sizeof(int));
    shadow->order.texts = Malloc(shadow->stats.texts, sizeof(int));
 
    // Add faces
    shadow->stats.faces = shadow->stats.verts - 1;
    shadow->faces = Malloc(shadow->stats.faces, sizeof(Face));   
    for( i=0 ; i<shadow->stats.faces ; i++ ) {
        shadow->faces[i].nodes = 3;   // shadow is build using triangles.
        shadow->faces[i].Node = Malloc(3,sizeof(FaceNode));
        shadow->faces[i].Node[0][0] = i+1;
        shadow->faces[i].Node[0][1] = 2;        // Light texture point
        shadow->faces[i].Node[0][2] = 1;
        shadow->faces[i].Node[1][0] = i+2;
        shadow->faces[i].Node[1][1] = 2;        // Light texture point
        shadow->faces[i].Node[1][2] = 1;
        shadow->faces[i].Node[2][0] = shadow->stats.verts; // last vertex
        shadow->faces[i].Node[2][1] = 1;        // Dark texture point
        shadow->faces[i].Node[2][2] = 1;
    }
    shadow->faces[i-1].Node[1][0] = 1;   // Last polygon uses first vertex 
    if ( Verbose>2 ) printf("Shadow done\n");

    return shadow;
}

void SaveObjFile(char *fname, ObjFile *obj) {
    FILE *fout;
    fout = stdout;
    int i;
    int nv=0, nt=0, nn=0;
    
    if ( fname && *fname ) {
        if ( Verbose ) printf("Writing file '%s'\n",fname);
        fout = fopen(fname,"w");
        if ( ! fout ) {
            fprintf(stderr,"Cannot open output file '%s'\n", fname);
            exit(5);
        }
    }
    fprintf(fout,"# file written by ObjTool\n\n");
        
    // Save Vertexes UseMtl: 0
    assert(obj->counts.verts);
    for( i=0 ; i<obj->stats.verts ; i++ ) {
        if ( Verbose>30 ) printf("vertex %d-%d\n", i, nv);
        if ( obj->counts.verts[i] > 0 ) {
            fprintf(fout,"v %f %f %f\n", obj->verts[i][0], obj->verts[i][1], obj->verts[i][2]);
            nv++;
        }
    }
    if ( Verbose>10 ) printf("    %d vertexes saved\n", nv);
    
    // Save Texture Coords
    if ( genTextCoords ) {
        for( i=0 ; i<obj->stats.verts ; i++ ) {
            if ( obj->counts.verts[i] > 0 ) {
                fprintf(fout,"vt %f %f\n", obj->verts[i][0], obj->verts[i][1]);
                nt++;
            }
        }
    }
    else {
        for( i=0 ; i<obj->stats.texts ; i++ ) {
            if ( obj->counts.texts[i] > 0 ) {
                fprintf(fout,"vt %f %f\n", obj->texts[i][0], obj->texts[i][1]);
                nt++;
            }
        }
    }
    if ( Verbose>10 ) printf("    %d texture vertexes saved\n", nt);
    
    // Save Norms
    for( i=0 ; i<obj->stats.norms ; i++ ) 
        if ( obj->counts.norms[i] > 0 ) {
            fprintf(fout,"vn %f %f %f\n", obj->norms[i][0], obj->norms[i][1], obj->norms[i][2]);
            nn++;
        }
    if ( Verbose>10 ) printf("    %d normals saved\n", nn);
    
    // Save faces
    int gidx=0,oidx=0,midx=0,lidx=0;
    int saved_faces = 0;
    for( i=0 ; i<obj->stats.faces ; i++ ) {
        int n;
        if ( Verbose>12) printf("  Working on face %d\n",i);
        if ( lidx<obj->stats.libs && i>=obj->libs[lidx].line ) {
            while( i<obj->stats.faces && obj->faces[i].nodes <= 0 ) i++;
            while( lidx<obj->stats.libs && i>=obj->libs[lidx].line ) lidx++;
            fprintf(fout,"mtllib %s\n",obj->libs[lidx-1].name);
        }
        if ( oidx<obj->stats.objs && i>=obj->objs[oidx].line ) {
            while( i<obj->stats.faces && obj->faces[i].nodes <= 0 ) i++;
            while( oidx<obj->stats.objs && i>=obj->objs[oidx].line ) oidx++;
            fprintf(fout,"o %s\n",obj->objs[oidx-1].name);
        }
        if ( gidx<obj->stats.grps && i>=obj->grps[gidx].line ) {
            while( i<obj->stats.faces && obj->faces[i].nodes <= 0 ) i++;
            while( gidx<obj->stats.grps && i>=obj->grps[gidx].line ) gidx++;
            fprintf(fout,"g %s\n",obj->grps[gidx-1].name);
        }
        if ( midx<obj->stats.mats && i>=obj->mats[midx].line ) {
            while( i<obj->stats.faces && obj->faces[i].nodes <= 0 ) i++;
            while( midx<obj->stats.mats && i>=obj->mats[midx].line ) midx++;
            fprintf(fout,"usemtl %s\n",obj->mats[midx-1].name);
        }
        if ( i>=obj->stats.faces ) break;
        if ( obj->faces[i].nodes <= 0 ) continue;  // face erased
        fprintf(fout,"f ");
        saved_faces += 1;
        for( n=0 ; n<obj->faces[i].nodes ; n++ ) {
            fprintf(fout,"%d", VertIndex(obj,obj->faces[i].Node[n][0],nv) );
            fprintf(fout,"/");
    
            if ( genTextCoords ) {
                fprintf(fout,"%d", VertIndex(obj,obj->faces[i].Node[n][0],nv) );
            }
            else {  
                if ( obj->faces[i].Node[n][1]!=NULL_IDX )
                    fprintf(fout,"%d", TextIndex(obj,obj->faces[i].Node[n][1],nt) );
            }
            if ( obj->faces[i].Node[n][2]!=NULL_IDX ) {
                fprintf(fout,"/");
                fprintf(fout,"%d", NormIndex(obj,obj->faces[i].Node[n][2],nn) );
            }
            fprintf(fout," ");
        }
        fprintf(fout,"\n");
    }
    if ( Verbose>9 ) printf("    %d faces saved\n", saved_faces);
    if ( fout!=stdout ) fclose(fout); 
}


void ExplodeOutputFile(char *OutputFile, ObjFile *obj) {

    // Split ObjFile *obj into parts.
    // Each part is saved to file named OutputFile_n_imat_iobj_igrp.obj

    if ( ! OutputFile ) {
        fprintf(stderr,"Explode requires the output file name using the -O option.\n");
        return;
    }
    assert(obj); 

    int i;
    int iface = 0;
    int imat = 0;
    int iobj = 0;
    int igrp = 0;
    #define MAX_FILE_NAME_LEN (480)
    char fname[MAX_FILE_NAME_LEN];
    unsigned int fnum = 0;
    
    do {
        // printf("    Adding face %d\n", iface);
        // Create new obj_part
        ObjFile *obj_part;
        obj_part = Malloc(1,sizeof(ObjFile));

        JoinObjFiles(1,obj,obj_part);
    
        SetUseCounters(obj_part);
        
        // Clear all previous faces
        for( i=0 ; i<iface ; i++ )
        obj_part->faces[i].nodes = 0;
            
        // Find the part end
        // Part end when at least one change detected from {imat,iobj,igrp}
        int part_cut = 0;
        while ( ! part_cut && iface<obj_part->stats.faces ) {
            //printf("  iface:%d\n", iface);
            //if ( iface>0 ) {        // ignore directive before the first face
                while ( imat<obj_part->stats.mats && iface==obj_part->mats[imat].line ) {
                    imat++;
                    if ( iface>0 ) part_cut = 1;
                }
                while ( iobj<obj_part->stats.objs && iface==obj_part->objs[iobj].line ) {
                    iobj++;
                    if ( iface>0 ) part_cut = 1;
                }
                while ( igrp<obj_part->stats.grps && iface==obj_part->grps[igrp].line ) {
                    igrp++;
                    if ( iface>0 )part_cut = 1;
                }
            // }
            if ( ! part_cut) iface++;
        }
        
        // End of part detected.
        // Clear all remaining faces.
        for( i=iface ; i<obj_part->stats.faces ; i++ )
            obj_part->faces[i].nodes = 0;

        SetUseCounters(obj_part);
        SetIndexs( obj_part );

        snprintf(fname,MAX_FILE_NAME_LEN,"%s_%u_%d_%d_%d.obj", OutputFile, fnum, imat, iobj, igrp);
        if (Verbose>2) printf("  Saving PartFile '%s'\n",fname); 
        //printf("    iface:%d\n", iface);
        SaveObjFile(fname, obj_part );

        Free(obj_part);

        fnum++;  // next file partition
            //iface++;
    } while ( iface < obj->stats.faces );   

}

void FreeObjFile(ObjFile *obj) {
    int i;
    if ( Verbose>4 ) printf("Free\n");
    
    if ( Verbose>5 ) printf("  free stats\n");
    for( i=0 ; i<obj->stats.faces ; i++ )
    Free(obj->faces[i].Node);
    for( i=0 ; i<obj->stats.grps ; i++ )
    Free(obj->grps[i].name);
    for( i=0 ; i<obj->stats.libs ; i++ )
    Free(obj->libs[i].name);
    for( i=0 ; i<obj->stats.mats ; i++ )
    Free(obj->mats[i].name);
    for( i=0 ; i<obj->stats.shds ; i++ )
    Free(obj->shds[i].name);
    for( i=0 ; i<obj->stats.objs ; i++ )
    Free(obj->objs[i].name);

    if ( Verbose>5 ) printf("  free faces\n");
    Free(obj->faces); 
    
    if ( Verbose>5 ) printf("  free texts\n");
    Free(obj->texts);
    
    if ( Verbose>5 ) printf("  free norms\n");
    Free(obj->norms);
    
    if ( Verbose>5 ) printf("  free grps\n");
    Free(obj->grps);
    
    if ( Verbose>5 ) printf("  free mats\n");
    Free(obj->mats);
    
    if ( Verbose>5 ) printf("  free objs\n");
    Free(obj->objs);
    
    if ( Verbose>5 ) printf("  free shds\n");
    Free(obj->shds);
    
    if ( Verbose>5 ) printf("  free lib\n");
    Free(obj->libs);
    
    if ( Verbose>5 ) printf("  free counts\n");   
    Free(obj->counts.verts);
    Free(obj->counts.norms);
    Free(obj->counts.texts);
    
    if ( Verbose>5 ) printf("free order\n");  
    Free(obj->order.verts);
    Free(obj->order.texts);
    Free(obj->order.norms);
    
    if ( Verbose>4 ) printf("free end\n");
}



void SetUseCounters(ObjFile *obj) {
    // Counts uses of each vertex, norm and texture coord.
    // results are stored in obj->counts structure
    assert(obj);
    assert(obj->counts.verts);
    int i;
    for( i=0 ; i<obj->stats.verts ; i++ )
        obj->counts.verts[i] = 0;
    for( i=0 ; i<obj->stats.norms ; i++ )
        obj->counts.norms[i] = 0;
    for( i=0 ; i<obj->stats.texts ; i++ )
        obj->counts.texts[i] = 0;
/*
    for( i=0 ; i<obj->stats.mats ; i++ )
        obj->counts.mats[i] = 0;
    for( i=0 ; i<obj->stats.objs ; i++ )
        obj->counts.objs[i] = 0;
    for( i=0 ; i<obj->stats.shads ; i++ )
        obj->counts.shds[i] = 0;
    for( i=0 ; i<obj->stats.grps ; i++ )
        obj->counts.grps[i] = 0;
*/
    int n;
    for( i=0 ; i<obj->stats.faces ; i++ ) {
        for( n=0 ; n<obj->faces[i].nodes ; n++ ) {
            int v = obj->faces[i].Node[n][0];
            if ( v<1 || v>obj->stats.verts ) {
                fprintf(stderr, "Vertice %d not found on face %d, node %d/%d.\n", v, i, n, obj->faces[i].nodes);
                continue;
            }
            obj->counts.verts[v-1] += 1;
            v = obj->faces[i].Node[n][1];
            if ( v!=NULL_IDX ) {
                if ( v<1 || v>obj->stats.texts ) {
                    fprintf(stderr, "Texture Coord %d not found on face %d, node %d.\n", v, i, n);
                    continue;
                }
                obj->counts.texts[v-1] += 1;
            }
            v = obj->faces[i].Node[n][2];
            if ( v!=NULL_IDX ) {
                if ( v<1 || v>obj->stats.norms ) {
                    fprintf(stderr, "Normal %d not found on face %d, node %d.\n", v, i, n);
                    continue;
                }
                obj->counts.norms[v-1] += 1;
            }
        }
    } 
}

static int CleanFaces(ObjFile *obj) {
    int i, n;
    int midx = 0, oidx = 0, gidx = 0;   // current indexes.
                        // 0(zero) means material/object/group not yet defined
                        // Faces without defined material/object/group always fail the test.
    int nRemFaces = 0;
    for( i=0 ; i<obj->stats.faces ; i++ ) {
        // Remove face if material matches/does not match selected Material
        if ( Material!=NULL ) {
            while ( midx<=obj->stats.mats && i==obj->mats[midx].line )     midx++;
            if ( Negate ) {
                if ( midx && strcmp(Material,obj->mats[midx-1].name)==0 ) {
                    obj->faces[i].nodes = 0;
                    nRemFaces++;
                    if ( Verbose>3 )
                        printf("  Face %d removed.\n", i);
                    continue;
                }
            }
            else {
                if ( ! midx || strcmp(Material,obj->mats[midx-1].name)!=0 ) {
                    obj->faces[i].nodes = 0;
                    nRemFaces++;
                    if ( Verbose>3 )
                        printf("  Face %d removed.\n", i);
                    continue;
                }
            }
        }
        if ( SelObjects>0 ) {
            // Remove face if object matches/does not match selected Object
            while ( oidx<=obj->stats.objs && i==obj->objs[oidx].line )     oidx++;
            // Current Test implementation is limited to 1 object ( SelectObject[0] )
            if ( Negate ) {
                if ( oidx && strcmp(SelectObject[0],obj->objs[oidx-1].name)==0 ) {
                    obj->faces[i].nodes = 0;
                    nRemFaces++;
                    if ( Verbose>3 )
                        printf("  Face %d removed.\n", i);
                    continue;
                }
            }
            else {
                if ( ! oidx || strcmp(SelectObject[0],obj->objs[oidx-1].name)!=0 ) {
                    obj->faces[i].nodes = 0;
                    nRemFaces++;
                    if ( Verbose>3 )
                        printf("  Face %d removed.\n", i);
                    continue;
                }
            }
        }
        if ( SelGroups>0 ) {
            // Remove face if group matches/does not match selected Group
            while ( gidx<=obj->stats.grps && i==obj->grps[gidx].line )     gidx++;
            // Current Test implementation is limited to 1 group ( SelectGroup[0] )
            if ( Negate ) {
                if ( gidx && strcmp(SelectGroup[0],obj->grps[gidx-1].name)==0 ) {
                    obj->faces[i].nodes = 0;
                    nRemFaces++;
                    if ( Verbose>3 )
                        printf("  Face %d removed.\n", i);
                    continue;
                }
            }
            else {
                if ( ! gidx || strcmp(SelectGroup[0],obj->grps[gidx-1].name)!=0 ) {
                    obj->faces[i].nodes = 0;
                    nRemFaces++;
                    if ( Verbose>3 )
                        printf("  Face %d removed.\n", i);
                    continue;
                }
            }
        }
        // Remove face if it includes inexistent vertexs.
        // Maybe try other options in the future:
        //     - Keep faces that have at least three selected vertexes, removing unselected vertexes.
        //     - Keep faces with all vertexes if they have at least one selected vertex. 
        //     - Adjust position of unselected vertexes so that they go to the border of the selected area.
        //       This is similar to the implemented SolidCut option.
        //       But SolidCut works only on vertexes and does not know any thing about faces.
        for( n=0 ; n<obj->faces[i].nodes ; n++ ) {
            int v = obj->faces[i].Node[n][0];
            if ( v<1 || v>obj->stats.verts  || obj->counts.verts[v-1]<0 ) {
                obj->faces[i].nodes = 0;
                nRemFaces++;
                if ( Verbose>=7 )
                    fprintf(stderr, "  Face %d removed.v=%d,count=%d\n", i, v, obj->counts.verts[v-1] );
                else if ( Verbose>3 )
                    printf("  Face %d removed.\n", i);
                break;
            }
        }
    } 
    if ( Verbose ) printf("Faces removed: %d/%d (%.2f%%)\n", nRemFaces, obj->stats.faces, 100.f * (float)nRemFaces / (float)obj->stats.faces );
    return nRemFaces;
}


void SetIndexs(ObjFile *obj) {
    // Assign an order number to each used vertex/texture coord/norm
    // Unused vertex/texture coord/norm will get no order number.
    // Order number are stored into obj->order.* arrays.
    int i;
    int nv=0, nt=0, nn=0;
    for( i=0 ; i<obj->stats.verts ; i++ ) 
        if ( obj->counts.verts[i]>0 )   obj->order.verts[i] = nv++;

    for( i=0 ; i<obj->stats.texts ; i++ ) 
        if ( obj->counts.texts[i]>0 )   obj->order.texts[i] = nt++;

    for( i=0 ; i<obj->stats.norms ; i++ ) 
        if ( obj->counts.norms[i]>0 )   obj->order.norms[i] = nn++;
}

void FilterVerts(ObjFile *obj, Vert min, Vert max) {
    // Keep only Verts inside AABB(min,max).
    // Unselected Verts are marked as unused (counts=NULL:IDX) or
    // translated to AABB border, depending on global SolidCut.
    int i, c;
    int nVertsCut = 0;
    for( i=0 ; i<obj->stats.verts ; i++ ) {
        int vertCut = 0;
        for( c=0 ; c<3 ; c++ ) {
            if ( obj->verts[i][c]<min[c] ) {
                vertCut = 1;
                if ( ! SolidCut ) {
                    obj->counts.verts[i] = NULL_IDX;
                    break;
                }
                else obj->verts[i][c] = min[c];
            }
            if ( obj->verts[i][c]>max[c] ) {
                vertCut = 1;
                if ( ! SolidCut ) {
                    obj->counts.verts[i] = NULL_IDX;
                    break;
                }
                else obj->verts[i][c] = max[c];
            }
        }
        nVertsCut += vertCut;
    }
    if ( Verbose>0 ) printf("Vertexes cropped out: %d/%d (%.2f%%)\n", nVertsCut, obj->stats.verts, 100.f * (float)nVertsCut / (float)obj->stats.verts );
}


void  vertCrossProduct(ObjFile *obj, int nFace, int n0, Vert cross) {
    // Cross product of edges in vertex (n0+1).
    // Cross=nullVec when edges are colinear.   
    
    Face face = obj->faces[nFace];
    int n1 = (n0+1) % face.nodes;
    int n2 = (n0+2) % face.nodes;
    int i0 = face.Node[n0][0];
    int i1 = face.Node[n1][0];
    int i2 = face.Node[n2][0];

    if ( i1<1 || i1>obj->stats.verts ) 
        return;
    if ( i2<1 || i2>obj->stats.verts ) 
        return; 

    vec3f v0, v1, v2;
    vec3f_copy(v0, obj->verts[i0-1]);
    vec3f_copy(v1, obj->verts[i1-1]);
    vec3f_copy(v2, obj->verts[i2-1]);
    
    //printf("        v0 %f,%f,%f\n", v0[0],v0[1],v0[2]);
    //printf("        v1 %f,%f,%f\n", v1[0],v1[1],v1[2]);
    //printf("        v2 %f,%f,%f\n", v2[0],v2[1],v2[2]);
    
    Vert A, B;
    vec3_subed(A, v1, v0);
    vec3_subed(B, v2, v1); 
    
    //printf("      A %f,%f,%f\n", A[0],A[1],A[2]);
    //printf("      B %f,%f,%f\n", B[0],B[1],B[2]);
    vec3_cross(cross, A, B);   
}

void ProcVerts(ObjFile *obj) {
    // Apply transformations to vertexs
    // Transformations orders is: Scale / Translation / Rotation
    int i, c;
    for( i=0 ; i<obj->stats.verts ; i++ ) {
        if ( obj->counts.verts[i] > 0 ) {
            for( c=0 ; c<3 ; c++ ) {
                obj->verts[i][c] *= Scale[c];
                obj->verts[i][c] += Translate[c];
            }
            float x,y,z;
            x = obj->verts[i][0];
            y = obj->verts[i][1];
            obj->verts[i][0] =  x*CosF(Rotate[2]) - y*SinF(Rotate[2]);
            obj->verts[i][1] =  x*SinF(Rotate[2]) + y*CosF(Rotate[2]);
            x = obj->verts[i][0];
            z = obj->verts[i][2];
            obj->verts[i][0] =  x*CosF(Rotate[1]) + z*SinF(Rotate[1]);
            obj->verts[i][2] = -x*SinF(Rotate[1]) + z*CosF(Rotate[1]);
            y = obj->verts[i][1];
            z = obj->verts[i][2];
            obj->verts[i][1] =  y*CosF(Rotate[0]) - z*SinF(Rotate[0]);
            obj->verts[i][2] =  y*SinF(Rotate[0]) + z*CosF(Rotate[0]);
        }
    }

    // Apply Rotations to Normals
    for( i=0 ; i<obj->stats.norms ; i++ ) {
        float x,y,z;
        x = obj->norms[i][0];
        y = obj->norms[i][1];
        obj->norms[i][0] =  x*CosF(Rotate[2]) - y*SinF(Rotate[2]);
        obj->norms[i][1] =  x*SinF(Rotate[2]) + y*CosF(Rotate[2]);
        x = obj->norms[i][0];
        z = obj->norms[i][2];
        obj->norms[i][0] =  x*CosF(Rotate[1]) + z*SinF(Rotate[1]);
        obj->norms[i][2] = -x*SinF(Rotate[1]) + z*CosF(Rotate[1]);
        y = obj->norms[i][1];
        z = obj->norms[i][2];
        obj->norms[i][1] =  y*CosF(Rotate[0]) - z*SinF(Rotate[0]);
        obj->norms[i][2] =  y*SinF(Rotate[0]) + z*CosF(Rotate[0]);
            /*
            // Moved to new section. 
        if ( InvertNormals & 0b0001 )
            for( c=0 ; c<3 ; c++)
                obj->norms[i][c] *= -1.;
        */
    }
    
    int nInvNormals = 0;
    if ( InvertNormals & 0b001 ) {  // -Ni
        for( i=0 ; i<obj->stats.norms ; i++ ) {
            nInvNormals++;
            for( c=0 ; c<3 ; c++)
            obj->norms[i][c] *= -1.;
        }
        if ( Verbose ) printf("Normals inverted: %d\n", nInvNormals);
    }
    
    int nInvFaces = 0;
    if ( InvertNormals & 0b110 ) {  // -Nc  or  -Ns
        for( i=0 ; i<obj->stats.faces ; i++ ) {
            vec3f polyCross;
            vec3f_set(polyCross, 0., 0., 0.);
            vec3f polyNormal;
            vec3f_set(polyNormal, 0., 0., 0.);
            for( int n=0 ; n<obj->faces[i].nodes ; n++ ) {
                Vert vertCross;
                vertCrossProduct(obj, i, n, vertCross);
                //printf("    vert %d Cross %f,%f,%f\n", n, vertCross[0], vertCross[1], vertCross[2]);
                vec3_add(polyCross, vertCross);
                int iNormal = obj->faces[i].Node[n][2];
                //printf("      iNormal %d\n", iNormal);
                vec3f vertNormal;
                vec3f_copy(vertNormal, obj->norms[iNormal-1]);
                if ( InvertNormals & 0b0010 ) {
                    // Compare and correct Normal at vertex.
                    float vDir = vec3_dot(vertNormal, vertCross);
                    if ( vDir < 0. ) {
                        // printf("  Invertendo normal %d\n", iNormal);
                        vec3_invert(obj->norms[iNormal-1]);
                        vec3f_copy(vertNormal, obj->norms[iNormal-1]);
                        nInvNormals++;
                        if ( Verbose>4 ) printf("  Normal %d inverted\n", iNormal);
                    }
                }
                vec3f_add(polyNormal, vertNormal);
            }
            //printf(" polyNormal %f,%f,%f\n", polyNormal[0], polyNormal[1], polyNormal[2]);
            //printf(" polyCross %f,%f,%f\n", polyCross[0], polyCross[1], polyCross[2]);
            float fDir = vec3_dot(polyNormal, polyCross);
            //printf("fdir = %f\n", fDir);
            if ( fDir < 0. ) {
                // poligono invertido.
                if ( InvertNormals & 0b0100 ) {
                    // Inverte ordem dos vertices da face
                    if ( Verbose>4 ) printf("  Face %d inverted\n", i);
                    int n1, n2;
                    for( n1=0, n2=obj->faces[i].nodes-1 ; n1<n2 ; n1++, n2-- ) {
                        int temp;
                        //printf("  invertendo ordem do poligono %d, n1:%d, n2:%d \n", i, n1, n2);
                        for( int r=0 ; r<3 ; r++ ) {
                            temp = obj->faces[i].Node[n1][r];
                            obj->faces[i].Node[n1][r] = obj->faces[i].Node[n2][r];
                            obj->faces[i].Node[n2][r] = temp;
                        }
                    }
                    nInvFaces++;
                }
            }   
        }
        if ( Verbose ) {
            if ( InvertNormals & 0b0100 ) 
                printf("Inverted faces: %d/%d (%.2f%%)\n", nInvFaces, obj->stats.faces, 100.f*(float)nInvFaces/(float)obj->stats.faces );       
            if ( InvertNormals & 0b0010 ) 
                printf("Inverted faces: %d/%d (%.2f%%)\n", nInvNormals, obj->stats.norms, 100.f*(float)nInvNormals/(float)obj->stats.norms );
        }
    }
}


void Usage() {
    fprintf(stderr,"Usage:\n\tObjTool [options] inputfile1.obj inputfile2.obj ...\n");
    fprintf(stderr,"\t\t-i                 file info\n");
    fprintf(stderr,"\t\t-xmin value        crop faces at xmin\n");
    fprintf(stderr,"\t\t-xmax value        crop faces at xmax\n");
    fprintf(stderr,"\t\t-ymin value        crop faces at ymin\n");
    fprintf(stderr,"\t\t-ymax value        crop faces at ymax\n");
    fprintf(stderr,"\t\t-zmin value        crop faces at zmin\n");
    fprintf(stderr,"\t\t-zmax value        crop faces at zmax\n");
    fprintf(stderr,"\t\t-o name            select object\n");
    fprintf(stderr,"\t\t-g name            select obj group\n");
    fprintf(stderr,"\t\t-m name            select material\n");
    fprintf(stderr,"\t\t-tx value          translate x\n");
    fprintf(stderr,"\t\t-ty value          translate y\n");
    fprintf(stderr,"\t\t-tz value      translate z\n");
    fprintf(stderr,"\t\t-sx value      scale x\n");
    fprintf(stderr,"\t\t-sy value      scale y\n");
    fprintf(stderr,"\t\t-sz value      scale z\n");
    fprintf(stderr,"\t\t-s value       global scale\n"); 
    fprintf(stderr,"\t\t-rx ang        rotate ang degrees around xx axis\n");
    fprintf(stderr,"\t\t-ry ang        rotate ang degrees around yy axis\n");
    fprintf(stderr,"\t\t-rz ang        rotate ang degrees around zz axis\n");
    fprintf(stderr,"\t\t-Ni            invert all normals\n");
    fprintf(stderr,"\t\t-Nc            invert wrong facing normals\n");
    fprintf(stderr,"\t\t-Ns            resort face vertexes to match poligon normal\n");
    fprintf(stderr,"\t\t-Nx            x can be a number [0,7] defining a bitwise combination of 'sci' options for processing normals\n");
    fprintf(stderr,"\t\t-R                 use relative coords\n");
    fprintf(stderr,"\t\t-M         do not output mtllib directives\n");
    fprintf(stderr,"\t\t-c                 solid cut\n");
    fprintf(stderr,"\t\t-n         negate face filter condition( -g,-m )\n");
    fprintf(stderr,"\t\t-fz x y        find Z coord for x,y vertical\n");
    fprintf(stderr,"\t\t-fy x z        find Y coord for x,z vertical\n");
    fprintf(stderr,"\t\t-O outfile         output to outfile (default: stdout)\n");
    fprintf(stderr,"\t\t-e         explode outfile into single objects\n");
    fprintf(stderr,"\t\t-Sy shadow_file    cast shadow in Y direction and output to shadow_file (default: no shadow ouput)\n");
    fprintf(stderr,"\t\t-Sz shadow_file    cast shadow in Z direction and output to shadow_file (default: no shadow ouput)\n");
    fprintf(stderr,"\t\t-S shadow_file     cast shadow in Y direction and output to shadow_file (default: no shadow ouput)\n");
    fprintf(stderr,"\t\t-v                 increase verbosity\n");
    fprintf(stderr,"\t\t-q                 quiet\n\n");
}

void InvalidOption(char *opt) {
    fprintf(stderr,"Invalid option: '%s'\n", opt);
    Usage();
    exit(1);
}

int GetOptions(int argc, char** argv) {
    int i = 1;
    while( i<argc ) {
        if ( argv[i][0]=='-' ) {
            switch (argv[i][1]) {
            case 'i':
                info = 1;
                break;
            case 'I':
                info = 2;
                break;
            case 'x':
                switch (argv[i][3]) {
                    case 'a':   // xmax
                        i++;
                        VMax[0] = atof(argv[i]);
                        break;
                    case 'i':   // xmin
                        i++;
                        VMin[0] = atof(argv[i]);
                        break;
                    default:
                        InvalidOption(argv[i]);
                }
                break;
            case 'y':
                switch (argv[i][3]) {
                    case 'a':   // ymax
                        i++;
                        VMax[1] = atof(argv[i]);
                        break;
                    case 'i':   // ymin
                        i++;
                        VMin[1] = atof(argv[i]);
                        break;
                    default:
                        InvalidOption(argv[i]);
                }
                break;
            case 'z':
                switch (argv[i][3]) {
                    case 'a':   // zmax
                        i++;
                        VMax[2] = atof(argv[i]);
                        break;
                    case 'i':   // zmin
                        i++;
                        VMin[2] = atof(argv[i]);
                        break;
                    default:
                        InvalidOption(argv[i]);
                }
                break;
            case 't':
                switch (argv[i][2]) {
                    case 'x':
                        i++;
                        Translate[0] = atof(argv[i]);
                        break; 
                    case 'y':
                        i++;
                        Translate[1] = atof(argv[i]);
                        break; 
                    case 'z':
                        i++;
                        Translate[2] = atof(argv[i]);
                        break; 
                    case ' ':
                    case '\0':
                        i++;
                        sscanf(argv[i],"%f,%f,%f", Translate, Translate+1, Translate+2);
                        break; 
                    default:
                        InvalidOption(argv[i]);
                }
                break;
            case 's':
                switch (argv[i][2]) {
                    case 'x':
                        i++;
                        Scale[0] = atof(argv[i]);
                        break; 
                    case 'y':
                        i++;
                        Scale[1] = atof(argv[i]);
                        break; 
                    case 'z':
                        i++;
                        Scale[2] = atof(argv[i]);
                        break; 
                    case ' ':
                    case '\0':
                        i++;
                        Scale[0] = Scale[1] = Scale[2] = atof(argv[i]);
                        break; 
                    default:
                        InvalidOption(argv[i]);
                }
                break;
                case 'r':
                switch (argv[i][2]) {
                    case 'x':
                        i++;
                        Rotate[0] = PI*atof(argv[i])/180.;
                        break;
                    case 'y':
                        i++;
                        Rotate[1] = PI*atof(argv[i])/180.;
                        break;
                    case 'z':
                        i++;
                        Rotate[2] = PI*atof(argv[i])/180.;
                        break;
                    case ' ':
                    case '\0':
                    default:
                            InvalidOption(argv[i]);
                }
                break;
            case 'm':
                i++;
                Material = argv[i];
                break;
            case 'N':
                switch (argv[i][2] ) {
                    case 'i':
                        InvertNormals = 1;
                        break;
                    case 'c':
                        InvertNormals = 2;
                        break;
                    case 's':
                            InvertNormals = 4;
                            break;
                    default:
                            InvalidOption(argv[i]);
                }           
                break;
            case 'f':
                // Find intersection
                switch (argv[i][2]) {
                    case 'z':
                        findIntersection = 1;
                        break;
                    case 'y':
                        findIntersection = 2;
                        break;
                    default:
                        InvalidOption(argv[i]);
                }
                i++;
                Vertical[0] = atof(argv[i]);
                i++;
                Vertical[findIntersection] = atof(argv[i]);
                break;
            case 'g':
                i++;
                SelectGroup[SelGroups++] = argv[i];
                break;
            case 'o':
                i++;
                SelectObject[SelObjects++] = argv[i];
                break;
            case 'O':
                i++;
                OutputFile = argv[i];
                break;
            case 'S':
                switch (argv[i][2]) {
                    case 'y':
                    case 'Y':
                        upAxis = 'y';
                        break;
                    case 'z':
                    case 'Z':
                        upAxis = 'z';
                        break;
                    default:
                        upAxis = 'y';
                }
                i++;
                ShadowOutputFile = argv[i];
                break;
            case 'c':
                SolidCut = 1;
                break;
            case 'R':
                Relative = 1;
                break;
            case 'v':
                Verbose += 1;
                break;
            case 'q':
                Verbose = -1;
                break;
            case 'n':
                Negate = 1;
                break;
            case 'e':
                Explode = 1;
                break;
            case 'X':
                genTextCoords = 1;
                break;
            default:
                InvalidOption(argv[i]);
            }
        }
        else break;
        i++;
    }
    return i;
}

//  Joins nObjs ObjFiles from ObjSet array.
//  Outputs whole ObjFile in buffer pointed by obj.
//  obj must be pre allocated.
//  If nObjs==1, copies ObjFile pointed by ObjSet to obj.
int JoinObjFiles(int nObjs, ObjFile ObjSet[], ObjFile *obj) {

    ResetObjStats( &(obj->stats) );

    int i;
    for( i=0 ; i<nObjs ; i++ ) {
        obj->stats.verts += ObjSet[i].stats.verts;
        obj->stats.faces += ObjSet[i].stats.faces;
        obj->stats.norms += ObjSet[i].stats.norms;
        obj->stats.texts += ObjSet[i].stats.texts;
        obj->stats.paras += ObjSet[i].stats.paras;
        obj->stats.grps  += ObjSet[i].stats.grps;
        obj->stats.mats  += ObjSet[i].stats.mats;
        obj->stats.objs  += ObjSet[i].stats.objs;
        obj->stats.libs  += ObjSet[i].stats.libs;
        obj->stats.shds  += ObjSet[i].stats.shds;
        
        int j;
        for( j=0 ; j<3 ; j++ ) {
            if ( ObjSet[i].stats.vmin[j] < obj->stats.vmin[j] )
                obj->stats.vmin[j] = ObjSet[i].stats.vmin[j];
            if ( ObjSet[i].stats.vmax[j] > obj->stats.vmax[j] )
                obj->stats.vmax[j] = ObjSet[i].stats.vmax[j];
        }
        for( j=0 ; j<2 ; j++ ) {
            if ( ObjSet[i].stats.tmin[j] < obj->stats.tmin[j] )
                obj->stats.tmin[j] = ObjSet[i].stats.tmin[j];
            if ( ObjSet[i].stats.tmax[j] > obj->stats.tmax[j] )
                obj->stats.tmax[j] = ObjSet[i].stats.tmax[j];
        }
        //PrintObjStats(&(obj->stats));
    }

    //printf("  Total\n");
    //PrintObjStats(&(obj->stats));

    AllocObjFile(obj);
    ObjStats counters;
    ResetObjStats(&counters);

    for( i=0 ; i<nObjs ; i++ ) {
        memcpy(obj->verts+counters.verts, ObjSet[i].verts, ObjSet[i].stats.verts*sizeof(Vert));
        memcpy(obj->norms+counters.norms, ObjSet[i].norms, ObjSet[i].stats.norms*sizeof(Norm));
        memcpy(obj->texts+counters.texts, ObjSet[i].texts, ObjSet[i].stats.texts*sizeof(Text));
        int j;
        for( j=0 ; j<ObjSet[i].stats.faces ; j++ ) {
            obj->faces[counters.faces+j].nodes = ObjSet[i].faces[j].nodes;
            obj->faces[counters.faces+j].Node = Malloc(ObjSet[i].faces[j].nodes, sizeof(FaceNode));
            memcpy(obj->faces[counters.faces+j].Node, ObjSet[i].faces[j].Node, 
                ObjSet[i].faces[j].nodes*sizeof(FaceNode) );
            int n;
            for( n=0 ; n<obj->faces[counters.faces+j].nodes ; n++ ) {
                obj->faces[counters.faces+j].Node[n][0] += counters.verts;
                obj->faces[counters.faces+j].Node[n][1] += counters.texts;
                obj->faces[counters.faces+j].Node[n][2] += counters.norms;
            }
        }
        for( j=0 ; j<ObjSet[i].stats.libs ; j++ ) {
            //printf("merging lib\n");
            //printf("merging lib %d %d: %s\n", i, j, ObjSet[i].libs[j].name);
            obj->libs[counters.libs+j].line = ObjSet[i].libs[j].line + counters.faces;
            obj->libs[counters.libs+j].name = strdup(ObjSet[i].libs[j].name);
        }
        for( j=0 ; j<ObjSet[i].stats.objs ; j++ ) {
            obj->objs[counters.objs+j].line = ObjSet[i].objs[j].line + counters.faces;
            obj->objs[counters.objs+j].name = strdup(ObjSet[i].objs[j].name);
        }
        for( j=0 ; j<ObjSet[i].stats.shds ; j++ ) {
            obj->shds[counters.shds+j].line = ObjSet[i].shds[j].line + counters.faces;
            obj->shds[counters.shds+j].name = strdup(ObjSet[i].shds[j].name);
        }
        for( j=0 ; j<ObjSet[i].stats.grps ; j++ ) {
            obj->grps[counters.grps+j].line = ObjSet[i].grps[j].line + counters.faces;
            obj->grps[counters.grps+j].name = strdup(ObjSet[i].grps[j].name);
        }
        for( j=0 ; j<ObjSet[i].stats.mats ; j++ ) {
            //printf("merging mat %d %d: %s\n", i, j, ObjSet[i].mats[j].name);
            obj->mats[counters.mats+j].line = ObjSet[i].mats[j].line + counters.faces;
            obj->mats[counters.mats+j].name = strdup(ObjSet[i].mats[j].name);
        }

        counters.verts += ObjSet[i].stats.verts;
        counters.norms += ObjSet[i].stats.norms;
        counters.texts += ObjSet[i].stats.texts;
        counters.faces += ObjSet[i].stats.faces;
        counters.libs  += ObjSet[i].stats.libs;
        counters.objs  += ObjSet[i].stats.objs;
        counters.shds  += ObjSet[i].stats.shds;
        counters.grps  += ObjSet[i].stats.grps;
        counters.mats  += ObjSet[i].stats.mats;
    }

//printf("merge done\n");
    return nObjs;
}

void genTextureCoords() {

}



void SaveImageMap(char *OutputFile, ObjFile *obj) {
    int nVerts = obj->stats.verts;
    float side = round(sqrt(nVerts));
    //int iSide = (int)side;
    
    printf("im size: %d %f\n", nVerts, side);
    // Find desired imagem geometry
    int xi, yi=0;
    for( xi=(int)side ; xi>=2 ; xi-- ) {
        yi = nVerts / xi;
        printf("%d*%d=%d  <> %d\n", xi, yi, xi*yi, nVerts);
        if ( yi*xi == nVerts ) break;
    }
    printf("xI*yI=%d*%d=%d  <> %d\n", xi, yi, xi*yi, nVerts);
    if ( xi*yi!=nVerts ) {
        printf("Erro de amostragem\n");
        return;
    }
    printf("xi yi %d %d\n", xi, yi);
    int maxS, minS;
    if ( xi>yi ) {
        maxS = xi;
        minS = yi;
    }
    else {
        maxS = yi;
        minS = xi;
    }
    int i;
    float delta[3];
    // float minV[3];
    // float maxV[3];
    float minD =  1.e8;
    int   minI = -1;
    float maxD = -1.e8;
    int   maxI = -1;
    int   medI = -1;
    for( i=0 ; i<3 ; i++ ) {
        delta[i] = obj->stats.vmax[i]-obj->stats.vmin[i];
        if (delta[i]<minD) {
            minD = delta[i];
            minI = i;       // UP Coord index
        }
        if (delta[i]>maxD) {
            maxD = delta[i];
            maxI = i;       // Future image X coord index. (Column number)
        }
        printf("i: %d, delta: %f, minD: %f, maxD:%f\n", i, delta[i], minD, maxD);
    }
    printf("im size: %f %f %f\n", delta[0], delta[1], delta[2]);
    printf("minS: %d maxS: %d, maxI: %d\n",minS, maxS, maxI);
    for( i=0 ; i<3 ; i++ ) {
        if ( i==minI )
            delta[i] /= 256.;
        else {
            if ( i==maxI )
                delta[i] /= maxS;   // Image X resolution
            else {
                medI = i;           // Image Y coord index (Number of row)
                delta[i] /= minS;   // Image Y resolution
            }
        }
    }
    assert(medI>=0);
    printf("im deltas: %f %f %f\n", delta[0], delta[1], delta[2]);
    FILE *fout = fopen(OutputFile, "w");
    fprintf(fout, "P5 %d %d 255\n", minS, maxS); 
    
    int flag[nVerts];
    for( i=0; i<nVerts ; i++ )
        flag[i] = 0;
    for( xi=0 ; xi<maxS ; xi++ ) {
        double x = obj->stats.vmin[maxI] + yi*delta[maxI];
        for( yi=0 ; yi<minS ; yi++ ) {
            double y = obj->stats.vmin[medI] + yi*delta[medI];
            int idx= getVertexIdx(obj, maxI, medI, x, y);
            if ( flag[idx]>0 ) printf("Erro de flag\n");
            flag[idx]++;
            float z = obj->verts[idx][minI];
            int zq = (int)((z-obj->stats.vmin[minI]) / delta[minI]);
            fputc(zq, fout);
        }
    }
    fclose(fout);
}
    
int main(int argc, char **argv) {

    int fi = GetOptions(argc,argv);
    if ( Verbose>=0 )
        fprintf(stderr, "# ObjTool\n# Compile date: %s\n", CDATE);
    
    //printf("# fi = %d, argc=%d\n",fi,argc);

    if ( fi>=argc ) {
        fprintf(stderr, "Al least one input file must be specified.\n");
        Usage();
        exit(1);
    }

    // Alloc and load obj files
    int NObjects = argc-fi;
    ObjFile *ObjSet;
    ObjSet = Malloc(NObjects, sizeof(ObjFile));

    int i = 0;
    while ( fi<argc) {
        if ( Verbose>0 ) fprintf(stderr,"Loading ObjFile '%s'\n", argv[fi]);
        LoadObjFile(argv[fi], ObjSet+i);
        SetUseCounters(ObjSet+i);
        i++;
        fi++;
    }

    // All files got.
    // Join all geometry into one ObjFile struct
    ObjFile obj;    
    JoinObjFiles(NObjects, ObjSet, &obj);
    
    SetUseCounters(&obj);
    if ( info ) {
        if ( Verbose>2 )
            PrintFullObjStats(&obj);
        else 
            if ( Verbose>1 )
                PrintLongObjStats(&obj);
            else
                PrintObjStats(&(obj.stats));
        exit(0);
    }
    if ( Verbose>1 ) {
        if ( Verbose>5 )
            PrintFullObjStats(&obj);
        else
            PrintObjStats(&(obj.stats));
    }
    
    if ( findIntersection==1 ) {
        double z = findVerticalIntersection(&obj, Vertical);
        printf("%lf\n", z);
        exit(0);
    }
    if ( findIntersection==2 ) {
        double y = findVerticalIntersection(&obj, Vertical);
        printf("%lf\n", y);
        exit(0);
    }
    FilterVerts( &obj, VMin, VMax );
    if ( Verbose>2 ) printf("Vertexes Filtered\n");
    //SetIndexs(&obj);
    
    CleanFaces(&obj);
    
    SetUseCounters(&obj);
    ProcVerts( &obj );
    if ( Verbose>2 ) printf("Vertexes processed\n");
    SetIndexs(&obj);
    if ( Verbose>2 ) printf("Indexes set\n");
    if ( Explode )
        ExplodeOutputFile(OutputFile, &obj);
    else {
        if ( Verbose>8 ) printf("Saving to file\n");
        if ( OutputFile && strrstr(OutputFile, ".pgm") ) 
            SaveImageMap(OutputFile, &obj);
        else
            SaveObjFile(OutputFile, &obj);
    }
    if ( ShadowOutputFile ) {
        if ( Verbose>3 ) printf("Creating Shadow Obj\n");
        //ObjFile *Shadow = CreateShadowObj(&obj);
        ObjFile *Shadow = CreateShadowObj(&obj);
        SetUseCounters(Shadow);
        SetIndexs(Shadow);
        if ( Verbose>4 ) PrintObjStats(&(Shadow->stats));
        SaveObjFile(ShadowOutputFile, Shadow);
        if ( Verbose>4 ) printf("saved\n");
        FreeObjFile(Shadow);
        if ( Verbose>4 ) printf("Shadow file freed\n");
    }
    if ( Verbose>4 ) printf("vai limpar verts\n");
    FreeObjFile(&obj);
    
    if ( Verbose>4 ) printf("limpou verts\n");
    return 0;
}
        


