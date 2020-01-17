/*
 *	ObjTool.h
 *
 *	Command line toolset to work with Wavefront 3D Obj files.
 *
 *	Uses OBJ file format as defined in http://paulbourke.net/dataformats/obj/
 * 	Not all OBJ file format specs are implemented.
 */

#define VERSION 0.13.4

#define _GNU_SOURCE
 	
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "linmath.h"

#define PI M_PI

#define MAX_LINE_LEN (5880)
#define MAX_FACE_NODES (380)

#define MAX_FLOAT (1.0e32)

#define NULL_IDX (-999999)

#define MEM_META_INFO_SIZE (16)
#define MEM_TAG (0x7a76b532)

typedef struct ObjStatsStruct {
	int verts;
	int faces;
	int norms;
	int texts;
	int paras;
	int objs;
	int grps;
	int mats;
	int libs;
	int shds;
	float vmin[3];
	float vmax[3];
	float tmin[2];
	float tmax[2];
} ObjStats;

typedef struct UseCountStruct {
	int *verts;
	int *norms;
	int *texts;
/*
	int *mats;
	int *shds;
	int *grps;
	int *objs;  */
} UseCounters;

typedef int FaceNode[3];	// VeterxIdx TextureCoordIdx NormalIdx

typedef float Vert[3];		// x y z

typedef float Vec2[2];
typedef float Vec3[3];
//typedef FaceNode *Face;

typedef struct FaceStruct {
	int nodes;		// Number of Vertexs
	FaceNode *Node;
} Face;

typedef struct ObjGroupStruct {
	int line;
	char *name;
} ObjGroup;

typedef float Norm[3];		// x y z
typedef float Text[2];		// u v

typedef struct ObjFileStruct {
	Vert *verts;
	Face *faces;
	Norm *norms;
	Text *texts;
	Vert *paras;
	ObjGroup *objs;
	ObjGroup *grps;
	ObjGroup *mats;
	ObjGroup *libs;
	ObjGroup *shds;
	ObjStats stats;
	UseCounters counts;
	UseCounters order;
} ObjFile;

// Globals
extern short int findIntersection;
extern short int Verbose;

/*
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
*/

#define SinF(x) ( sinf(x) )
#define CosF(x) ( cosf(x) )

#define AngNormalize(x) ( (x)>PI ? (x-2*PI) : ( (x)<-PI ? (x+2*PI) : (x) ) )

int JoinObjFiles(int nObjs, ObjFile ObjSet[], ObjFile *obj);
void SetUseCounters(ObjFile *obj);
void SetIndexs(ObjFile *obj);
int getVertexIdx(ObjFile *obj, int xCoord, int yCoord, double x, double y);
double findVerticalIntersection(ObjFile *obj, Vert coords);

