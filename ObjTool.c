/*
 *	ObjTool.c
 *
 *	Command line toolset to work with Wavefront 3D Obj files.
 *
 *	Uses OBJ file format as defined in http://paulbourke.net/dataformats/obj/
 * 	Not all OBJ file format specs are implemented.
 */

#define VERSION 0.13.2
#define _GNU_SOURCE
 	
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

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

Vert VMin = { -1e30, -1e30, -1e30 };
Vert VMax = {  1e30,  1e30,  1e30 };
Vert Translate = { 0., 0., 0. };
Vert Scale = { 1., 1., 1. };
Vert Rotate = { 0., 0., 0. };
char *Material = NULL;
char *OutputFile = NULL;
char *ShadowOutputFile = NULL;
short int SolidCut = 0;
short int Relative = 0;
short int Verbose = 0;
short int info = 0;
short int Negate = 0;
short int Explode = 0;
char *SelectGroup[20];
char *SelectObject[20];
int SelGroups = 0;
int SelObjects = 0;;

#define SinF(x) ( sinf(x) )
#define CosF(x) ( cosf(x) )

int JoinObjFiles(int nObjs, ObjFile ObjSet[], ObjFile *obj);
void SetUseCounters(ObjFile *obj);
void SetIndexs(ObjFile *obj);

void *Malloc(int n, size_t dim) {
	long size = n*dim;
	void *buf = malloc(size+MEM_META_INFO_SIZE);
	if ( ! buf ) {
		fprintf(stderr,"Cannot alloc memory\n");
		exit(3);
	}
	*((long *)buf) = MEM_TAG;
	*((long *)(buf+8)) = (long)(size+MEM_META_INFO_SIZE);
	memset(buf+MEM_META_INFO_SIZE,0,size);
	return buf+MEM_META_INFO_SIZE;
}

void Free(void *p) {
	long *tag_p = p-MEM_META_INFO_SIZE;
	if ( *tag_p == MEM_TAG ) {
		*tag_p = 0L;
		free(p-MEM_META_INFO_SIZE);
	}
//printf("Memory freed\n");
//else printf("##Memory not freed\n");
}

char *StrDup(char *strin) {
	char *strout = strdup(strin);
	char *strptr = strout;
	while( *strptr ) {
		if ( *strptr=='\n' ) break;
		if ( *strptr=='\r' ) break;
		if ( *strptr=='\t' ) break;
		strptr++;
	}
	*strptr = 0;
	return strout;
}

int ReadVec2f(char *str, float *data) {
	int res;
	res = sscanf(str,"%f %f", data+0, data+1);
	return res;
}

int ReadVec3f(char *str, float *data) {
	int res;
	res = sscanf(str,"%f %f %f", data+0, data+1, data+2);
	return res;
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



void ResetObjStats(ObjStats *os) {
	os->verts = 0;
	os->faces = 0;
	os->norms = 0;
	os->texts = 0;
	os->paras = 0;
	os->grps  = 0;
	os->mats  = 0;
	os->objs  = 0;
	os->libs  = 0;
	os->shds  = 0;
	os->vmin[0] =  MAX_FLOAT;
	os->vmin[1] =  MAX_FLOAT;
	os->vmin[2] =  MAX_FLOAT;
	os->vmax[0] = -MAX_FLOAT;
	os->vmax[1] = -MAX_FLOAT;
	os->vmax[2] = -MAX_FLOAT;

	os->tmin[0] =  MAX_FLOAT;
	os->tmin[1] =  MAX_FLOAT;
	os->tmax[0] = -MAX_FLOAT;
	os->tmax[1] = -MAX_FLOAT;

}

void PrintObjStats(ObjStats *os) {
printf("printing stats\n");
	printf("Vertices: %d\n", os->verts);
	printf("  Min: %.3f %.3f %.3f\n", os->vmin[0], os->vmin[1], os->vmin[2]);
	printf("  Max: %.3f %.3f %.3f\n", os->vmax[0], os->vmax[1], os->vmax[2]);
	printf("Normals: %d\n", os->norms);
	printf("Texture Coords: %d\n", os->texts);
	if ( os->texts > 0 ) {
	    printf("  Min: %.3f %.3f\n", os->tmin[0], os->tmin[1]);
	    printf("  Max: %.3f %.3f\n", os->tmax[0], os->tmax[1]);
	}
	printf("Faces: %d\n", os->faces);
	printf("Objects: %d\n", os->objs);
	printf("ObjGroups: %d\n", os->grps);
	printf("UseMtl: %d\n", os->mats);
	printf("MtlLib:	%d\n", os->libs);
	printf("ShadGrp:	%d\n", os->shds);
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

int GetObjStats(char *fname, ObjStats *os) {
	int nl = 0;
	int i;
	float vtx[3];
	ResetObjStats(os);

	FILE  *fobj;
	fobj = fopen(fname,"r");
	if ( ! fobj ) {
		fprintf(stderr,"Cannot open input file '%s'\n",fname);
		return 0;
	}
	char line[MAX_LINE_LEN+1];
	while( fgets(line,MAX_LINE_LEN,fobj) ) {
		nl++;
		char *lptr = line;
		while( *lptr==' ' || *lptr=='\t' )	lptr++;
		switch (*lptr) {
			case 0:
			case '\n':
			case '\r':
			case '#': continue;
			case 'f':
			case 'F':
				os->faces += 1;
				continue;
			case 'm':
				os->libs += 1;
				continue;
			case 'o':
				os->objs +=1;
				continue;
			case 'g':
				os->grps +=1;
				continue;
			case 'u':
				os->mats += 1;
				continue;
			case 's':
				os->shds += 1;
				continue;
			case 'v':
				switch ( *(lptr+1) ) {
					case ' ':
						ReadVec3f(lptr+2,vtx);
						for( i=0 ; i<3 ; i++ ) {
						    if ( vtx[i]<os->vmin[i] )	os->vmin[i] = vtx[i];
						    if ( vtx[i]>os->vmax[i] )	os->vmax[i] = vtx[i];
						}
						os->verts += 1;
						break;
					case 'n':
						os->norms += 1;
						break;
					case 't':
						ReadVec3f(lptr+3,vtx);
						for( i=0 ; i<2 ; i++ ) {
						    if ( vtx[i]<os->tmin[i] )	os->tmin[i] = vtx[i];
						    if ( vtx[i]>os->tmax[i] )	os->tmax[i] = vtx[i];
						}
						os->texts += 1;
						break;
					case 'p':
						os->paras += 1;
						break;
					default:
						fprintf(stderr,"Invalid line %d:%s\n", nl, line);
						break;
				}
				continue;
			default:
				fprintf(stderr,"Invalid line %d:%s\n", nl, line);
				break;
		}
		break;
	}
	fclose(fobj);
	//printf("%d %d %d %d %d\n",os->faces,os->verts,os->norms,os->texts,os->paras);
	return 1;
}


char *ReadFaceNode(char *lptr, FaceNode node) {
        // Read %d/%d/%d set
	// VeterxIdx/TextureCoordIdx/NormalIdx
	while( *lptr==' ' || *lptr=='\t' ) lptr++;
	int i;
	for( i=0 ; i<3 ; i++ ) node[i] = NULL_IDX;
	for( i=0 ; i<3 ; i++ ) {
		if ( *lptr=='/' ) {
			lptr++;
			continue;
		}
		node[i] = atoi(lptr);
		if ( node[i]==0 )
			fprintf(stderr,"Invalid Node Coord %d, %d: %s\n", i, node[i], lptr);
		while( *lptr!=' ' && *lptr!='\t' && *lptr!='/' && *lptr ) lptr++;
		if ( *lptr!='/' ) return lptr;
		lptr++;
	}	
	return lptr;
}

int VerticalCrossSegment(Vert *p, Vert v1, Vert v2) {
	// Test if vector UP from p intersect segment v1,v2
	// p,v1,v2 are 2d points

	double dx = v2[0] - v1[0];
	if ( fabs(dx)<0.00001 ) 	// vertical
		return 0;
	int res = 0;
	double dy = v2[1] - v1[1];
	double m = dy / dx;		
	double y0 = v1[1]-m*(v1[0]-(*p)[0]);
	if ( y0>(*p)[1] ) {
		if ( dy>0 ) {
			if ( y0-v1[1]>.1 && y0-v2[1]<-.1 ) res = 1;
		} else {
			if ( y0-v1[1]<-.1 && y0-v2[1]>.1 ) res = 1;
		}
	}
	//printf( "        v1 %f %f, v2 %f %f, res:%d\n", v1[0],v1[1],v2[0],v2[1],res );
	return res;	
}

int PointInPolygon(int v, Vert *verts, FaceNode *node, int n) {
	// Test if vertex p is outside polygon vert[]
	Vert *p = verts+v;
	int cross = 0; 
	if ( n<3 ) return cross;	// polygon must be closed	
	int i;
	for( i=1 ; i<n ; i++ )
		cross += VerticalCrossSegment(p,verts[node[i-1][0]],verts[node[i][0]]);
	cross += VerticalCrossSegment(p,verts[node[0][0]],verts[node[n-1][0]]);
//	printf("    PointInPoly n:%d vert %d, point p %f %f: cross: %d\n",n,v,(*p)[0],(*p)[1], cross);
	return cross % 2;
}

ObjFile *CreateShadowObj_v1(ObjFile *obj) {
  //  printf("Creating shadow\n");
    int i, si;	// face indexs
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
  
    shadow->stats.norms = 1;	// Only one norm, facing up.

    shadow->stats.texts = obj->stats.texts;

	shadow->norms = Malloc(1, sizeof(Norm));
	shadow->norms[0][0] = 0.;
	shadow->norms[0][1] = 0.;
	shadow->norms[0][2] = 1.;	// Up vector

	
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
		    fprintf(stderr,"Vertice %d not found on face %d, node %d/%d.\n", v, i, n, obj->faces[i].nodes);
		    break;
		}
		// Test if vertex is outside polygon already draw
		int pp;
		for( pp=0 ; pp<shadow->stats.faces ; pp++ ) {
		    if ( PointInPolygon( v-1, shadow->verts, shadow->faces[pp].Node, shadow->faces[pp].nodes ) ) {
		    	inside += 1;
		    	break;
		    }
		}
	    }
	    if ( inside<=n ) {	// Add face
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
fprintf(stderr,"Creating shadow\n");
  // nova versao por convex hull
  //  printf("Creating shadow\n");
    int i, si;	// face indexs
/*
    double max[2];
    double min[2];
    for( i=0 ; i<2 ; i++ ) {
	max[i] = -MAX_FLOAT;
	min[i] =  MAX_FLOAT;
    }
*/
    ObjFile *shadow = Malloc(1,sizeof(ObjFile));

    shadow->stats.verts = 0;
    shadow->verts = Malloc(obj->stats.verts, sizeof(Vert));

    shadow->stats.faces = 0;
  
    shadow->stats.norms = 1;	// Only one norm, facing up.

    shadow->stats.texts = 0;

	shadow->norms = Malloc(1, sizeof(Norm));
	shadow->norms[0][0] = 0.;
	shadow->norms[0][1] = 0.;
	shadow->norms[0][2] = 1.;	// Up vector

	
	shadow->texts = 0;
	shadow->grps = 0 ;
	shadow->mats = 0;
	shadow->objs = 0;
	shadow->libs = 0;
	shadow->shds = 0;

	shadow->counts.verts = 0;
	shadow->counts.norms = 0;
	shadow->counts.texts = 0;

	shadow->order.verts = obj->order.verts;
	shadow->order.texts = obj->order.texts;
	shadow->order.norms = obj->order.norms ;

   // Find z min
   float zmin = MAX_FLOAT;
   for( i=0 ; i<obj->stats.verts ; i++ )
	if ( obj->verts[i][2]<zmin )
		zmin = obj->verts[i][2];


   // Find leftest vertex
   float xmin = MAX_FLOAT;
   float ymin = MAX_FLOAT;
   for( i=0 ; i<obj->stats.verts ; i++ )
	if ( obj->verts[i][0]<xmin ) {
		xmin = obj->verts[i][0];
		ymin = obj->verts[i][1];
	}
   
   // Find next vertex, with greatest slope
   float xmax, ymax;
   float xp = xmin;
   float yp = ymin;
   while (1) {
	double mmax = -MAX_FLOAT;
	ymax = xmax = -MAX_FLOAT;
	for( i=0 ; i<obj->stats.verts ; i++ ) {
	    float x = obj->verts[i][0];
	    float y = obj->verts[i][1];
	    if ( x<=xp ) continue;
	    double m = (y-yp)/(x-xp);
	    if ( m>mmax ) {
		mmax = m;
		xmax = x;
		ymax = y;
	    }
	}
printf("mmax: %f\n",mmax);
	if ( mmax<=-MAX_FLOAT/2. ) break;
	// new border vertex found
	xp = xmax;
	yp = ymax;
	si = shadow->stats.verts;
	shadow->verts[si][0] = xp;
	shadow->verts[si][1] = yp;
	shadow->stats.verts = si+1;
    }
   
   // smachdown zz coords to zmin
   for( i=0 ; i<shadow->stats.verts ; i++ )
	shadow->verts[i][2] = zmin;
printf("shadow done\n");
  return shadow;
    for( i=0 ; i<obj->stats.faces ; i++ ) {
	    int n;
	    short int inside = 0;
//printf("Testing face %d, %d nodes\n", i, obj->faces[i].nodes);
	    for( n=0 ; n<obj->faces[i].nodes ; n++ ) {
		//printf(" Testing node %d\n",n);
		int v = obj->faces[i].Node[n][0];
		if ( v<1 || v>obj->stats.verts ) {
		    fprintf(stderr,"Vertice %d not found on face %d, node %d/%d.\n", v, i, n, obj->faces[i].nodes);
		    break;
		}
		// Test if vertex is outside polygon already draw
		int pp;
		for( pp=0 ; pp<shadow->stats.faces ; pp++ ) {
		    if ( PointInPolygon( v-1, shadow->verts, shadow->faces[pp].Node, shadow->faces[pp].nodes ) ) {
		    	inside += 1;
		    	break;
		    }
		}
	    }
	    if ( inside<=n ) {	// Add face
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

void AllocObjFile(ObjFile *obj) {
	assert(obj!=NULL);
	obj->faces = Malloc(obj->stats.faces, sizeof(Face));
	obj->verts = Malloc(obj->stats.verts, sizeof(Vert));
	obj->norms = Malloc(obj->stats.norms, sizeof(Norm));
	obj->texts = Malloc(obj->stats.texts, sizeof(Text));
	obj->grps  = Malloc(obj->stats.grps,  sizeof(ObjGroup));
	obj->mats  = Malloc(obj->stats.mats,  sizeof(ObjGroup));
	obj->objs  = Malloc(obj->stats.objs,  sizeof(ObjGroup));
	obj->libs  = Malloc(obj->stats.libs,  sizeof(ObjGroup));
	obj->shds  = Malloc(obj->stats.shds,  sizeof(ObjGroup));

	obj->counts.verts = Malloc(obj->stats.verts, sizeof(int));
	obj->counts.norms = Malloc(obj->stats.norms, sizeof(int));
	obj->counts.texts = Malloc(obj->stats.texts, sizeof(int));

	obj->order.verts = Malloc(obj->stats.verts, sizeof(int));
	obj->order.texts = Malloc(obj->stats.texts, sizeof(int));
	obj->order.norms = Malloc(obj->stats.norms, sizeof(int));
}

void LoadObjFile(char *fname, ObjFile *obj) {
	int nl = 0;
	int i;
	int res = GetObjStats( fname, &(obj->stats) );
	if ( ! res ) return;

	AllocObjFile(obj);
	
//printf("face %ld, Vert %ld  Norm %ld, Text %ld\n", sizeof(Face), sizeof(Vert), sizeof(Norm), sizeof(Text) );
	ObjStats counters;
	ResetObjStats(&counters);
	FILE  *fobj;
	fobj = fopen(fname,"r");
	char line[MAX_LINE_LEN+1];
	while( fgets(line,MAX_LINE_LEN,fobj) ) {
		nl++;
//printf("line %d(%d):%s.\n",nl,strlen(line),line);
		char *lptr = line;
		while( *lptr==' ' || *lptr=='\t' )	lptr++;
		switch (*lptr) {
			case 0:
			case '\n':
			case '\r':
			case '#': continue;
			case 'f':
			case 'F':
				if ( counters.faces>=obj->stats.faces ) {
					fprintf(stderr,"Error: %d/%d faces read\n",counters.faces,obj->stats.faces);
					exit(7);
				}
				lptr++;
				int fv = 0;
				FaceNode tmp_face_nodes[MAX_FACE_NODES]; 
				while ( lptr && *lptr && *lptr!='\n' && *lptr!='\r' ) {
					lptr = ReadFaceNode(lptr,tmp_face_nodes[fv]);
					while( *lptr==' ' || *lptr=='\t' )	lptr++;
					if ( tmp_face_nodes[fv][0]<0 )
						tmp_face_nodes[fv][0] += counters.verts+1;
					if ( tmp_face_nodes[fv][1]!=NULL_IDX && tmp_face_nodes[fv][1]<0 )
						tmp_face_nodes[fv][1] += counters.texts+1;
					if ( tmp_face_nodes[fv][2]!=NULL_IDX && tmp_face_nodes[fv][2]<0 )
						tmp_face_nodes[fv][2] += counters.norms+1;
					fv++;
					if ( fv>=MAX_FACE_NODES ) {
					    fprintf(stderr,"Too many nodes in face %d.\n", counters.faces);
					    break;
					}
				}
				obj->faces[counters.faces].nodes = fv;
				obj->faces[counters.faces].Node = Malloc(fv,sizeof(FaceNode));
//printf("size=%d\n", fv*sizeof(FaceNode) );
				memcpy(obj->faces[counters.faces].Node,tmp_face_nodes,fv*sizeof(FaceNode));
				counters.faces += 1;
				continue;
			case 'm':		// mtllib
				obj->libs[counters.libs].name = StrDup(lptr+7);
				obj->libs[counters.libs].line = counters.faces;
				counters.libs += 1;
				continue;
			case 'o':
				obj->objs[counters.objs].name = StrDup(lptr+2);
				obj->objs[counters.objs].line = counters.faces;
				counters.objs += 1;
				continue;
			case 's':
				obj->shds[counters.shds].name = StrDup(lptr+2);
				obj->shds[counters.shds].line = counters.faces;
				counters.shds += 1;
				continue;
			case 'g':
				obj->grps[counters.grps].name = StrDup(lptr+2);
				obj->grps[counters.grps].line = counters.faces;
				counters.grps += 1;
				continue;
			case 'u':		// usemtl
				obj->mats[counters.mats].name = StrDup(lptr+7);
				obj->mats[counters.mats].line = counters.faces;
				counters.mats += 1;
				continue;
			case 'v':
				switch ( *(lptr+1) ) {
					case ' ':		// Vertex
						ReadVec3f(lptr+2,obj->verts[counters.verts]);
						for( i=0 ; i<3 ; i++ ) {
						    if ( obj->verts[counters.verts][i]<obj->stats.vmin[i] )
							obj->stats.vmin[i] = obj->verts[counters.verts][i];
						    if ( obj->verts[counters.verts][i]>obj->stats.vmax[i] )
							obj->stats.vmax[i] = obj->verts[counters.verts][i];
						}
						counters.verts += 1;
						break;
					case 'n':		// Normal
						ReadVec3f(lptr+2,obj->norms[counters.norms]);
						counters.norms += 1;
						break;
					case 't':		// Texture Vertex
						ReadVec2f(lptr+2,obj->texts[counters.texts]);
						for( i=0 ; i<2 ; i++ ) {
						    if ( obj->texts[counters.texts][i]<obj->stats.tmin[i] )
							obj->stats.tmin[i] = obj->texts[counters.texts][i];
						    if ( obj->texts[counters.texts][i]>obj->stats.tmax[i] )
							obj->stats.tmax[i] = obj->texts[counters.texts][i];
						}
						counters.texts += 1;
						//os->texts += 1;
						break;
					case 'p':
						//os->paras += 1;
						break;
					default:
						fprintf(stderr,"Invalid line %d:%s\n", nl, line);
						break;
				}
				continue;
			default:
				fprintf(stderr,"Invalid line %d:%s\n", nl, line);
				break;
		}
		break;
	}
	fclose(fobj);


}

void SaveObjFile(char *fname, ObjFile *obj) {
	FILE *fout;
	fout = stdout;
	int i;
	int nv=0, nt=0, nn=0;
	if ( fname && *fname ) {
		//printf("vai abrir ficheiro %s\n",fname);
		fout = fopen(fname,"w");
		if ( ! fout ) {
		    fprintf(stderr,"Cannot open output file '%s'\n", fname);
		    exit(5);
		}
	        //printf("    abriu\n");
	}
	for( i=0 ; i<obj->stats.verts ; i++ ) 
	    if ( obj->counts.verts[i] > 0 ) {
		fprintf(fout,"v %f %f %f\n", obj->verts[i][0],obj->verts[i][1],obj->verts[i][2]);
		nv++;
	    }
	//printf("    gravou verts\n");
	for( i=0 ; i<obj->stats.texts ; i++ ) 
	    if ( obj->counts.texts[i] > 0 ) {
		fprintf(fout,"vt %f %f %f\n", obj->texts[i][0],obj->texts[i][1],obj->texts[i][2]);
		nt++;
	    }

	//printf("    gravou text verts\n");
	for( i=0 ; i<obj->stats.norms ; i++ ) 
	    if ( obj->counts.norms[i] > 0 ) {
		fprintf(fout,"vn %f %f %f\n", obj->norms[i][0],obj->norms[i][1],obj->norms[i][2]);
		nn++;
	    }
	//printf("    gravou norms\n");
	int gidx=0,oidx=0,midx=0,lidx=0;
	for( i=0 ; i<obj->stats.faces ; i++ ) {
	    int n, ci;

	    while ( lidx<obj->stats.libs && i==obj->libs[lidx].line ) {
		    fprintf(fout,"mtllib %s\n",obj->libs[lidx].name);
		    lidx++;
	    }
	    while ( oidx<obj->stats.objs && i==obj->objs[oidx].line ) {
		for( ci=obj->objs[oidx].line ; oidx<obj->stats.objs-1 && ci<obj->objs[oidx+1].line ; ci++ ) 
		    if ( obj->faces[ci].nodes > 0 ) {
			fprintf(fout,"o %s\n",obj->objs[oidx].name);
			break;
		    }
		oidx++;
	    }
	    while ( gidx<obj->stats.grps && i==obj->grps[gidx].line ) {
		for( ci=obj->grps[gidx].line ; gidx<obj->stats.grps-1 && ci<obj->grps[gidx+1].line ; ci++ )
		    if ( obj->faces[ci].nodes > 0 ) {
		    	fprintf(fout,"g %s\n",obj->grps[gidx].name);
			break;
		    }
		gidx++;
	    }
	    while ( midx<obj->stats.mats && i==obj->mats[midx].line ) {
		for( ci=obj->mats[midx].line ; midx<obj->stats.mats-1 && ci<obj->mats[midx+1].line ; ci++ ) 
		    if ( obj->faces[ci].nodes > 0 ) {
			fprintf(fout,"usemtl %s\n",obj->mats[midx].name);
			break;
		    }
		midx++;
	    }

	    if ( obj->faces[i].nodes <= 0 ) continue;  // face erased
	    fprintf(fout,"f "); 
	    for( n=0 ; n<obj->faces[i].nodes ; n++ ) {
		fprintf(fout,"%d", VertIndex(obj,obj->faces[i].Node[n][0],nv) );
		fprintf(fout,"/");
		
		if ( obj->faces[i].Node[n][1]!=NULL_IDX )
			fprintf(fout,"%d", TextIndex(obj,obj->faces[i].Node[n][1],nt) );
		fprintf(fout,"/");

		if ( obj->faces[i].Node[n][2]!=NULL_IDX )
			fprintf(fout,"%d", NormIndex(obj,obj->faces[i].Node[n][2],nn) );

		fprintf(fout," ");
	    }
	    fprintf(fout,"\n");
	}
	if ( fout!=stdout ) fclose(fout); 
}


void ExplodeOutputFile(char *OutputFile, ObjFile *obj) {

	// Split ObjFile *obj into parts.
	// Each part is saved to file named OutputFile_n_imat_iobj_igrp.obj

	if ( ! OutputFile ) {
	    fprintf(stderr,"Explode requires the output file name using the -o option");
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
printf("    Adding face %d\n", iface);
	    ObjFile *obj_part;
	    obj_part = Malloc(1,sizeof(ObjFile));

	    JoinObjFiles(1,obj,obj_part);
	
	    SetUseCounters(obj_part);
	    
	    for( i=0 ; i<iface ; i++ )
		obj_part->faces[i].nodes = 0;
	    int part_cut = 0;
	    while ( ! part_cut ) {
		while ( imat<obj_part->stats.mats && iface==obj_part->mats[imat].line ) {
		    imat++;
		    part_cut = 1;
	    	}
		while ( iobj<obj_part->stats.objs && iface==obj_part->objs[iobj].line ) {
		    iobj++;
		    part_cut = 1;
	    	}
		while ( igrp<obj_part->stats.grps && iface==obj_part->grps[igrp].line ) {
		    igrp++;
		    part_cut = 1;
	    	}
		iface++;
	    }
	    for( i=iface ; i<obj_part->stats.faces ; i++ )
		obj_part->faces[i].nodes = 0;

	    SetUseCounters(obj_part);
	    SetIndexs( obj_part );

	    snprintf(fname,MAX_FILE_NAME_LEN,"%s_%u_%d_%d_%d.obj", OutputFile, fnum, imat, iobj, igrp);
printf("  Saving PartFile '%s'\n",fname); 	    
	    SaveObjFile(fname, obj_part );

	    Free(obj_part);

	    fnum++;  // next file partition

	} while ( iface < obj->stats.faces );	

}

void FreeObjFile(ObjFile *obj) {
    int i;
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

    Free(obj->faces); 
    Free(obj->verts); 
    Free(obj->norms); 
    Free(obj->texts);
    Free(obj->grps);  
    Free(obj->mats);
    Free(obj->objs);
    Free(obj->shds);
    Free(obj->libs);
    Free(obj->counts.verts);
    Free(obj->counts.norms);
    Free(obj->counts.texts);
    Free(obj->order.verts);
    Free(obj->order.texts);
    Free(obj->order.norms);
}


void SetUseCounters(ObjFile *obj) {
	// Counts uses of each vertex, norm and texture coord.
	// results are stored in obj->counts structure
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
		    fprintf(stderr,"Vertice %d not found on face %d, node %d/%d.\n", v, i, n, obj->faces[i].nodes);
		    continue;
		}
		obj->counts.verts[v-1] += 1;
		v = obj->faces[i].Node[n][1];
		if ( v!=NULL_IDX ) {
		    if ( v<1 || v>obj->stats.texts ) {
		    	fprintf(stderr,"Texture Coord %d not found on face %d, node %d.\n", v, i, n);
		    	continue;
		    }
		    obj->counts.texts[v-1] += 1;
		}
		v = obj->faces[i].Node[n][2];
		if ( v!=NULL_IDX ) {
		    if ( v<1 || v>obj->stats.norms ) {
		    	fprintf(stderr,"Normal %d not found on face %d, node %d.\n", v, i, n);
		    	continue;
		    }
		    obj->counts.norms[v-1] += 1;
		}
	    }
	} 
}



void CleanFaces(ObjFile *obj) {
	int i, n;
	int midx = 0, oidx = 0, gidx = 0;	// current indexes.
						// 0(zero) means material/object/group not yet defined
						// Faces without defined material/object/group always fail the test. 
	for( i=0 ; i<obj->stats.faces ; i++ ) {
	    // Remove face if material matches/does not match selected Material
	    if ( Material!=NULL ) {
		while ( midx<obj->stats.mats && i==obj->mats[midx].line )     midx++;
		if ( Negate ) {
			if ( midx && strcmp(Material,obj->mats[midx-1].name)==0 ) {
			    obj->faces[i].nodes = 0;
			    continue;
			}
		}
		else {
			if ( ! midx || strcmp(Material,obj->mats[midx-1].name)!=0 ) {
			    obj->faces[i].nodes = 0;
			    continue;
			}
		}
	    }
	    if ( SelObjects>0 ) {
	        // Remove face if object matches/does not match selected Object
		while ( oidx<obj->stats.objs && i==obj->objs[oidx].line )     oidx++;
		// Current Test implementation is limited to 1 object ( SelectObject[0] )
		if ( Negate ) {
			if ( oidx && strcmp(SelectObject[0],obj->objs[oidx-1].name)==0 ) {
			    obj->faces[i].nodes = 0;
			    continue;
			}
		}
		else {
			if ( ! oidx || strcmp(SelectObject[0],obj->objs[oidx-1].name)!=0 ) {
			    obj->faces[i].nodes = 0;
			    continue;
			}
		}
	    }
	    if ( SelGroups>0 ) {
	        // Remove face if group matches/does not match selected Group
		while ( gidx<obj->stats.grps && i==obj->grps[gidx].line )     gidx++;
		// Current Test implementation is limited to 1 group ( SelectGroup[0] )
		if ( Negate ) {
			if ( gidx && strcmp(SelectGroup[0],obj->grps[gidx-1].name)==0 ) {
			    obj->faces[i].nodes = 0;
			    continue;
			}
		}
		else {
			if ( ! gidx || strcmp(SelectGroup[0],obj->grps[gidx-1].name)!=0 ) {
			    obj->faces[i].nodes = 0;
			    continue;
			}
		}
	    }
	    // Remove face if it includes inexistent vertexs
	    for( n=0 ; n<obj->faces[i].nodes ; n++ ) {
		int v = obj->faces[i].Node[n][0];
		if ( v<1 || v>obj->stats.verts  || obj->counts.verts[v-1]<0 ) {
		    obj->faces[i].nodes = 0;
		    if ( Verbose>=2 )
			fprintf(stderr,"face %d removed.v=%d,count=%d\n", i, v, obj->counts.verts[v-1] );
		    break;
		}
	    }
	} 
}


void SetIndexs(ObjFile *obj) {
	// Assign an order number to each used vertex/texture coord/norm
	// Unused vertex/texture coord/norm will get no order number.
	int i;
	int nv=0, nt=0, nn=0;
	for( i=0 ; i<obj->stats.verts ; i++ ) 
	    if ( obj->counts.verts[i]>0 )	obj->order.verts[i] = nv++;

	for( i=0 ; i<obj->stats.texts ; i++ ) 
	    if ( obj->counts.texts[i]>0 )	obj->order.texts[i] = nt++;

	for( i=0 ; i<obj->stats.norms ; i++ ) 
	    if ( obj->counts.norms[i]>0 )	obj->order.norms[i] = nn++;
}

void FilterVerts(ObjFile *obj, Vert min, Vert max) {
	int i, c;
	for( i=0 ; i<obj->stats.verts ; i++ ) {
	    for( c=0 ; c<3 ; c++ ) {
		if ( obj->verts[i][c]<min[c] ) {
		    if ( ! SolidCut ) {
			obj->counts.verts[i] = NULL_IDX;
			break;
		    }
		    else obj->verts[i][c] = min[c];
		}
		if ( obj->verts[i][c]>max[c] ) {
		    if ( ! SolidCut ) {
			obj->counts.verts[i] = NULL_IDX;
			break;
		    }
		    else obj->verts[i][c] = max[c];
		}
	    }
	}
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
	}
}


void Usage() {
    fprintf(stderr,"Usage:\n\tObjTool [options] inputfile1.obj inputfile2.obj ...\n");
    fprintf(stderr,"\t\t-i                 file info\n");
    fprintf(stderr,"\t\t-xmin value        define xmin\n");
    fprintf(stderr,"\t\t-xmax value        define xmax\n");
    fprintf(stderr,"\t\t-ymin value        define ymin\n");
    fprintf(stderr,"\t\t-ymax value        define ymax\n");
    fprintf(stderr,"\t\t-zmin value        define zmin\n");
    fprintf(stderr,"\t\t-zmax value        define zmax\n");
    fprintf(stderr,"\t\t-o name            select object\n");
    fprintf(stderr,"\t\t-g name 	   select obj group\n");
    fprintf(stderr,"\t\t-m name	           select material\n");
    fprintf(stderr,"\t\t-tx value	   translate x\n");
    fprintf(stderr,"\t\t-ty value	   translate y\n");
    fprintf(stderr,"\t\t-tz value	   translate z\n");
    fprintf(stderr,"\t\t-sx value	   scale x\n");
    fprintf(stderr,"\t\t-sy value	   scale y\n");
    fprintf(stderr,"\t\t-sz value	   scale z\n");
    fprintf(stderr,"\t\t-s value	   global scale\n"); 
    fprintf(stderr,"\t\t-rx value	   rotate around xx axis\n");
    fprintf(stderr,"\t\t-ry value	   rotate around yy axis\n");
    fprintf(stderr,"\t\t-rz value	   rotate around zz axis\n");


    fprintf(stderr,"\t\t-R                 use relative coords\n");
    fprintf(stderr,"\t\t-M		   Do not output mtllib directives\n");
    fprintf(stderr,"\t\t-c                 solid cut\n");
    fprintf(stderr,"\t\t-n		   negate face filter condition( -g,-m )\n");
    fprintf(stderr,"\t\t-O outfile         output to outfile (default: stdout)\n");
    fprintf(stderr,"\t\t-e		   Explode outfile into single objects\n");
    fprintf(stderr,"\t\t-S shadow_file     shadow output to shadow_file (default: no shadow ouput)\n");
    
}

void InvalidOption(char *opt) {
    fprintf(stderr,"Invalid option: '%s'\n", opt);
    Usage();
}



int GetOptions(int argc, char** argv) {
    int i = 1;
    while( i<argc ) {
	if ( argv[i][0]=='-' ) {
	    switch (argv[i][1]) {
		case 'i':
		    info = 1;
		    break;
		case 'x':
		    switch (argv[i][3]) {
			case 'a':	// xmax
			    i++;
			    VMax[0] = atof(argv[i]);
			    break;
			case 'i':	// xmin
			    i++;
			    VMin[0] = atof(argv[i]);
			    break;
			default:
			    InvalidOption(argv[i]);
		    }
		    break;
		case 'y':
		    switch (argv[i][3]) {
			case 'a':	// ymax
			    i++;
			    VMax[1] = atof(argv[i]);
			    break;
			case 'i':	// ymin
			    i++;
			    VMin[1] = atof(argv[i]);
			    break;
			default:
			    InvalidOption(argv[i]);
		    }
		    break;
		case 'z':
		    switch (argv[i][3]) {
			case 'a':	// zmax
			    i++;
			    VMax[2] = atof(argv[i]);
			    break;
			case 'i':	// zmin
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
			    sscanf(argv[i],"%f,%f,%f",Translate,Translate+1,Translate+2);
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
		case 'n':
		    Negate = 1;
		    break;
		case 'e':
		    Explode = 1;
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

int JoinObjFiles(int nObjs, ObjFile ObjSet[], ObjFile *obj) {

	ResetObjStats( &(obj->stats) );

	int i;
	for( i=0 ; i<nObjs ; i++ ) {
		obj->stats.verts += ObjSet[i].stats.verts;
		obj->stats.faces += ObjSet[i].stats.faces;
		obj->stats.norms += ObjSet[i].stats.norms;
		obj->stats.texts += ObjSet[i].stats.texts;
		obj->stats.paras += ObjSet[i].stats.paras;
		obj->stats.grps += ObjSet[i].stats.grps;
		obj->stats.mats += ObjSet[i].stats.mats;
		obj->stats.objs += ObjSet[i].stats.objs;
		obj->stats.libs += ObjSet[i].stats.libs;
		obj->stats.shds += ObjSet[i].stats.shds;
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
		counters.libs += ObjSet[i].stats.libs;
		counters.objs += ObjSet[i].stats.objs;
		counters.shds += ObjSet[i].stats.shds;
		counters.grps += ObjSet[i].stats.grps;
		counters.mats += ObjSet[i].stats.mats;
	}

//printf("merge done\n");
	return nObjs;
}


int main(int argc, char **argv) {

	fprintf(stderr,"# ObjTool Compile date: %s\n", CDATE);

	if ( argc<=1 ) {
		Usage();
		exit(0);
	}
	int fi = GetOptions(argc,argv);
	//printf("# fi = %d, argc=%d\n",fi,argc);

	if ( fi>=argc ) {
	    fprintf(stderr,"Al least one input file must be specified.\n");
	    Usage();
	}


	// Alloc and load obj files
	int NObjects = argc-fi;
	ObjFile *ObjSet;
	ObjSet = Malloc(NObjects, sizeof(ObjFile));

	int i = 0;
	while ( fi<argc) {
		fprintf(stderr,"Loading ObjFile '%s'\n", argv[fi]);
		LoadObjFile(argv[fi], ObjSet+i);
		SetUseCounters(ObjSet+i);
		i++;
		fi++;
	}

	// All files got
	// Join all geometry into one objfile	
	ObjFile obj;	
	JoinObjFiles(NObjects,ObjSet,&obj);
	
	SetUseCounters(&obj);
	if ( info ) {
		    PrintObjStats(&(obj.stats));
		    exit(0);
        }
	if ( Verbose ) {
		    if ( Verbose>3 )
			PrintFullObjStats(&obj);
		    else
			PrintObjStats(&(obj.stats));
	}

	FilterVerts( &obj, VMin, VMax );

	//SetIndexs(&obj);
	CleanFaces(&obj);

	SetUseCounters(&obj);
	ProcVerts( &obj );
	SetIndexs(&obj);
	if ( Explode )
		ExplodeOutputFile(OutputFile, &obj);
	else
		SaveObjFile(OutputFile, &obj);

	if ( ShadowOutputFile ) {
printf("shadow\n");
			ObjFile *Shadow = CreateShadowObj(&obj);
PrintObjStats(&(Shadow->stats));
			SaveObjFile(ShadowOutputFile, Shadow);
printf("saved\n");
			FreeObjFile(Shadow);
	}
	FreeObjFile(&obj);


	return 0;
}
		


