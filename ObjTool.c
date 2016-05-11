
/*
 *	ObjTool.c
 *
 *	Command line toolset to work with Wavefront 3D Obj files.
 *
 *	Copyright (C) 2014, Paulo Neto, pneto@magnocomp.com
 *
 */

#define VERSION 0.12
 	
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_LINE_LEN (5880)
#define MAX_FACE_NODES (380)

#define MAX_FLOAT (1.0e32)

#define NULL_IDX (-999999)

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
char *Material = NULL;
char *OutputFile = NULL;
char *ShadowOutputFile = NULL;
short int SolidCut = 0;
short int Relative = 0;
short int Verbose = 0;
short int info = 0;
char *SelectGroup[20];
int Groups = 0;



void *Malloc(int n, size_t dim) {
	long size = n*dim;
	void *buf = malloc(size+8);
	if ( ! buf ) {
		fprintf(stderr,"Cannot alloc memory\n");
		exit(3);
	}
	*((long *)buf) = MEM_TAG;
	memset(buf+8,0,size);
	return buf+8;
}

void Free(void *p) {
	long *tag_p = p-8;
	if ( *tag_p == MEM_TAG ) {
		*tag_p = 0L;
		free(p-8);
//printf("Memory freed\n");
	}
//else printf("##Memory not freed\n");
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
printf("Creating shadow\n");
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
void LoadObjFile(char *fname, ObjFile *obj) {
	int nl = 0;
	int i;
	int res = GetObjStats( fname, &(obj->stats) );
	if ( ! res ) return;
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
			case 'm':
				obj->libs[counters.libs].name = strdup(lptr+7);
				obj->libs[counters.libs].line = counters.faces;
				counters.libs += 1;
				continue;
			case 'o':
				obj->objs[counters.objs].name = strdup(lptr+2);
				obj->objs[counters.objs].line = counters.faces;
				counters.objs += 1;
				continue;
			case 's':
				obj->shds[counters.shds].name = strdup(lptr+2);
				obj->shds[counters.shds].line = counters.faces;
				counters.shds += 1;
				continue;
			case 'g':
				obj->grps[counters.grps].name = strdup(lptr+2);
				obj->grps[counters.grps].line = counters.faces;
				counters.grps += 1;
				continue;
			case 'u':	// usemtl
				obj->mats[counters.mats].name = strdup(lptr+7);
				obj->mats[counters.mats].line = counters.faces;
				counters.mats += 1;
				continue;
			case 'v':
				switch ( *(lptr+1) ) {
					case ' ':
						ReadVec3f(lptr+2,obj->verts[counters.verts]);
						for( i=0 ; i<3 ; i++ ) {
						    if ( obj->verts[counters.verts][i]<obj->stats.vmin[i] )
							obj->stats.vmin[i] = obj->verts[counters.verts][i];
						    if ( obj->verts[counters.verts][i]>obj->stats.vmax[i] )
							obj->stats.vmax[i] = obj->verts[counters.verts][i];
						}
						counters.verts += 1;
						break;
					case 'n':
						ReadVec3f(lptr+2,obj->norms[counters.norms]);
						counters.norms += 1;
						break;
					case 't':
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
printf("vai abrir ficheiro %s\n",fname);
	if ( fname && *fname ) {
		fout = fopen(fname,"w");
		if ( ! fout ) {
		    fprintf(stderr,"Cannot open output file '%s'\n", fname);
		    exit(5);
		}
	}
	printf("    abriu\n");
	for( i=0 ; i<obj->stats.verts ; i++ ) 
	    if ( obj->counts.verts[i] > 0 ) {
		fprintf(fout,"v %f %f %f\n", obj->verts[i][0],obj->verts[i][1],obj->verts[i][2]);
		nv++;
	    }
	printf("    gravou verts\n");
	for( i=0 ; i<obj->stats.texts ; i++ ) 
	    if ( obj->counts.texts[i] > 0 ) {
		fprintf(fout,"vt %f %f %f\n", obj->texts[i][0],obj->texts[i][1],obj->texts[i][2]);
		nt++;
	    }

	printf("    gravou text verts\n");
	for( i=0 ; i<obj->stats.norms ; i++ ) 
	    if ( obj->counts.norms[i] > 0 ) {
		fprintf(fout,"vn %f %f %f\n", obj->norms[i][0],obj->norms[i][1],obj->norms[i][2]);
		nn++;
	    }
	printf("    gravou norms\n");
	int gidx=0,oidx=0,midx=0,lidx=0;
	for( i=0 ; i<obj->stats.faces ; i++ ) {
	    int n;

	    while ( lidx<obj->stats.libs && i==obj->libs[lidx].line ) {
		    fprintf(fout,"mtllib %s\n",obj->libs[lidx].name);
		    lidx++;
	    }
	    while ( oidx<obj->stats.objs && i==obj->objs[oidx].line ) {
		if ( oidx==obj->stats.objs-1 || obj->objs[oidx].line!=obj->objs[oidx+1].line )
		    fprintf(fout,"o %s\n",obj->objs[oidx].name);
		oidx++;
	    }
	    while ( gidx<obj->stats.grps && i==obj->grps[gidx].line ) {
		if ( gidx==obj->stats.grps-1 || obj->grps[gidx].line!=obj->grps[gidx+1].line )
		    fprintf(fout,"g %s\n",obj->grps[gidx].name);
		gidx++;
	    }
	    while ( midx<obj->stats.mats && i==obj->mats[midx].line ) {
		if ( midx==obj->stats.mats-1 || obj->mats[midx].line!=obj->mats[midx+1].line )
		    fprintf(fout,"usemtl %s\n",obj->mats[midx].name);
		midx++;
	    }


	    if ( obj->faces[i].nodes <= 0 ) continue;
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
printf("fechando\n");
	if ( fout!=stdout ) fclose(fout); 
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
	int midx = 0;
	for( i=0 ; i<obj->stats.faces ; i++ ) {
	    if ( Material!=NULL ) {
		while ( midx<obj->stats.mats && i==obj->mats[midx].line )     midx++;
		if ( ! midx || strncmp(Material,obj->mats[midx-1].name,strlen(Material))!=0 ) {
		    obj->faces[i].nodes = 0;
		    continue;
		}
	    }
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
	int i, c;
	for( i=0 ; i<obj->stats.verts ; i++ ) {
	    if ( obj->counts.verts[i] > 0 ) {
	    	for( c=0 ; c<3 ; c++ ) {
		    obj->verts[i][c] *= Scale[c];
		    obj->verts[i][c] += Translate[c];
		}
	    }
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
    fprintf(stderr,"\t\t-g name 	   select obj group\n");
    fprintf(stderr,"\t\t-m name	           select material\n");
    fprintf(stderr,"\t\t-tx value	   translate x\n");
    fprintf(stderr,"\t\t-ty value	   translate y\n");
    fprintf(stderr,"\t\t-tz value	   translate z\n");
    fprintf(stderr,"\t\t-sx value	   scale x\n");
    fprintf(stderr,"\t\t-sy value	   scale y\n");
    fprintf(stderr,"\t\t-sz value	   scale z\n");
    fprintf(stderr,"\t\t-s value	   global scale\n"); 

    fprintf(stderr,"\t\t-r                 use relative coords\n");
    fprintf(stderr,"\t\t-c                 solid cut\n");
    fprintf(stderr,"\t\t-o outfile         output to outfile (default: stdout)\n");
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
		case 'm':
		    i++;
		    Material = argv[i];
		    break;
/*
		case 'g':
		    i++;
		    Group = argv[i];
		    break;
*/
		case 'o':
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
		case 'r':
		    Relative = 1;
		    break;
		case 'v':
		    Verbose += 1;
		    break;
		case 'g':
		    i++;
		    SelectGroup[Groups++] = argv[i];
		default:
		    InvalidOption(argv[i]);
	    }
	}
	i++;
    }
    return i;
}


int main(int argc, char **argv) {

printf("Compile date: %s\n", CDATE);

	ObjFile obj;

	if ( argc<=1 ) {
		Usage();
		exit(0);
	}
	int fi = GetOptions(argc,argv)-1;
printf("fi = %d, argc=%d\n",fi,argc);
	if ( fi<argc) {
		LoadObjFile(argv[fi], &obj);
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
		SaveObjFile(OutputFile, &obj);
printf("shadow\n");
		if ( ShadowOutputFile ) {

			ObjFile *Shadow = CreateShadowObj(&obj);
PrintObjStats(&(Shadow->stats));
			SaveObjFile(ShadowOutputFile, Shadow);
printf("saved\n");
			FreeObjFile(Shadow);
		}
		FreeObjFile(&obj);
	}
	else {
	    fprintf(stderr,"Al least one input file must be specified.\n");
	    Usage();
	}
	return 0;
}
		


