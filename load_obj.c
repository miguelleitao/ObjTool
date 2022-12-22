/*
 *	load_obj.c
 *
 *	Part of the ObjTool command line toolset.
 */


#define _GNU_SOURCE
 	 	
#include "ObjTool.h"


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
    if ( ! p ) {
        if ( Verbose>3 ) fprintf(stderr, "    No need to free unallocated memory\n");
        return;
    }
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

static int ReadVec2f(char *str, float *data) {
	int res;
	res = sscanf(str,"%f %f", data+0, data+1);
	return res;
}

static int ReadVec3f(char *str, float *data) {
	int res;
	res = sscanf(str,"%f %f %f", data+0, data+1, data+2);
	return res;
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
		while( *lptr && *lptr!=' ' && *lptr!='\t' && *lptr!='/' ) lptr++;
		if ( *lptr!='/' ) break;
		lptr++;
		if ( !*lptr || *lptr==' ' || *lptr=='\t' || *lptr=='\n' || *lptr=='\r' ) break;
	}
	return lptr;
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

