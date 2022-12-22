

CFLAGS=-Wall -O2 -I linmath -DCDATE="\"`date`\""

all: ObjTool

ObjTool: ObjTool.o intersection.o load_obj.o
	$(CC) -o $@ $^ -lm

ObjTool.o: ObjTool.h
intersection.o: ObjTool.h
load_obj.o: ObjTool.h

ObjTool.h: linmath/linmath.h
	touch $@

linmath:
	git submodule init
	git submodule update

install: ObjTool
	cp $^ /usr/local/bin

shadow.png:
	convert -size 256x256 -define gradient:vector=128,64,128,0 gradient:'graya(0,0.5)-none' $@

clean:
	$(RM) *.o ObjTool

push:
	git add .
	@#ObjTool.c intersection.c load_obj.c ObjTool.h Makefile README.md
	git commit -m "update"
	git push


