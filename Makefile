

CFLAGS=-Wall -O2 -I linmath -DCDATE="\"`date`\""

all: ObjTool

ObjTool: ObjTool.o intersection.o load_obj.o
	$(CC) -o $@ $^ -lm

ObjTool.o: ObjTool.c linmath/linmath.h
	$(CC) -c ${CFLAGS} $<

linmath:
	git submodule init
	git submodule update

install: ObjTool
	cp $^ /usr/local/bin

clean:
	$(RM) *.o ObjTool

push:
	git add ObjTool.c intersection.c load_obj.c ObjTool.h Makefile README.md
	git commit -m "update"
	git push


