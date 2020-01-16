

CFLAGS=-Wall -O2 -I linmath -DCDATE="\"`date`\""

all: ObjTool

ObjTool: ObjTool.o intersection.o
	$(CC) -o $@ $^ -lm

linmath:
	git submodule init
	git submodule update

install: ObjTool
	cp $^ /usr/local/bin

clean:
	$(RM) *.o ObjTool

push:
	git add ObjTool.c Makefile README.md
	git commit -m "update"
	git push


