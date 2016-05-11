

CFLAGS=-Wall -O2 -g -DCDATE="\"`date`\""


all: ObjTool

ObjTool: ObjTool.o
	cc -o $@ -lm $^
	git add ObjToll.c
	git commitg
	git push


