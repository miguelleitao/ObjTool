

CFLAGS=-Wall -O2 -g -DCDATE="\"`date`\""


all: ObjTool

ObjTool: ObjTool.o
	cc -o $@ -lm $^
	git add ObjTool.c
	git commitg
	git push


