

CFLAGS=-Wall -O2 -I linmath -DCDATE="\"`date`\""

all: ObjTool

ObjTool: ObjTool.o
	$(CC) -o $@ $^ -lm

install: ObjTool
	cp $^ /usr/local/bin

clean:
	$(RM) *.o ObjTool

push:
	git add ObjTool.c Makefile README.md
	git commit -m "update"
	git push


