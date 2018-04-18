#MACHINEFLAGS = -n32 -mips3
MACHINEFLAGS = 
CFLAGS = $(MACHINEFLAGS)
LFLAGS = -lm $(MACHINEFLAGS)
OBJLIST = readPDBrecs.o geom3d.o utility.o select.o \
	 atomprops.o stdconntable.o abin.o

CC = gcc

.c.o:
	${CC} -c $*.c $(CFLAGS)

bndlst: bndlst.o $(OBJLIST)
	${CC} -o $@ bndlst.o $(OBJLIST) $(LFLAGS)

clean:
	@rm -f *.o *.ckp

install: bndlst
	mv bndlst /local/bin

# DO NOT DELETE THIS LINE -- make depend uses it
abin.o:	./abin.h	./geom3d.h	./utility.h	abin.c
atomprops.o:	./atomprops.h	atomprops.c
bndlst.o:	./abin.h	./atomprops.h	./bndlst.h	./geom3d.h	./parse.h \
		./readPDBrecs.h	./select.h	./stdconntable.h	./utility.h \
		bndlst.c
geom3d.o:	./geom3d.h	geom3d.c
readPDBrecs.o:	./geom3d.h	./readPDBrecs.h	./utility.h	readPDBrecs.c
select.o:	./abin.h	./atomprops.h	./geom3d.h	./parse.h	./select.h \
		./stdconntable.h	./utility.h	select.c
stdconntable.o:	./stdconntable.h	stdconntable.c
utility.o:	utility.c
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
