EXEC	= main 
CXX	= gcc
OBJECTS	= io.o memcpy.o alloc.o
LIBS	= -lm
MAIN	= main.c
OPTIONS	= -O3 -Wall

main:	clean alloc memcpy_d io  
	$(CXX) $(OPTIONS) -o $(EXEC) $(OBJECTS) $(MAIN) $(LIBS)

io:
	$(CXX) $(OPTIONS) -c io.c

alloc:
	$(CXX) $(OPTIONS) -c alloc.c

memcpy_d:
	$(CXX) $(OPTIONS) -c memcpy.c

edit_header:
	$(CXX) -Wall -o edit_header edit_header.c

copy_dca:
	scp -P2222 rtm_migrated_otimizado.su username@FQDN:/home/username/public_html/rtm/
#scp -P2222 ../fortran/out.su  username@FQDN:/home/username/public_html/rtm/

clean:
	clear
	rm -f $(OBJECTS)
	rm -f *~
	rm -f fdbackward.bin
	rm -f fdforward.bin
	rm -f cshot.su
	rm -f rtm_migrated.su
