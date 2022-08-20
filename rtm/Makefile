EXEC	= main
CXX	= gcc
OBJECTS	= io.o alloc.o
LIBS	= -lm
MAIN	= main.c
OPTIONS	= -O3 -Wall -fopenmp

main:	clean backup alloc io  
	$(CXX) $(OPTIONS) -o $(EXEC) $(OBJECTS) $(MAIN) $(LIBS)

exec:	clean alloc io  
	$(CXX) $(OPTIONS) -o $(EXEC) $(OBJECTS) $(MAIN) $(LIBS)
	./exec_tests.sh

io:
	$(CXX) $(OPTIONS) -c io.c

alloc:
	$(CXX) $(OPTIONS) -c alloc.c

edit_header:
	$(CXX) -Wall -o edit_header edit_header.c

copy_to_dca:
	scp -P2222 rtm_migrated_otimizado.su username@FQDN:/home/username/public_html/rtm/
#scp -P2222 fdbackward.bin username@FQDN:/home/username/public_html/rtm/
#scp -P2222 fdforward.bin username@FQDN:/home/username/public_html/rtm/
#scp -P2222 ../fortran/out.su  username@FQDN:/home/username/public_html/rtm/

copy_from_dca:
	wget SERVER/~username/rtm/rtm_migrated_otimizado.su
	wget SERVER/~username/rtm/fdbackward.bin
	wget SERVER/~username/rtm/fdforward.bin

clean:
	clear
	rm -f $(OBJECTS)
	rm -f *~
	rm -f fdbackward.bin
	rm -f fdforward.bin
	rm -f cshot.su
	rm -f rtm_migrated.su

backup:
	cp alloc.c  ../backup/
	cp alloc.h  ../backup/
	cp io.c	    ../backup/
	cp io.h     ../backup/
	cp segy.h   ../backup/
	cp main.c   ../backup/
	cp Makefile ../backup/

gcov:	clean clean_gcov
	$(CXX) -c -fprofile-arcs -ftest-coverage alloc.c
	$(CXX) -c -fprofile-arcs -ftest-coverage io.c
	$(CXX) -fopenmp $(OBJECTS) $(LIBS) -fprofile-arcs -ftest-coverage $(MAIN)
	time ./a.out -p 24

clean_gcov:
	clear
	rm -f a.out
	rm -f alloc.gcno
	rm -f alloc.gcda
	rm -f alloc.c.gcov
	rm -f io.gcno
	rm -f io.gcda
	rm -f main.gcno
	rm -f main.gcda
	rm -f main.c.gcov
	rm -f gmon.out

gprof: clean clean_gprof
	$(CXX) -c -pg alloc.c
	$(CXX) -c -pg io.c
	$(CXX) -fopenmp -pg $(OBJECTS) $(LIBS) $(MAIN)
	time ./a.out -p 24

clean_gprof:
	clear
	rm -f a.out
	rm -f gmon.out

