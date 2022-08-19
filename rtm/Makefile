EXEC	= main 
CXX	= gcc
OBJECTS	= io.o
LIBS	= -lm
MAIN	= main.c

main: clean io
	$(CXX) -Wall -o $(EXEC) $(OBJECTS) $(MAIN) $(LIBS)
io:
	$(CXX) -c io.c

edit_header:
	$(CXX) -Wall -o edit_header edit_header.c

copy_dca:
	scp -P2222 rtm_migrated.su username@FQDN:/home/username/public_html/rtm/
#scp -P2222 ../fortran/out.su  username@FQDN:/home/username/public_html/rtm/

clean:
	clear
	rm -f $(OBJECTS)
	rm -f *~
	rm -f fdbackward.bin
	rm -f fdforward.bin
	rm -f cshot.su
	rm -f rtm_migrated.su
