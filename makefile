CFLAGS = -lm -lpthread -g -Wall -W -Wextra -O3 --fast-math 
prog.exe: main.o qr_razlog.o synchronize.o get_time.o
	g++ $(CFLAGS) main.o qr_razlog.o synchronize.o get_time.o -o prog.exe -lrt -pg
main.o: main.cpp qr_razlog.h get_time.h
	g++ $(CFLAGS) -c main.cpp     
qr_razlog.o: qr_razlog.cpp qr_razlog.h get_time.h
	g++ $(CFLAGS) -c qr_razlog.cpp  
synchronize.o: synchronize.cpp synchronize.h
	g++ $(CFLAGS) -c synchronize.cpp 
get_time.o: get_time.cpp get_time.h
	g++ $(CFLAGS) -c get_time.cpp 
clean:
	rm -rf *.o prog.exe