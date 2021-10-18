#
#	Code written by krishna pande, 180101038
#	BTP
#

CC=g++ -std=c++11
OBJECT=main

$(OBJECT): system_model.o main.o application_model.o
	$(CC) system_model.o application_model.o main.o -o $(OBJECT)

main.o : main.cpp system_model.h
	$(CC) -c main.cpp

system_model.o: system_model.cpp system_model.h application_model.h
	$(CC) -c system_model.cpp system_model.h application_model.h

application_model.o: application_model.cpp application_model.h
	$(CC) -c application_model.cpp application_model.h 

clean:
	@rm -f $(OBJECT)  *.o