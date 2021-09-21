CC = gcc
CFLAGS = -Wall -Werror
TARGET1 = ABC
TARGET2 = ABC.so
OBJECTS = standard_ABC.o

all : $(TARGET1) $(TARGET2)

$(TARGET1) : $(TARGET2)
	$(CC) $(CFLAGS) -o $@ main.c $^ -lm -Wl,-rpath=$(shell pwd)

$(TARGET2) : $(OBJECTS)
	$(CC) $(CFLAGS) -fPIC -shared -o $@ $^

clean :
	rm -f *.o $(TARGET1) $(TARGET2)
