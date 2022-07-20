TARGET=SLRS
CC=g++ --std=c++11 -no-pie
CPPFLAGS = -g -Wall -O3

SRC=./src
INC=-I ./include
OBJS=$(addsuffix .o, $(basename $(wildcard $(SRC)/*.cpp)))

$(TARGET):$(OBJS)
        $(CC) -o $@ $^ ./lib/lp/liblpsolve55.a -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lm -ldl -lz

%.o: %.cpp
        $(CC) $(CPPFLAGS) -c $< -o $@ $(INC)-I $(BAMTOOLS_HOME_INCLUDE)/ -lz 
    
all: SLRS
        rm -f *.o

.PHONY: clean

clean:
        rm -f *.o
        rm SLRS
