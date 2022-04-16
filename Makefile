# -*- Mode: Makefile; -*-

CCPLUS=g++
MPICC = mpic++
CFLAGS= -O3 -Wall -DRK2 -DAUSM -DTIMING
BINS=./bin/euler
INCLUDE= -I./include
OUPUT_DIR:= build/
SOURCES= $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, %.o, $(SOURCES))
OBJ_WITH_DIR= $(addprefix $(OUPUT_DIR), $(OBJS))
$(shell if [ ! -e $(OUPUT_DIR) ];then mkdir -p $(OUPUT_DIR); fi)

all: $(BINS)

$(BINS):$(OBJS)
	$(CCPLUS) $(OBJ_WITH_DIR) -o $@

%.o:%.cpp
	$(CCPLUS) -c $^ $(INCLUDE) $(CFLAGS) -o $(OUPUT_DIR)$@

$(OUPUT_DIR):
	mkdir $(OUPUT_DIR)

clean:
	 rm -rf $(BINS) $(OUPUT_DIR)
