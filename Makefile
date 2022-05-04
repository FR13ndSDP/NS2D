# -*- Mode: Makefile; -*-

CCPLUS=g++
MPICC = mpic++
CFLAGS= -O3 -Wall -DRK2
BINS=./bin/solver
INCLUDE= -I./include
OUPUT_DIR:= build
SOURCES= $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, %.o, $(SOURCES))
OBJ_WITH_DIR= $(addprefix $(OUPUT_DIR)/, $(OBJS))
$(shell if [ ! -e $(OUPUT_DIR) ];then mkdir -p $(OUPUT_DIR); fi)

all: $(BINS)

$(BINS):$(OBJ_WITH_DIR)
	$(CCPLUS) $(OBJ_WITH_DIR) -o $@

$(OUPUT_DIR)/%.o:%.cpp
	$(CCPLUS) -c $^ $(INCLUDE) $(CFLAGS) -o $@

clean:
	 rm -rf $(BINS) $(OUPUT_DIR)
