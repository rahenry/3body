NAME = 3imb
CC = g++
CFLAGS = -std=c++0x -O3 -MMD -MP -Wall -Wno-unused-variable -fopenmp -I/usr/local/gsl/1.15-gcc/include/
#LDFLAGS = -O3 -DNDEBUG -msse2 -fopenmp -L/usr/local/gsl/1.15-gcc/lib -L/usr/local/lapack/3.5.0-gcc-4.8.2/lib -L/usr/local/blas/1.0.248-gcc-4.8.2/lib/
#LDFLAGS = -O3 -DNDEBUG -msse2 -fopenmp -L/usr/local/gsl/1.15-gcc/lib
#
# works up til libblas
#LDFLAGS = -O3 -DNDEBUG -msse2 -fopenmp -L/usr/local/gsl/1.15-gcc/lib -L/usr/local/lapack/3.4.2/lib
LDFLAGS = -O3 -DNDEBUG -msse2 -fopenmp -L/usr/local/gsl/1.15-gcc/lib -L/usr/local/lapack/3.4.2/lib -L/usr/local/blas/1.0.248/lib

LIBS = -lgsl -llapack -lblas -lgslcblas -lgfortran

SRCDIR = src
OBJDIR = build
SRC := $(wildcard $(SRCDIR)/*.cpp) 
OBJ := $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o) 


all: $(SRC) $(NAME)

$(NAME): $(OBJ) 
	$(CC) $(LDFLAGS) $(OBJ) -o $@ $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJDIR)/* $(NAME)

$(OBJDIR):
	mkdir -p $(OBJDIR)


-include $(OBJ:.o=.d)
# ....
