FC = gfortran
# NOTE: reduction(.or.:var) requires OpenMP 4.5+ (gfortran 6+)
FFLAGS = -Wall -Wextra -O2 -pedantic -fopenmp
THREADS ?= 1

OBJECTS = grid_module.o main.o
EXECUTABLE = critical_cell

# Use /usr/bin/time for more detailed output, or just 'time' for the shell default
TIME = /usr/bin/time -p

debug: FFLAGS += -g -fcheck=bounds
debug: $(EXECUTABLE)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	@echo "Linking executable..."
	@$(TIME) $(FC) $(FFLAGS) -o $@ $^

grid_module.o: grid_module.f90
	@echo "Compiling grid_module..."
	@$(TIME) $(FC) $(FFLAGS) -c grid_module.f90

main.o: main.f90 grid_module.o
	@echo "Compiling main..."
	@$(TIME) $(FC) $(FFLAGS) -c main.f90

clean:
	rm -f *.o *.mod $(EXECUTABLE)

run: $(EXECUTABLE)
	@echo "Running execution timer with $(THREADS) threads..."
	@OMP_NUM_THREADS=$(THREADS) time ./$(EXECUTABLE)

rebuild: clean all

pipeline: clean run visualise clean

.PHONY: all debug clean run rebuild visualise pipeline
