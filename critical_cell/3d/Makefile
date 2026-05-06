FC = gfortran
FFLAGS = -Wall -Wextra -O2 -pedantic

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
	@echo "Running execution timer..."
	@time ./$(EXECUTABLE)

rebuild: clean all

pipeline: clean run visualise clean

.PHONY: all debug clean run rebuild visualise pipeline
