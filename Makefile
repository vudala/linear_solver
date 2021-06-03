CC = gcc
LFLAGS = -lm
OUTPUT = labSisLin 
OBJS = utils.o SistemasLineares.o

.PHONY: clean purge all $(OUTPUT)

$(OUTPUT) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS)$<

clean:
	@rm -f *~
	@rm -f *.o

purge: clean
	@rm -f $(OUTPUT)

run: $(OUTPUT) 
	./$(OUTPUT) < sistemas.dat
