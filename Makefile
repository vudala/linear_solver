CC = gcc -g -std=c11
LFLAGS = -lm
OUTPUT = labSisLin 
OBJS = utils.o SistemasLineares.o

.PHONY: clean purge all

$(OUTPUT) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

clean:
	@rm -f *~ *.bak

purge: clean
	@rm -f *.o core a.out
	@rm -f $(OUTPUT)

run: $(OUTPUT) 
	./$(OUTPUT) < input.in
