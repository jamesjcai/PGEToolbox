SOURCES =  hweStrata.c get_optval.c gammln.c gammp.c gcf.c gser.c pchisq.c

OBJECTS = $(SOURCES:.c=.obj)
PROGRAM = hweStrata

CC = cl 

$(PROGRAM):  $(OBJECTS)  
	 $(CC) -O -o $(PROGRAM) $(OBJECTS) -lm
clean:
	rm -f $(OBJECTS)

