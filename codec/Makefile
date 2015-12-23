CC=gcc
CFLAGS=`pkg-config libosmocore --cflags` -Wall -O2 -march=native
LDLIBS=`pkg-config libosmocore --libs` -lm

OBJS=ir77_ambe_decode

all: $(OBJS)

ir77_ambe_decode: ir77_ambe_decode.o ambe.o ecc.o frame.o math.o synth.o tables.o tone.o

clean:
	rm -f *.o $(OBJS)
