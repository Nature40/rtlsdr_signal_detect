CC=gcc
CFLAGS=-g -O2 -ffast-math -mcpu=cortex-a7 -mfloat-abi=hard -mfpu=neon-vfpv4 -Wall -fPIC
LDFLAGS=-lfftw3f -lm -lc
LDFLAGS+=`mysql_config --cflags --libs`

rtlsdr_signal_detect: rtlsdr_signal_detect.c
	${CC} ${CFLAGS} ${LDFLAGS} $< /usr/local/lib/libliquid.a -o $@
