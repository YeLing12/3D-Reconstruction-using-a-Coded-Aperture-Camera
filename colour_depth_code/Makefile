
ifeq (, $(findstring CYGWIN, $(shell uname -s)))
  EXE =
  CORE = core*
else
  EXE = .exe
  CORE = *.stackdump
endif


MRF       = ./MRFx.x


CC        = g++
CFLAGS    = -Wall -O2 -I. -I$(MRF)
LDFLAGS   = -L$(MRF)
LDLIBS    = -lm -lmrf


TARGET    = main
BIN      := $(addsuffix $(EXE), $(TARGET))

SRCS     := $(wildcard *.c)
OBJS     := $(patsubst %.c, %.o, $(SRCS))


.PHONY: all clean


all: $(TARGET)

$(TARGET): $(OBJS)

clean:
	$(RM) $(BIN) *.o $(CORE)
