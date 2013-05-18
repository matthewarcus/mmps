# $Revision: 1.7 $
# MMPS main makefile

EXES := project addgrid stars testimage circularize combine
OBJS := matrix.o transform.o equations.o utils.o image.o

all: $(EXES)

ppms:
	cd images; mogrify -format ppm `ls *.jpg *.gif *.png 2> /dev/null || true`

test: project images/earth.ppm
	./project perspective -lat 20 -long 30 -x 6 -sun -date 172 -time 3.83 -f images/earth.ppm | display

test2: project images/earth.ppm
	./project gnomonic -lat 52.2 -dlat 52.2 -grid \
	-color darkyellow -tropics -temporaryhours \
        -color darkred -analemma \
	-color orange -altitudes \
        -color red -dateline 91.583 -datetime 91.583 \
	-f images/earth.ppm | display

test3: project images/earth.ppm
	./project mercator -tilt 90 -grid -gridx 20 -gridy 20 -color 255:0:0 -grid -scale 2 -f images/earth.ppm | display -rotate 270 -

startest: stars
	cat starlist | ./stars -ra 5.24 -dec -8.2 | display

startest2: stars
	cat starlist | ./stars hammer -grid -gridx 20 -gridy 20 -w 1024 -adjust | display

images/%.ppm: images/%.jpg
	convert -format ppm $< $@

release:
	$(MAKE) -f Makefile.release $@

.PHONY: startest startest2 test ppms release


# Generic after here
CC 	 := g++

OPTIMIZE := -O3
WARN     := -Wall -Wshadow

ifeq (EXTRAWARNINGS,true)
WARN += -W -Wold-style-cast
WARN += -Woverloaded-virtual -Wsign-promo -Wundef -Wcast-align
WARN += -Wwrite-strings -Wconversion
#WARN += -Weffc++
endif

DEFS 	 :=

$(EXES): %: %.o $(OBJS)

CFLAGS = -g -MMD -MP -ansi $(WARN) $(OPTIMIZE) $(DEFS)
LDFLAGS = -lm

ifeq ($(PROFILE),true)
CFLAGS += -pg
LDFLAGS += -pg
endif

# Look for files in parent directory so we can compile in a subdirectory
vpath %.cpp ..
vpath %.h ..

# Include auto-generated dependency files
-include *.d

# Rule for all of our executables
$(EXES):
	$(CC) $(LDFLAGS) -o $@ $^

# And a rule for .cpp files, just in case we don't have one.
%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

%.a:
	$(AR) -crs $@ $^

clean:
	rm -f $(EXES)
	rm -f *.o *.d *.a

purge:	clean
	rm -f *~ a.out gmon.out foo.* *.tgz image0*.ppm

.PHONY: all clean purge
