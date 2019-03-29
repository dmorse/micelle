# ------------------------------------------------------------------
# Users may need to change the following Makefile macros:
#  - SRC_DIR, if this directory is moved or copied to another location.
#  - BIN_DIR, to change the directory for executables

# Directory containing Simpatico library source files
BLD_DIR_REL=../../../..

#-----------------------------------------------------------------------
# Include master makefiles

include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

#-------------------------------------------------------------------
# Major targets

all: $(mcMd_modules_micelle_analyzers_OBJS)

clean:
	-rm -f $(mcMd_modules_micelle_analyzers_OBJS)
	-rm -f $(mcMd_modules_micelle_analyzers_OBJS:.o=.d)

-include $(mcMd_modules_micelle_analyzers_OBJS:.o=.d)

#-------------------------------------------------------------------
