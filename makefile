BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_modules_micelle_OBJS)

show:
	-echo $(mcMd_modules_micelle_)

clean:
	rm -f $(mcMd_modules_micelle_OBJS) $(mcMd_modules_micelle_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_modules_micelle_OBJS:.o=.d)

-include $(mcMd_modules_micelle_OBJS:.o=.d)

