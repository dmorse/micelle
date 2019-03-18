include $(SRC_DIR)/mcMd/modules/micelle/mcMoves/sources.mk
include $(SRC_DIR)/mcMd/modules/micelle/analyzers/sources.mk

mcMd_modules_micelle_= \
    mcMd/modules/micelle/MicelleMcModule.cpp \
    $(mcMd_modules_micelle_mcMoves_) \
    $(mcMd_modules_micelle_analyzers_) \

mcMd_modules_micelle_SRCS=\
         $(addprefix $(SRC_DIR)/, $(mcMd_modules_micelle_))
mcMd_modules_micelle_OBJS=\
         $(addprefix $(BLD_DIR)/, $(mcMd_modules_micelle_:.cpp=.o))

