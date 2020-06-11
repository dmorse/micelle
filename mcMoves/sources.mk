include $(SRC_DIR)/mcMd/modules/micelle/mcMoves/semigrand/sources.mk
include $(SRC_DIR)/mcMd/modules/micelle/mcMoves/micelle/sources.mk

mcMd_modules_micelle_mcMoves_=\
  $(mcMd_modules_micelle_mcMoves_semigrand_) \
  $(mcMd_modules_micelle_mcMoves_micelle_) \
  mcMd/modules/micelle/mcMoves/MicelleMcMoveFactory.cpp

mcMd_modules_micelle_mcMoves_SRCS=\
  $(addprefix $(SRC_DIR)/, $(mcMd_modules_micelle_mcMoves_))
mcMd_modules_micelle_mcMoves_OBJS=\
  $(addprefix $(BLD_DIR)/, $(mcMd_modules_micelle_mcMoves_:.cpp=.o))

