#
# adapted from: https://fortran-lang.org/learn/building_programs/project_make
#
# Disable the default rules
.SUFFIXES:
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
SHELL=/bin/bash

# Project name
NAME := cales

TARGET := $(NAME)
INPUT_FILE := input.nml

ROOT_DIR := .
SRCS_DIR := $(ROOT_DIR)/src
BUILD_DIR := $(ROOT_DIR)/build
EXE_DIR := $(BUILD_DIR)
CONFIG_DIR := $(ROOT_DIR)/configs
LIBS_DIR := $(ROOT_DIR)/dependencies
LIBS :=
INCS :=

DEFINES :=

EXE := $(EXE_DIR)/$(TARGET)

# Configuration settings
FC := mpifort
FFLAGS :=
AR := ar rcs
LD := $(FC)
RM := rm -f
GD := $(ROOT_DIR)/src/.gen-deps.awk
CPP := -cpp

# Edit build.conf file desired
-include $(ROOT_DIR)/build.conf
-include $(CONFIG_DIR)/compilers.mk
-include $(CONFIG_DIR)/flags.mk
-include $(CONFIG_DIR)/libs.mk

# List of all source files
SRCS_INC := $(wildcard $(SRCS_DIR)/*.h90)
SRCS := $(filter-out $(SRCS_INC), $(wildcard $(SRCS_DIR)/*.f90))

# Add source directory to search paths
vpath % .:$(SRCS_DIR)
vpath % $(patsubst -I%,%,$(filter -I%,$(INCS)))

# Define a map from each file name to its object file
obj = $(BUILD_DIR)/$(notdir $(src)).o
$(foreach src, $(SRCS), $(eval $(src) := $(obj)))

# Create lists of the build artefacts in this project
OBJS := $(patsubst $(SRCS_DIR)/%, $(BUILD_DIR)/%.o, $(SRCS))
DEPS := $(BUILD_DIR)/.depend.mk

# Declare all public targets
.PHONY: all clean allclean libs libsclean
all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) $^ $(LIBS) $(INCS) -o $(EXE)

# Create object files from Fortran source
$(OBJS): $(BUILD_DIR)/%.o: $(SRCS_DIR)/%
	$(FC) $(FFLAGS) $(CPP) $(DEFINES) $(INCS) $(FFLAGS_MOD_DIR) $(BUILD_DIR) -c -o $@ $<

# Process the Fortran source for module dependencies
# Create the build directory automatically if missing
$(DEPS):
	@mkdir -p $(BUILD_DIR) 
	@echo '# This file contains the module dependencies' > $(DEPS)
	@$(foreach file, $(SRCS), $(GD) $(file) >> $(DEPS))

# Define all module interdependencies
-include $(DEPS)
$(foreach dep, $(OBJS), $(eval $(dep): $($(dep))))

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(BUILD_DIR)/*.{i,mod,smod,d,o} $(EXE) $(DEPS)

allclean:
	@make libsclean
	@make clean
#
# rules for building the external libraries (compile with 'make libs'):
#
include $(LIBS_DIR)/external.mk

# Rules to generate config files if missing
$(CONFIG_DIR)/compilers.mk:
	@echo "Generating $(CONFIG_DIR)/compilers.mk from $(CONFIG_DIR)/compilers.mk.example..."
	cp $(CONFIG_DIR)/compilers.mk.example $(CONFIG_DIR)/compilers.mk
$(CONFIG_DIR)/flags.mk:
	@echo "Generating $(CONFIG_DIR)/flags.mk from $(CONFIG_DIR)/flags.mk.example..."
	cp $(CONFIG_DIR)/flags.mk.example $(CONFIG_DIR)/flags.mk
$(CONFIG_DIR)/libs.mk:
	@echo "Generating $(CONFIG_DIR)/libs.mk from $(CONFIG_DIR)/libs.mk.example..."
	cp $(CONFIG_DIR)/libs.mk.example $(CONFIG_DIR)/libs.mk
$(ROOT_DIR)/build.conf:
	@echo "Generating $(ROOT_DIR)/build.conf from $(ROOT_DIR)/build.conf.example..."
	cp $(ROOT_DIR)/build.conf.example $(ROOT_DIR)/build.conf
