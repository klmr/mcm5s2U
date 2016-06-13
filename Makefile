ROOT_DIR := $(shell dirname $(realpath $(lastword ${MAKEFILE_LIST})))

report/%.html: code/%.r
	Rscript -e "ezknitr::ezspin('$<', out_dir = '${@D}', fig_dir = 'figures/$*/', wd = '${ROOT_DIR}', doc = '^## ', move_intermediate_file = TRUE)"
