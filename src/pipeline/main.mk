.PHONY: MAKE_DIRS all


SCRIPTS = $(wildcard ../scripts/*.sh)
TARGET_SCRIPTS = $(patsubst ../scripts/%, ${TARGET_DIR}/%, $(SCRIPTS))
PLGD_PM = $(wildcard ./pipeline/Plgd/*.pm)
TARGET_PLGD_PM = $(patsubst ./pipeline/Plgd/%, ${TARGET_DIR}/Plgd/%, $(PLGD_PM))

all: ${TARGET_DIR}/necat.pl ${TARGET_DIR}/Necat.pm ${TARGET_DIR}/necat.sh \
               ${TARGET_PLGD_PM} \
	       ${TARGET_SCRIPTS}

${TARGET_DIR}/necat.pl: pipeline/necat.pl
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR}/ ; fi
	cp -pf pipeline/necat.pl ${TARGET_DIR}/necat.pl
	chmod +x ${TARGET_DIR}/necat.pl

	
${TARGET_DIR}/Necat.pm: pipeline/Necat.pm
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR}/ ; fi
	cp -pf $< $@

${TARGET_DIR}/necat.sh: pipeline/necat.sh
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR} ; fi
	cp -pf pipeline/necat.sh ${TARGET_DIR}/necat.sh
	chmod +x ${TARGET_DIR}/necat.sh

${TARGET_DIR}/renecat.sh: pipeline/renecat.sh
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR} ; fi
	cp -pf pipeline/renecat.sh ${TARGET_DIR}/renecat.sh
	chmod +x ${TARGET_DIR}/renecat.sh

$(TARGET_PLGD_PM):${TARGET_DIR}/Plgd/% : ./pipeline/Plgd/% 
	@if [ ! -e ${TARGET_DIR}/Plgd ] ; then mkdir -p ${TARGET_DIR}/Plgd ; fi
	cp -pf  $< $@

$(TARGET_SCRIPTS):${TARGET_DIR}/% : ../scripts/%
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR} ; fi
	cp -pf  $^ $@
	chmod +x $@
	

clean:
	rm -f ${TARGET_DIR}/necat.pl
	rm -f ${TARGET_DIR}/necat.sh
	rm -f ${TARGET_SCRIPTS}
	rm -f ${TARGET_PLGD_PM}
	rm -rf ${TARGET_DIR}/Plgd
