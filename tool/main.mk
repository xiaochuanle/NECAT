.PHONY: all pigz

all: pigz

pigz: ${TARGET_DIR}/pigz

${TARGET_DIR}/pigz: ${TARGET_DIR}/../pigz-2.4/pigz
	cp -pf $^ $@

${TARGET_DIR}/../pigz-2.4/pigz: ../tool/pigz-2.4.tar.gz
	tar -xzf ../tool/pigz-2.4.tar.gz -C ${TARGET_DIR}/.. || exit 255
	cd ${TARGET_DIR}/../pigz-2.4 && make CC=${CC} LDFLAGS="${LDFLAGS}"

clean:
	rm -f ${TARGET_DIR}/pigz
	rm -rf ${TARGET_DIR}/../pigz-2.4
