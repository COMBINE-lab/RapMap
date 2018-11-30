HEADERS=(EnumCoder.hpp HitManager.hpp IndexHeader.hpp JFRaw.hpp PairAlignmentFormatter.hpp RapMapConfig.hpp RapMapFileSystem.hpp RapMapIndex.hpp RapMapSAIndex.hpp RapMapUtils.hpp SACollector.hpp SASearcher.hpp ScopedTimer.hpp SingleAlignmentFormatter.hpp SparseHashSerializer.hpp SpinLock.hpp)

SRCS=(EnumCoder.cpp EnumCoderTest.cpp HitManager.cpp RapMap.cpp RapMapFileSystem.cpp RapMapIndex.cpp RapMapIndexer.cpp RapMapMapper.cpp RapMapSAIndex.cpp RapMapSAIndexer.cpp RapMapSAMapper.cpp RapMapUtils.cpp UtilTest.cpp) 

for h in ${HEADERS[*]};
do
copyright-header --add-path ../include/$h \
                 --license GPL3 \
                 --copyright-holder 'Rob Patro' \
                 --copyright-holder 'Avi Srivastava' \
                 --copyright-holder 'Hirak Sarkar' \
                 --copyright-software 'RapMap' \
                 --copyright-software-description "Rapid and accurate mapping of short reads to transcriptomes using quasi-mapping." \
                 --copyright-year 2015 \
                 --copyright-year 2016 \
		         --word-wrap 80 \
                 --output-dir .
done

for h in ${SRCS[*]};
do
copyright-header --add-path ../src/$h \
                 --license GPL3 \
                 --copyright-holder 'Rob Patro' \
                 --copyright-holder 'Avi Srivastava' \
                 --copyright-holder 'Hirak Sarkar' \
                 --copyright-software 'RapMap' \
                 --copyright-software-description "Rapid and accurate mapping of short reads to transcriptomes using quasi-mapping." \
                 --copyright-year 2015 \
                 --copyright-year 2016 \
		         --word-wrap 80 \
                 --output-dir .
done
