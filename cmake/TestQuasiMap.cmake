set(QUASI_INDEX_CMD ${CMAKE_BINARY_DIR}/rapmap quasiindex -t transcripts.fasta -i sample_quasi_index)
execute_process(COMMAND ${QUASI_INDEX_CMD}
                WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
                RESULT_VARIABLE QUASI_INDEX_RESULT
                )

if (QUASI_INDEX_RESULT)
    message(FATAL_ERROR "Error running ${QUASI_INDEX_COMMAND}")
endif()

set(MAP_COMMAND ${CMAKE_BINARY_DIR}/rapmap quasimap -t 2 -i sample_quasi_index -1 reads_1.fastq -2 reads_2.fastq -o sample_quasi_map.sam)
execute_process(COMMAND ${MAP_COMMAND}
	            WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
                RESULT_VARIABLE QUASI_MAP_RESULT
                )
if (QUASI_MAP_RESULT)
    message(FATAL_ERROR "Error running ${QUASI_MAP_RESULT}")
endif()

if (EXISTS ${TOPLEVEL_DIR}/sample_data/sample_quasi_map.sam)
    message("RapMap (quasi) ran successfully")
else()
    message(FATAL_ERROR "RapMap (quasi-index & map) failed to produce output")
endif()
