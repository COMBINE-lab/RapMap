set(QUASI_INDEX_CMD_PH ${CMAKE_BINARY_DIR}/rapmap quasiindex --perfectHash -t transcripts.fasta -i sample_quasi_index_ph)
execute_process(COMMAND ${QUASI_INDEX_CMD_PH}
                WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
                RESULT_VARIABLE QUASI_INDEX_RESULT_PH
                )

if (QUASI_INDEX_RESULT_PH)
    message(FATAL_ERROR "Error running ${QUASI_INDEX_COMMAND_PH}")
endif()

set(MAP_COMMAND_PH ${CMAKE_BINARY_DIR}/rapmap quasimap -t 2 -i sample_quasi_index_ph -1 reads_1.fastq -2 reads_2.fastq -o sample_quasi_map_ph.sam)
execute_process(COMMAND ${MAP_COMMAND_PH}
	            WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
                RESULT_VARIABLE QUASI_MAP_RESULT_PH
                )
if (QUASI_MAP_RESULT_PH)
    message(FATAL_ERROR "Error running ${QUASI_MAP_RESULT_PH}")
endif()

if (EXISTS ${TOPLEVEL_DIR}/sample_data/sample_quasi_map_ph.sam)
    message("RapMap (quasi, perfect-hash) ran successfully")
else()
    message(FATAL_ERROR "RapMap (quasi-index & map) failed to produce output")
endif()
