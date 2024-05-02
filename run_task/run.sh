o2-analysis-lf-epvector -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-pid-tof-base -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-ft0-corrected-table -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-timestamp -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-event-selection -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-multiplicity-table -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-centrality-table -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-lf-nuclei-spectra -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 --aod-writer-keep AOD/NUCLEITABLE/0,AOD/NUCLEITABLEFLOW/0 | \
o2-analysis-lf-flow-qc -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-track-propagation -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-trackselection -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 | \
o2-analysis-qvector-table -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000
