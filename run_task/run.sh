o2-analysis-lf-nuclei-flow -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000


/data/lbariogl/flow/LHC23_PbPb_pass2_new
o2-aod-merger --input dir_7/input_data.txt --output AO2D_dir7.root
