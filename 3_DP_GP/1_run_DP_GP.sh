# Run the DP_GP analysis 
# Severin Einspanier 
# Perform DP_GP
DP_GP_cluster.py -i in_data/LA1282_LFC_LRT.txt \
    -o out_data/1282_LTR/2025_02_02_plowered_infmockLRT_MAP \
    -p png pdf -n 2000 --plot -a 1 -c MAP

# extract the plots only 

DP_GP_cluster.py -i in_data/LA1282_LFC_LRT.txt \
    -o out_data/1282_LTR/2025_02_03_EDIT_FOR_PDF \
    -p png pdf -n 2000 --plot -a 1 -c MAP

# get PDF for 1809

DP_GP_cluster.py -i in_data/LA1809_LFC_LRT.txt \
    -o out_data/1809_LTR/2025_02_03_EDIT_FOR_PDF \
    -p png pdf -n 2000 --plot -a 1 -c MAP
# new degs

DP_GP_cluster.py -i in_data/LA1282_48hpi_DEGs.txt \
    -o out_data/1282_48degs/2025_02_15_48erDEGs \
    -p png pdf -n 2000 --plot -a 1 -c MAP
