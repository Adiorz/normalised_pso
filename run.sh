#!/bin/bash

#./main ~/workspace/my_pso_freq/data/displacement.csv "16" 5 50; PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py found_spectrum.log "pdfs/channel_16.pdf"
#num_modes=12
num_modes=2
num_iter=20
pdfs="pdfs_${num_modes}_modes_04_mod_A20^2"
logs="logs_${num_modes}_modes_04_mod_A20^2"

echo "$pdfs"
echo "$logs"
mkdir -p "$pdfs"
mkdir -p "$logs"

i=2
echo "PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py $logs/04_mod_A20_FRF_log_channel_$i.log $pdfs/channel_$i.pdf"
#gdb -ex run --args ./main_m data/04_mod_A20_FRF_displacement.csv "$logs/04_mod_A20_FRF_log_channel_$i.log" "$i" 1001 "$num_modes" "$num_iter"
./main data/04_mod_A20_FRF_displacement.csv "$logs/04_mod_A20_FRF_log_channel_$i.log" "$i" 1001 "$num_modes" "$num_iter"
PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py "$logs/04_mod_A20_FRF_log_channel_$i.log" "$pdfs/channel_$i.pdf"
#./main_m data/04_mod_A20_FRF_displacement.csv "$logs/04_mod_A20_FRF_log_channel_$i.log" "$i" 1001 "$num_modes" "$num_iter"

#for i in {2..16}
#do
#    echo "$i"
#    ./main data/04_mod_A20_FRF_displacement.csv "$logs/04_mod_A20_FRF_log_channel_$i.log" "$i" 1001 "$num_modes" 50; PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py "$logs/04_mod_A20_FRF_log_channel_$i.log" "$pdfs/channel_$i.pdf"
#done
