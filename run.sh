#!/bin/bash

#./main ~/workspace/my_pso_freq/data/displacement.csv "16" 5 50; PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py found_spectrum.log "pdfs/channel_16.pdf"
#num_modes=12
num_modes=5
num_iter=20
num_particles=100

c=$1

data_file="data/04_mod_A20_FRF_displacement.csv"

pdfs_dir="pdfs_${num_modes}_modes_04_mod_A20^2"
logs_dir="logs_${num_modes}_modes_04_mod_A20^2"
pdf_file="channel_${c}.pdf"
log_file="04_mod_A20_FRF_log_channel_${c}.log"

python_path="PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages"

echo "$pdfs_dir"
echo "$logs_dir"
mkdir -p "$pdfs_dir"
mkdir -p "$logs_dir"

rm -f "$logs_dir/04_mod_A20_FRF_log_channel_$c.log"
rm -f "$logs_dir/04_mod_A20_FRF_log_channel_$c.log_main"

#echo "$python_path" ./plot.py "$logs_dir/$log_file" "$pdfs_dir/pdf_file"
#gdb -ex run --args ./main "$data_file" "$logs_dir//$log_file" "$c" 1001 "$num_modes" "$num_iter" "$num_particles"
./main "$data_file" "$logs_dir/$log_file" "$c" 1001 "$num_modes" "$num_iter" "$num_particles"

PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py "$logs_dir/${log_file}_main" "$pdfs_dir/$pdf_file"
#PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py "$logs_dir/$log_file" "$pdfs_dir/pdf_file"

#./main_m "$data_file" "$logs_dir/$log_file" "$c" 1001 "$num_modes" "$num_iter" "$num_particles"

#for c in {2..16}
#do
#    echo "$c"
#    ./main "$data_file $logs_dir/$log_file" "$c" 1001 "$num_modes" "$num_iter" "$num_particles" 50; PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py "$logs_dir/$log_file" "$pdfs_dir/pdf_file"
#done
echo "$python_path ./plot.py $logs_dir/${log_file}_main $pdfs_dir/$pdf_file"
