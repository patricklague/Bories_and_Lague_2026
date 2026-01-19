#!/usr/bin/env bash
# average_profiles.sh
# Usage: ./average_profiles.sh dens_mono_frame_*.dat > dens_mono_avg.dat

files=("$@")
nfiles=${#files[@]}

awk -v n=$nfiles '
{
    z = $1
    val = $2
    sum[z] += val
    count[z]++
}
END {
    for (z in sum) {
        avg = sum[z] / count[z]
        printf "%7.2f\t%.5f\n", z, avg
    }
}' "${files[@]}" | sort -n


