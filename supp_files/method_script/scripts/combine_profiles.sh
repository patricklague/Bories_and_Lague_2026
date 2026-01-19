#!/usr/bin/env bash
# combine_profiles.sh
# Usage: ./combine_profiles.sh dens_mono_avg.dat dens_nonmono_avg.dat > dens_all_avg.dat

mono=$1
nonmono=$2

awk '
NR==FNR {mono[$1]=$2; next}
{
    z = $1
    dens_mono = (z in mono ? mono[z] : 0.0)
    dens_nonmono = $2
    printf "%7.2f\t%.5f\t%.5f\n", z, dens_mono, dens_nonmono
}' "$mono" "$nonmono" | sort -n

