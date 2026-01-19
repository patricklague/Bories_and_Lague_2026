#!/usr/bin/env python3
"""
combine_profiles.py
Combine deux profils de densité (mono et nonmono) avec un merge "outer".
Toutes les valeurs de z sont conservées, les valeurs manquantes sont remplacées par 0.

Usage: python combine_profiles.py dens_mono_avg.dat dens_nonmono_avg.dat > dens_all_avg.dat
"""

import sys


def read_profile(filename):
    """Lit un fichier de profil et retourne un dictionnaire {z: density}."""
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                z = float(parts[0])
                density = float(parts[1])
                data[z] = density
    return data


def combine_profiles(mono_file, nonmono_file):
    """Combine deux profils avec un merge outer."""
    mono_data = read_profile(mono_file)
    nonmono_data = read_profile(nonmono_file)
    
    # Récupérer toutes les valeurs de z (union des deux ensembles)
    all_z = sorted(set(mono_data.keys()) | set(nonmono_data.keys()))
    
    # Combiner avec 0.0 pour les valeurs manquantes
    for z in all_z:
        dens_mono = mono_data.get(z, 0.0)
        dens_nonmono = nonmono_data.get(z, 0.0)
        print(f"{z:7.2f}\t{dens_mono:.10f}\t{dens_nonmono:.10f}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python combine_profiles.py <mono_file> <nonmono_file>", file=sys.stderr)
        print("Example: python combine_profiles.py dens_mono_avg.dat dens_nonmono_avg.dat", file=sys.stderr)
        sys.exit(1)
    
    mono_file = sys.argv[1]
    nonmono_file = sys.argv[2]
    
    combine_profiles(mono_file, nonmono_file)


if __name__ == "__main__":
    main()