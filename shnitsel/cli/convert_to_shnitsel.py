

import argparse
import sys

def main():
    args = argparse.ArgumentParser(sys.argv[0], f"{sys.argv[0]} <input_path> [OPTIONS]", 
                                   description="Script to read in an individual trajectory or a directory containing multiple " \
                                   "sub-trajectories and convert them into a shnitsel-style ")

if __name__ == "__main__":
    main()