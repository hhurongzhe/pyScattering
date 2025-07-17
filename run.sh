#!/bin/bash
#SBATCH --partition c128m1024
#SBATCH -o output/log-%j.out
#SBATCH -J test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16


PYTHON_BIN="/tns/hrz/local_venv/bin/python3.12"
export PYTHONPATH="/tns/hrz/local_venv/lib/python3.12/site-packages"

$PYTHON_BIN -u examples/renorm_stochastic_single_channel.py
# $PYTHON_BIN -u examples/renorm_stochastic_coupled_channel.py
# $PYTHON_BIN -u examples/renorm_stochastic_couple_channel_E0.py
