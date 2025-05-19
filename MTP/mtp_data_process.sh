#!/bin/bash

# 将 MTP 构型的 DFT 和 预测数据（能量、受力、应力）进行后处理


get_xyz_efs.sh train.xyz dft
echo -e "\n"
get_xyz_efs.sh train_predict.xyz predict

echo -e "\nPredict and DFT data combined to *_train.out"

paste -d " " energy_predict.dat energy_dft.dat > energy_train.out
paste -d " " force_predict.dat force_dft.dat > force_train.out
paste -d " " virial_predict.dat virial_dft.dat > stress_train.out
