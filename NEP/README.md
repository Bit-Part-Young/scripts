# NEP

NEP 相关脚本

## 脚本列表

```bash
count_xyz.sh                         # 统计 xyz 文件构型帧数及总原子数
get_xyz_efs.sh                       # 获取 NEP xyz 文件中的能量、力、位力数据
thermo_info_gpumd.sh                 # 从 GPUMD 输出的 thermo.out 文件中格式化输出热力学信息

get_mp_surface_gb.py                 # 从 MP 获取 表面性质 和 晶界 数据并保存为 json 文件
json2surfaces.py                     # 从 json 文件中提取 MP 表面构型并保存为 xyz 文件
json2gb.py                           # 从 json 文件中提取 MP 晶界构型并保存为 xyz 文件

xdatcar2xyz.py                       # 将 VASP XDATCAR 文件转换为 NEP xyz 文件，并指定间隔抽样
vaspout2xyz.py                       # 将 VASP 输出文件 OUTCAR/vasprun.xml 转换为 xyz 文件
xyz2cfg.py                           # 将 NEP xyz 文件转换成 MTP cfg
outcar2xyz_singleframe.sh            # 将 OUTCAR（静态计算，单帧） 转换为 NEP xyz 文件
outcar2xyz_multipleframes.sh         # 将 OUTCAR（弛豫计算，多帧） 转换为 NEP xyz 文件

time_consuming_gpumd.sh              # gpumd 程序运行耗时监控
time_consuming_nep.sh                # nep 程序运行耗时监控

relax.py                             # 结构弛豫 ASE 实现
eos_cal_calorine.py                  # 使用 GPUMD & calorine 进行 EOS 计算
nep_pca_plot.py                      # NEP 未归一化 & 归一化 描述符 PCA 2 维绘制

relax_nep.py                         # NEP 势函数 结构优化
elastic_nep.py                       # NEP 势函数 弹性常数计算
eos_nep.py                           # NEP 势函数 EOS 计算

fps_select.py                        # 基于最远点采样选择结构
perturb_hiphive.py                   # 用 dpdata 生成随机 cell & position 扰动结构
perturb_dpdata.py                    # 使用 hiphive 生成扰动结构（原子位置 扰动）


rmse_cal.py                          # 计算势函数预测与DFT 计算的能量、力、应力、位力指标的 RMSE
plot_nep_loss.py                     # 绘制 NEP 训练 loss 演化与能量、力、virial/应力 RMSE 指标 DFT计算值 NEP 预测值对比图及其对应误差直方图
plt_nep_train.py                     # 绘制 NEP 训练 loss 演化与能量、力、virial/应力 RMSE 指标 DFT计算值 NEP 预测值对比图

json2xyz.py                          # 将 json 构型及其数据文件转换为 xyz 文件
son2db.py                            # 将 json 构型及其数据文件转换为 ASE db 格式
json2df.py                           # 将 json 构型及其数据文件转换为 Pandas DataFrame 格式
```
