# HPC

超算用脚本

## 脚本列表

```bash
sruns                                # 根据超算类型申请计算节点
time_cost                            # 计算 HPC 提交任务运行耗时
slurm_generation.py                  # 生成不同平台的 VASP/LAMMPS/Python/Bash 任务 Slurm 提交脚本
gpu_slurm_generation.py              # 生成 GPU 任务的 Slurm 提交脚本（课题组服务器平台）
job_history.sh                       # 查看提交至 Slurm 队列系统的 job 历史
running_jobid_info.sh                # 获取正在运行的任务 JobId 信息
running_jobid_info2.sh               # 获取正在运行的任务 JobId 信息（针对 sacct 无法使用的情况）
completed_jobid_info.sh              # 获取已完成任务 JobId 信息
```
