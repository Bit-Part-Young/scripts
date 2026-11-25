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
hpcsync                              # 本地目录 <-> 服务器目录同步（link/where/push/pull），映射记录在 .hpcsync
rsync_upload / rsync_download        # 手写本地与远程路径的上传/下载
```

## hpcsync

本地目录与服务器目录的对应关系写在本地目录的 `.hpcsync`（一行 `SERVER:REMOTE_DIR`），子目录自动沿用最近祖先的映射。
`pull` 默认排除 `hpcsync.exclude` 中的 VASP 大文件，`--all` 取消。不使用 `--delete`。

```bash
cd ~/scripts/proj/DFT-calculations/elastic-202611
hpcsync link Master:~/work/proj/elastic-202611   # 每个项目目录登记一次
hpcsync where refine                             # -> Master:work/proj/elastic-202611/refine
hpcsync push refine                              # 上传
hpcsync pull -n refine                           # 预演
hpcsync pull refine                              # 下载
```
