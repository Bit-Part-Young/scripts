#!/bin/bash

# 对个人常用的 Git 仓库进行 pull 操作
# 在不同电脑、服务器上使用该脚本时
# 需自定义修改 pull_folders 和 pull_path

pull_folders=(
  "obsidian-md"
  "cms-scripts"
  "dotfiles"
  "dotfiles2"
  "reference-notes"
  "VASP-Exercise"
  "atomate-test"
  "Ti-project"
  "Nb-Si-data"
  "draft-Nb-Si"
  "group-meeting"
  "hexo-butterfly-demo"
  "mkdocs-demo"
)

for folder in "${pull_folders[@]}"; do
  pull_path="${HOME}/scripts/${folder}"

  echo "Starting git pull in ${folder}..."

  if [[ ${folder} == "group-meeting" ]]; then
    pull_path="${HOME}/Documents/${folder}"
  fi

  if [[ ${folder} == "hexo-butterfly-demo" || ${folder} == "mkdocs-demo" ]]; then
    pull_path="${HOME}/scripts/sites/${folder}"
  fi

  # 在子 Shell 中执行，并放到后台运行
  (
    cd ${pull_path}
    # 获取当前仓库 branch 名称
    current_branch=$(git branch --show-current)

    git pull origin ${current_branch}

    echo "Completed git pull in ${folder}"
  ) &

done

wait # 等待所有后台任务完成

echo
echo "All git pulls are completed."
