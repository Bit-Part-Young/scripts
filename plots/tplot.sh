#!/bin/bash

# 使用 gnuplot 简单快速绘图


#------------------------- Get help ------------------------
usage() {
    script_name=$(basename $0)
    cat <<EOF

Plot data simply and quickly with gnuplot.

Usage: ${script_name} [OPTIONS]

Options:
    -i, --input FILE          input data filename
    -o, --output FILE         output figure filename (optional, default: output.png)
    -d, --delimiter STR       data file delimiter (e.g., ",", "\\t" for tab, " "; optional, default: whitespace)
    -pstyle, --pstyle STYLE   plot style: lines, points, linespoints (optional, default: lines)
    -x, --xcol N              x axis data column number (optional)
    -y, --ycol N              y axis data column numbers (comma-separated; optional, default: 2)
    -l, --label STR           legend labels (comma-separated; optional)
    -xr, --xrange MIN:MAX     x axis range (optional)
    -yr, --yrange MIN:MAX     y axis range (optional)
    -xl, --xlabel STR         x axis label (optional)
    -yl, --ylabel STR         y axis label (optional)
    -xlog, --xlog             set x logscale (optional)
    -xtics, --xtics N         x axis major tic interval (optional)
    -ytics, --ytics N         y axis major tic interval (optional)
    -xmtics, --xmtics N       x axis minor tic interval (optional)
    -ymtics, --ymtics N       y axis minor tic interval (optional)
    -h, --help                show this help message and exit

Examples:
    ${script_name} -i data.txt -y 2 -o plot.png
    ${script_name} -i data.txt -x 1 -y 2
    ${script_name} -i data.csv -d "," -y 2,3,4 -l "Col2,Col3,Col4" -pstyle lp
EOF
    exit 1
}


#------------------------- Default values ------------------------
INPUT_FILE=""
OUTPUT_FILE="output.png"
PSTYLE="lines"
DELIMITER=""
LABELS=()  # 图例标签数组
XCOL=""  # 空值表示使用行号
YCOLS=()
XRANGE=""
YRANGE=""
XLABEL=""
YLABEL=""
XLOG=false
XTICS=""
YTICS=""
XMTICS=""
YMTICS=""


#------------------------- Parse command line arguments ------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -pstyle|--pstyle)
            PSTYLE="$2"
            shift 2
            ;;
        -d|--delimiter)
            DELIMITER="$2"
            shift 2
            ;;
        -x|--xcol)
            XCOL="$2"
            shift 2
            ;;
        -y|--ycol)
            # 支持逗号分隔的多个列号，或单个列号
            IFS=',' read -ra COLS <<< "$2"
            for col in "${COLS[@]}"; do
                YCOLS+=("${col}")
            done
            shift 2
            ;;
        -l|--label)
            # 支持逗号分隔的多个标签，或单个标签
            IFS=',' read -ra LABEL_ARRAY <<< "$2"
            for label in "${LABEL_ARRAY[@]}"; do
                LABELS+=("${label}")
            done
            shift 2
            ;;
        -xr|--xrange)
            XRANGE="$2"
            shift 2
            ;;
        -yr|--yrange)
            YRANGE="$2"
            shift 2
            ;;
        -xl|--xlabel)
            XLABEL="$2"
            shift 2
            ;;
        -yl|--ylabel)
            YLABEL="$2"
            shift 2
            ;;
        -xlog|--xlog)
            XLOG=true
            shift
            ;;
        -xtics|--xtics)
            XTICS="$2"
            shift 2
            ;;
        -ytics|--ytics)
            YTICS="$2"
            shift 2
            ;;
        -xmtics|--xmtics)
            XMTICS="$2"
            shift 2
            ;;
        -ymtics|--ymtics)
            YMTICS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done


# 检查必需参数
if [[ -z "${INPUT_FILE}" ]]; then
    echo "Error: input filename is required!"
    usage
fi

# 如果没有指定y列，使用默认值2
if [[ ${#YCOLS[@]} -eq 0 ]]; then
    YCOLS=(2)
fi

# 检查输入文件是否存在
if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "Error: input file '${INPUT_FILE}' does not exist!"
    exit 1
fi

# 解析范围参数
XMIN=""
XMAX=""
YMIN=""
YMAX=""

if [[ -n "${XRANGE}" ]]; then
    IFS=':' read -r XMIN XMAX <<< "${XRANGE}"
fi

if [[ -n "${YRANGE}" ]]; then
    IFS=':' read -r YMIN YMAX <<< "${YRANGE}"
fi


#------------------------- Generate gnuplot script ------------------------
if hostname | grep -q sjtu; then
  config_path="~/yangsl/scripts/cms-scripts/plots"
else
  config_path="~/scripts/cms-scripts/plots"
fi

gnuplot_fn=".tmp.gnu"

cat > ${gnuplot_fn} <<EOF
set loadpath "${config_path}"
load "config.gnu"

set output "${OUTPUT_FILE}"
EOF


# 设置数据文件分隔符
if [[ -n "${DELIMITER}" ]]; then
    # 处理特殊字符转义序列
    if [[ "${DELIMITER}" == "\\t" ]]; then
        echo -e "\nset datafile separator \"\\t\"" >> ${gnuplot_fn}
    elif [[ "${DELIMITER}" == "\\s" ]]; then
        echo -e "\nset datafile separator \" \"" >> ${gnuplot_fn}
    else
        # 转义引号以便在 gnuplot 中正确使用
        DELIMITER_ESC=$(printf '%s' "${DELIMITER}" | sed 's/"/\\"/g')
        echo -e "\nset datafile separator \"${DELIMITER_ESC}\"" >> ${gnuplot_fn}
    fi
fi

if [[ -n "${XLABEL}" ]]; then
    echo -e "\nset xlabel \"${XLABEL}\"" >> ${gnuplot_fn}
fi

if [[ -n "${YLABEL}" ]]; then
    echo -e "\nset ylabel \"${YLABEL}\"" >> ${gnuplot_fn}
fi

if [[ "${XLOG}" == true ]]; then
    echo -e "\nset logscale x" >> ${gnuplot_fn}
fi

if [[ -n "${XRANGE}" ]]; then
    echo -e "\nset xrange [${XMIN}:${XMAX}]" >> ${gnuplot_fn}
fi

if [[ -n "${YRANGE}" ]]; then
    echo -e "\nset yrange [${YMIN}:${YMAX}]" >> ${gnuplot_fn}
fi

if [[ -n "${XTICS}" ]]; then
    echo -e "\nset xtics ${XTICS}" >> ${gnuplot_fn}
fi

if [[ -n "${YTICS}" ]]; then
    echo -e "\nset ytics ${YTICS}" >> ${gnuplot_fn}
fi

if [[ -n "${XMTICS}" ]]; then
    echo -e "\nset mxtics ${XMTICS}" >> ${gnuplot_fn}
fi

if [[ -n "${YMTICS}" ]]; then
    echo -e "\nset mytics ${YMTICS}" >> ${gnuplot_fn}
fi


# 确定 x 轴数据源：若未指定 XCOL，则使用行号 ($0+1)
if [[ -z "${XCOL}" ]]; then
    XDATA="(\$0+1)"
else
    XDATA="${XCOL}"
fi

# 生成 plot 命令，支持标签和样式
if [[ ${#YCOLS[@]} -eq 1 ]]; then
    # 只有一个 y 列
    YCOL="${YCOLS[0]}"
    LABEL_STR=""
    if [[ ${#LABELS[@]} -gt 0 ]]; then
        LABEL_STR=" t \"${LABELS[0]}\""
    fi
    echo -e "\nplot \"${INPUT_FILE}\" u ${XDATA}:${YCOL} w ${PSTYLE} ${LABEL_STR}" >> ${gnuplot_fn}
else
    # 多个 y 列，使用反斜杠换行
    echo -e "\nplot \\" >> ${gnuplot_fn}
    for i in "${!YCOLS[@]}"; do
        YCOL="${YCOLS[${i}]}"
        LABEL_STR=""
        if [[ $i -lt ${#LABELS[@]} ]]; then
            # 使用提供的标签
            LABEL_STR="t \"${LABELS[${i}]}\""
        else
            # 如果没有提供标签或标签数量不足，使用默认标签
            LABEL_STR="t \"Col ${YCOL}\""
        fi

        if [[ $i -lt $((${#YCOLS[@]} - 1)) ]]; then
            # 不是最后一个，需要逗号和反斜杠
            echo "     \"${INPUT_FILE}\" u ${XDATA}:${YCOL} w ${PSTYLE} ${LABEL_STR},\\" >> ${gnuplot_fn}
        else
            # 最后一个，不需要逗号
            echo "     \"${INPUT_FILE}\" u ${XDATA}:${YCOL} w ${PSTYLE} ${LABEL_STR}" >> ${gnuplot_fn}
        fi
    done
fi

echo "" >> ${gnuplot_fn}


gnuplot ${gnuplot_fn}


echo -e "\n${OUTPUT_FILE} is generated!"

echo -e "\nNote: The plot settings are saved in ${gnuplot_fn}, you can further modify it!"
