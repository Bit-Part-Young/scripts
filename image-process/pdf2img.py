"""
pdf 转图片；工具：https://github.com/Belval/pdf2image

"""

import os
import sys

from pdf2image import convert_from_bytes, convert_from_path
from pdf2image.exceptions import (
    PDFInfoNotInstalledError,
    PDFPageCountError,
    PDFSyntaxError,
)


def pdf2img(
    pdf_file: str,
    output_file: str,
    dpi: int = 300,
    single_file: bool = False,
    output_folder: str = "./fig-pdf",
):
    """单个多页面 pdf 文件转图片"""

    os.makedirs(output_folder, exist_ok=True)

    images = convert_from_path(
        pdf_path=pdf_file,
        dpi=dpi,
        fmt="jpg",  # jpg 大小比 png 小
        single_file=single_file,
        output_folder=output_folder,
        output_file=output_file,
    )

    return images


if __name__ == "__main__":
    pdf_file = sys.argv[1]
    output_folder = sys.argv[2]

    pdf2img(
        pdf_file=pdf_file,
        output_file=pdf_file[:-4],
        output_folder=output_folder,
    )

    print("\npdf to image conversion is done.")
