"""tiff 格式图片压缩测试；效果一般"""

from PIL import Image


def compress_tiff(
    input_fn: str,
    compression: str = "tiff_adobe_deflate",
):
    with Image.open(input_fn) as img:
        img.save(
            f"{input_fn[:-4]}_{compression}.tif",
            compression=compression,
        )


if __name__ == "__main__":
    copmress_method_list = [
        "tiff_adobe_deflate",  # 压缩效果最好 减小 25%
        "tiff_deflate",  # 同上
        "tiff_lzw",  # 其次
        "tiff_packbits",  # 基本不压缩
        "tiff_lzmam",  # 基本不压缩
    ]

    input_path = "GND.tif"
    for copmress_method in copmress_method_list:
        compress_tiff(
            input_fn=input_path,
            compression=copmress_method,
        )
        print(f"tiff compress method: {copmress_method} done!")
