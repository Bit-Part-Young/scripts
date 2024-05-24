"""
图片多余空白裁切，ChatGPT 生成
"""

from PIL import Image, ImageChops


def trim_whitespace(
    input_image: str,
    output_image: str,
):
    """图片多余空白裁切"""

    # 加载图像
    image = Image.open(input_image)

    # 将图像转换为灰度图像
    grayscale_image = image.convert("L")

    # 反转图像的颜色
    inverted_image = ImageChops.invert(grayscale_image)

    # 计算边界框
    bbox = inverted_image.getbbox()

    # 裁切图像
    cropped_image = image.crop(bbox)

    # 保存裁切后的图像
    cropped_image.save(output_image)
