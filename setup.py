from setuptools import setup, find_packages

setup(
    name='pymulti',  # 你的包名
    version='0.1.0',  # 初始版本号
    author='Zixu Wang',  # 作者名
    author_email='wangzx2022@sjtu.edu.cn',  # 作者邮箱
    description='PyMulti is developed in order to make Multi Code available and convenient on Python and Pytorch.',  # 简短描述
    long_description=open('README.md').read(),  # 从README文件读取长描述
    long_description_content_type='text/markdown',  # README文件的格式
    url='https://github.com/xuxu20040407/PyMulti',  # 项目主页，通常是GitHub仓库地址
    packages=find_packages(),  # 自动查找并包括所有包
    install_requires=[
        'numpy',
        'os',
        'multiprocessing',
        'subprocess',
        'shutil'
    ],
)