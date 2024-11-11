Pymulti 是基于Multi程序设计的接口，旨在方便批处理与数据库建立。

# 文件结构
- pymulti： 源代码
  - __init__.py： 初始化文件
  - multi_1D.py： 有关1D的方法
  - multi_2D.py： 有关2D的方法
  - process_2D.py： 有关2D数据的处理方法
  - sample_method.py： 有关采样的方法
- source： 初始化算例的源文件
  - 1D： 初始化1D算例的源文件
    - fort.12：输入模板
    - multi：运行程序
  - 2D： 初始化2D算例的源文件
    - multi2d：运行程序
    - User.r： 输入模板
  - tables： 运行1D程序必备的文件
- main.py： 示例代码