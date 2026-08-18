"""PyInstaller 打包入口(顶层脚本:无相对导入)。

直接以 blastprime/app.py 作为 PyInstaller 入口会被当作顶层脚本执行,
其中 `from .config import ...` 等相对导入会报
'attempted relative import with no known parent package';
本文件以包内 main() 方式启动,打包/源码运行语义一致。

Top-level entry for PyInstaller (no relative imports). Pointing the
packager at blastprime/app.py executes it as a top-level script, which
breaks its relative imports; this file starts main() through the package
so packaged and source runs behave identically.
"""
from blastprime.app import main

if __name__ == "__main__":
    main()
