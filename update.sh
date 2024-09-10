# 有时候需要开启代理
git config https.proxy 127.0.0.1:1086

# 关于添加子模块的参考：https://knightyun.github.io/2021/03/21/git-submodule

# 添加 htslib 作为子模块
git submodule add git@github.com:samtools/htslib.git

# 并初始化更新(含递归)
# git clone下来的工程中带有submodule时，初始的时候，
# submodule的内容并不会自动下载下来的，此时，只需执行如下
# 命令：`git submodule update --init --recursive`
# 即可将子模块内容下载下来后工程才不会缺少相应的文件。
git submodule update --init --recursive

# 更新子模块（二选一）
## 方法 1：进入子模块目录执行 `git pull`; 
## 方法 2：直接在父目录下执行 `git submodule update --remote` 即可
git submodule update --remote --recursive





