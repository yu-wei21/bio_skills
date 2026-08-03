# 项目规则

## 运行环境

所有 R / Python 包及脚本的运行环境统一为：

```
/home2/wei_yu/.conda/envs/seurat
```

运行 R 脚本时使用该环境中的 R（例如 `Rscript` 通过该 conda 环境调用），安装或使用任何 R / Python 包均以该环境为准，避免使用系统默认环境或其它 conda 环境。
