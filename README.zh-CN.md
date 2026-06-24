# 欢迎使用 IPS！

<p align="center">
  <img src="assets/ips-logo.png" alt="IPS logo" width="420">
</p>

<p align="right">
  <a href="README.md">English</a> | <a href="README.zh-CN.md">中文</a>
</p>

# 快速开始

开始之前，请确认已经安装 Python 3.5 或更高版本，以及 2020 年 3 月 3 日或更新版本的 LAMMPS。

1. 确认 Python 环境

   ```shell
   python --version
   ```

   如果 Python 环境安装正确，你会看到类似下面的输出：

   ```shell
   $ python --version
   Python 3.10.6
   ```

2. 确认 LAMMPS 环境

   ```shell
   lmp -h
   ```

   如果没有报错，说明 LAMMPS 环境正常。

3. 将 LAMMPS 数据文件重命名为 `data.lmps`，并放在项目根目录下。

4. 将相关输入文件和设置文件放入 `./config` 文件夹。

5. 在项目根目录执行以下命令安装依赖：

   ```shell
   pip install -r requirements.txt
   ```

   然后执行：

   ```shell
   python Simulator.py
   ```

   程序将开始运行。

# 配置指南

配置文件使用 YAML 格式。

请注意，YAML 对空格有严格要求。例如，下面的格式是错误的：

```yaml
minfile:"config/min.in" # 错误格式：minfile: 和 "config/min.in" 之间没有空格
```

正确格式如下：

```yaml
minfile: "config/min.in" # 正确格式：minfile: 和 "config/min.in" 之间有一个空格
```

请按照 YAML 格式规范进行编辑。

## 配置参数

### 基本程序设置

`minfile`：
- 说明：能量最小化输入文件路径
- 值类型：字符串
- 示例：
  ```yaml
  minfile: "config/min.in"
  ```

`mdfile`：
- 说明：分子动力学输入文件路径
- 值类型：字符串
- 示例：
  ```yaml
  mdfile: "config/md.in"
  ```

`bondsNum`：
- 说明：当前已经形成的键数。设为 `0` 表示开始新的计算；如果设为大于 `0` 的值，程序会基于 `outPath` 目录中的数据继续计算。
- 示例：
  ```yaml
  bondsNum: 0
  ```

`bondLimits`：
- 说明：计算过程中允许形成的最大键数。达到该上限后，能量最小化计算将停止。
- 示例：
  ```yaml
  bondLimits: 500
  ```

`tryNum`：
- 说明：在当前成键循环中寻找有效成键原子对的尝试次数。如果没有找到有效原子对，程序会继续尝试，直到找到可成键原子对或达到最大尝试次数。
- 示例：
  ```yaml
  tryNum: 1
  ```

`trylimit`：
- 说明：每个成键循环中寻找有效成键原子对的最大尝试次数。如果超过该限制仍未成功，整个计算将终止。
- 示例：
  ```yaml
  trylimit: 140
  ```

`lammps`：
- 说明：从命令行执行 LAMMPS 的命令。
- 示例：
  ```yaml
  lammps: "mpirun -np 32 lmp_oneapi"
  ```

`outPath`：
- 说明：每轮迭代计算结果的输出目录。
- 示例：
  ```yaml
  outPath: "outdirs"
  ```

`lmprundir`：
- 说明：LAMMPS 运行目录。
- 示例：
  ```yaml
  lmprundir: "lammps"
  ```

### 成键条件设置

这些设置用于定义成键约束。

`link`：
- 说明：定义参与成键的原子类型。每一项由成对原子类型组成，原子类型代码必须在分子体系设置中定义。支持多个 link 组以适配不同成键场景。
- 示例：
  ```yaml
  link:
    -
      - nh,n
      - c,c_1
    -
      - tnv,tns
      - c,c_1
  ```

`deletelink`：
- 说明：成键时需要删除的原子类型，例如离去基团。
- 示例：
  ```yaml
  deletelink:
    - hn
    - cl
    - thn
  ```

`cutoff`：
- 说明：两个原子能够成键的最大距离，单位为 Å。
- 示例：
  ```yaml
  cutoff: 3.5
  ```

`charge`：
- 说明：赋给两个成键原子的虚拟电荷。第一个值和第二个值分别对应 `link` 中相同位置的原子类型。
- 示例：
  ```yaml
  charge:
    - 0
    - 0
  ```

`intra`：
- 说明：第一个参数表示是否检查两个潜在成键原子是否属于同一分子内的相邻原子；第二个参数表示该邻接关系允许间隔的最大原子数。
- 示例：
  ```yaml
  intra:
    - true
    - 1
  ```

### 平衡计算设置

当启用平衡计算时，使用以下设置。

`equilibrium`：
- 说明：是否在模拟过程中执行平衡计算。
- 值类型：布尔值
- 示例：
  ```yaml
  equilibrium: true
  ```

当 `equilibrium: true` 时，需要设置以下参数：

`upperPlate`：
- 说明：上层板的原子类型名称。
- 值类型：字符串
- 示例：
  ```yaml
  upperPlate: "up"
  ```

`upperSolvent`：
- 说明：上层溶剂对应的原子类型编号列表。
- 值类型：整数列表
- 示例：
  ```yaml
  upperSolvent:
    - 11
    - 12
  ```

`upperMonomer`：
- 说明：上层单体的 LAMMPS 数据文件路径。
- 值类型：字符串
- 示例：
  ```yaml
  upperMonomer: "config/tmc.lmps"
  ```

`upperFactor`：
- 说明：比较上层单体数量和溶剂数量时使用的比例因子。单体数量会先乘以该因子再进行比较。
- 值类型：整数
- 示例：
  ```yaml
  upperFactor: 97
  ```

`lowerPlate`：
- 说明：下层板的原子类型名称。
- 值类型：字符串
- 示例：
  ```yaml
  lowerPlate: "be"
  ```

`lowerSolvent`：
- 说明：下层溶剂对应的原子类型编号列表。
- 值类型：整数列表
- 示例：
  ```yaml
  lowerSolvent:
    - 13
    - 14
  ```

`lowerMonomer`：
- 说明：下层单体的 LAMMPS 数据文件路径。
- 值类型：字符串
- 示例：
  ```yaml
  lowerMonomer: "config/mpd.lmps"
  ```

`lowerFactor`：
- 说明：比较下层单体数量和溶剂数量时使用的比例因子。单体数量会先乘以该因子再进行比较。
- 值类型：整数
- 示例：
  ```yaml
  lowerFactor: 258
  ```

### 原子类型信息

这些参数定义 `data.lmps` 中 `atoms`、`bonds`、`angles`、`dihedrals` 和 `impropers` 对应的类型编号。对于 improper，第一个原子是中心原子。

示例：

```yaml
atoms:
  c: 1
  c3: 2
  ca: 3
  ha: 4

bonds:
  c3,c3: 1
  c3,hc: 2
  c3,n: 3
  c3,n3: 4

angles:
  c3,c3,hc: 1
  c3,c3,n3: 2
  c3,c3,n: 3
  hc,c3,hc: 4
  hc,c3,n3: 5

dihedrals:
  ca,ca,c,cl: 1
  ca,ca,c_1,n: 2
  o,c_1,n,c3: 3
  ca,c_1,n,c3: 4
  n3,c3,c3,n: 5

impropers:
  ca,o,c_1,n: 1
  ca,o,c,cl: 2
  ca,ca,ca,ha: 3
```

# 引用

如果你在研究工作中使用本软件，请引用：

```bibtex
@article{ips_2026,
  title = {{IPS}: {An} interfacial polymerization simulator for modeling the highly-cross-linked polymer membranes},
  author = {Liu, Gan and Xu, Zhaoqin and Wei, Mingjie and Huang, Jun and Du, Yan and Wang, Yong},
  journal = {Advanced Membranes},
  volume = {9},
  pages = {100252},
  year = {2026},
  doi = {10.1016/j.advmem.2026.100252}
}
```

# 许可证

本项目基于 GNU General Public License v3.0 (GPL-3.0) 授权发布。详情请见 [LICENSE](LICENSE)。
