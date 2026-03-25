# ICON 深大气 (`ldeepatmo`) 代码要点梳理（聚焦垂直坐标与垂直平流）

## 1) 开关与总体设计

- `dynamics_nml` 中通过 `ldeepatmo` 控制是否启用深大气修正；默认值是 `.FALSE.`。
- 启用后，并不是“只改一个方程项”，而是通过一组预计算的**度量修正因子**（`deepatmo_*`）统一进入动力核与平流算子。
- 代码中还做了配置约束：例如要求 `iadv_rhotheta=2`，并限制与平面/环面模式、部分 testcase 的组合。

## 2) 垂直坐标相关：几何高度与位势高度变换

`src/shared/mo_deepatmo.f90` 给出统一的高度变换接口 `deepatmo_htrafo`，支持两种方向：

- `z2zgpot`: 几何高度 `z` → 位势高度 `z_gpot`  
  \[
  z_{gpot} = \frac{z}{1+z/a}
  \]
- `zgpot2z`: 位势高度 `z_gpot` → 几何高度 `z`  
  \[
  z = \frac{z_{gpot}}{1-z_{gpot}/a}
  \]

其中 `a` 是地球半径。该模块同时有 OpenMP/OpenACC 实现，以便在 CPU/GPU 上都能调用。

## 3) 深大气度量因子的核心思想

在 `mo_vertical_grid.f90` 中，ICON 会基于层顶/层底半层面半径（`r_u`, `r_l`）和层中心半径（`r_mc`）构造深大气修正：

- 水平梯度修正：`deepatmo_gradh_mc = a / r_mc`
- 水平散度修正：`deepatmo_divh_mc`（与 `a/r` 接近但不完全相同）
- 垂直散度上下界面权重：`deepatmo_divzL_mc`, `deepatmo_divzU_mc`
- 控制体体积修正：`deepatmo_vol_mc`
- 逆半径：`deepatmo_invr_mc`, `deepatmo_invr_ifc`

这些因子都是**按垂直层的一维数组**（仅随 `jk` 变），这是一个关键工程设计：

- 注释明确说明：为节省内存并控制几何复杂性，深大气修正忽略 `flat_height` 以下坐标面的地形依赖，避免把这些因子变成 3D 场。

## 4) 你关心的“每个单元控制体随高度膨胀”如何进入方程

核心不是直接改变网格连通关系，而是通过 `deepatmo_vol_mc` 与散度系数进入离散守恒式：

- `deepatmo_vol_mc(jk)` 随半径增大而变化，表示同一“地表归一化面积”对应的层体积比例变化。
- 在质量/空气质量诊断中，单元质量按
  \[
  \text{airmass} = \rho \cdot \Delta z \cdot \text{deepatmo\_vol\_mc}(jk)
  \]
  计算，从而体现高空控制体膨胀。

## 5) 垂直平流离散：深大气如何改写通量散度

在 `mo_advection_stepping.f90` 里，垂直质量与示踪物更新都使用：

\[
\Delta(\cdot) \propto F_{k+1/2}\,\text{deepatmo\_divzL\_mc}(k)
           - F_{k-1/2}\,\text{deepatmo\_divzU\_mc}(k)
\]

对应代码层面体现为：

- 先更新中间 `rhodz`（考虑垂直通量贡献）；
- `vert_adv` 中按上下界面通量与 `deepatmo_divz{L,U}_mc` 组合更新 tracer；
- 最后除以更新后的 `rhodz_new`，保证守恒形式。

这相当于把“球壳几何中上下界面面积与体积比”的效应压缩进两个垂直系数，避免每步显式做复杂几何积分。

## 6) 设计总结（面向实现）

- **统一入口开关**：`ldeepatmo` 在 namelist 控制，默认关闭。
- **统一几何变换模块**：`mo_deepatmo` 管理 `z` 与 `z_gpot` 互转。
- **统一度量因子注入**：`deepatmo_*` 在垂直网格阶段预计算，然后在动力与平流中复用。
- **守恒优先**：通过“通量散度系数 + 体积修正”表达控制体膨胀，而不是引入昂贵 3D 几何因子场。
- **性能导向**：保留 1D 垂直系数（而非 3D）是明显的内存/算力折中。

## 7) 额外提示：配置语义演进

`doc/Namelist_overview/incompatible_changes.tex` 指出：

- 早期一些深大气相关细分开关已移除；
- 现在 `dynamics_nml/ldeepatmo=.TRUE.` 隐含采用非传统项与非常量重力等约定语义。

这意味着：当前代码路径倾向“以 `ldeepatmo` 为总开关 + 内置物理假设”，而非让用户单独拼装多个子开关。
