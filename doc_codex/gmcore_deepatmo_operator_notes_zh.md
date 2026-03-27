# GMCORE 深大气算子修改说明（当前实现，重点：垂直平流）

本文档对应当前 GMCORE 工作树中的深大气相关算子修改，重点解释：

- `deepwater` 打开后，哪些算子进入深大气路径
- 为什么垂直平流要改成带 \((r/a)^2\) 的通量散度
- 水平梯度、PGF、旋度等配套几何项目前是怎样落到代码里的

这里的目标是说明“当前代码实现为什么这么写”，不是完整复述论文推导。

## 1. 总开关与子开关

当前实现把深大气当作一组联动开关处理：

- `deepwater = .true.`：启用深大气框架
- 若用户没有单独指定子开关，则自动打开：
  - `use_mesh_change`
  - `use_vert_nct`
  - `use_hor_nct`
- `deepwater = .false.` 时，上述三个子开关都会被关掉

代码位置：

- [`src/utils/namelist_mod.F90:423`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/utils/namelist_mod.F90#L423)

这意味着当前代码里“几何尺度变化”主要由 `deepwater .and. use_mesh_change` 这一路条件控制。

## 2. 局地半径场 `rdp*` 的含义

深大气几何修正建立在局地半径场 `rdp*` 上：

- `rdp`：格心半径
- `rdp_lon` / `rdp_lat`：水平边上的半径
- `rdp_lev`：半层面半径
- `rdp_lev_lon` / `rdp_lev_lat`：半层面的水平边半径
- `rdp_vtx`：顶点半径

这些量由 `calc_rdp` 统一计算，再插值/平均到不同位置：

- [`src/dynamics/operators_mod.F90:1235`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L1235)

实现上，`rdp` 不是简单的 `a + z`，而是通过柱质量和密度积分得到的等效半径。随后：

- `average_run(rdp, rdp_lon)`
- `average_run(rdp, rdp_lat)`
- `interp_run(rdp, rdp_lev)`

见：

- [`src/dynamics/operators_mod.F90:1284`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L1284)

## 3. 为什么垂直平流必须带 \((r/a)^2\)

### 3.1 物理原因

浅大气近似默认同一柱体上下界面的水平面积不变，因此垂直通量散度常写成：

\[
-\left(F_{k+1/2} - F_{k-1/2}\right)
\]

但深大气里，柱体随半径增大而膨胀。若水平面积按参考半径 \(a\) 归一化，则真实球壳面积满足：

\[
A(r) \propto r^2
\]

因此保守的垂直通量散度应写成：

\[
- \frac{F_{k+1/2} A_{k+1/2} - F_{k-1/2} A_{k-1/2}}{A_k}
\]

若以参考半径 \(a\) 归一化，

\[
A_k \propto \left(\frac{r_k}{a}\right)^2
\]

就得到当前代码使用的形式：

\[
- \frac{
F_{k+1/2}\left(r_{k+1/2}/a\right)^2 -
F_{k-1/2}\left(r_{k-1/2}/a\right)^2
}{
\left(r_k/a\right)^2
}
\]

这一步的本质不是“额外乘了个经验因子”，而是把球壳面积随高度变化显式写进了垂直守恒式。

### 3.2 工程含义

如果不加这组面积权重：

- 高层球壳面积变大这一事实不会进入垂直散度
- 质量/示踪物柱积分守恒会和深大气几何不一致
- 即使水平算子改成深大气形式，垂直输送仍然会保留浅大气体积假设

所以从一致性上说，`use_mesh_change` 一旦负责“几何尺度变化”，垂直平流就不能只改水平、不改竖直。

## 4. 垂直平流在代码里的三处核心落点

## 4.1 质量垂直平流：`adv_run_mass`

质量更新原来的浅大气形式是：

\[
m^{n+1} = m^{*} - \Delta t \left(mfz_{k+1/2} - mfz_{k-1/2}\right)
\]

当前深大气实现改成：

\[
m^{n+1} = m^{*} - \Delta t \,
\frac{
mfz_{k+1/2}(r_{k+1/2}/a)^2 -
mfz_{k-1/2}(r_{k-1/2}/a)^2
}{
(r_k/a)^2
}
\]

代码位置：

- [`src/dynamics/adv/adv_mod.F90:272`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_mod.F90#L272)
- [`src/dynamics/adv/adv_mod.F90:281`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_mod.F90#L281)

实现要点：

- 先算水平通量并更新 `m_new`
- 再算垂直通量 `mfz`
- 深大气路径下，对上下界面通量分别乘 `(rdp_lev/radius)^2`
- 最后再除以本层 `(rdp_lev/radius)^2`

这就是标准的“界面面积加权后，再除以控制体参考面积”的通量形式。

## 4.2 示踪物垂直平流：`adv_run_tracers`

示踪物更新不是直接对 `q` 做差分，而是先对示踪物质量通量 `qmfz` 做守恒更新，再除以更新后的质量 `m_new`：

\[
q^{n+1} =
\frac{
m^{xy} q^{xy} - \Delta t \, \nabla_z \left(qmfz\right)_{\text{deep}}
}{
m^{n+1}
}
\]

其中深大气垂直散度写成：

\[
\nabla_z(qmfz)_{\text{deep}}
=
\frac{
qmfz_{k+1/2}(r_{k+1/2}/a)^2 -
qmfz_{k-1/2}(r_{k-1/2}/a)^2
}{
(r_k/a)^2
}
\]

代码位置：

- [`src/dynamics/adv/adv_mod.F90:361`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_mod.F90#L361)
- [`src/dynamics/adv/adv_mod.F90:371`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_mod.F90#L371)

这一步和质量更新保持完全同构，因此质量和示踪物的垂直输送是一致的。

## 4.3 半层量垂直平流：`calc_adv_lev`

`nh_mod` 里的 `calc_adv_lev` 用于半层量（例如 `w_lev`, `gz_lev` 一类）平流倾向计算。当前实现先去掉水平散度部分，再去掉垂直散度部分。

浅大气形式：

\[
dq/dt \leftarrow dq/dt -
\left[
qmfz_k - qmfz_{k-1} - q_k (mfz_k - mfz_{k-1})
\right]
\]

深大气形式改成：

\[
dq/dt \leftarrow dq/dt -
\frac{
qmfz_k A_k - qmfz_{k-1} A_{k-1}
- q_k \left(mfz_k A_k - mfz_{k-1} A_{k-1}\right)
}{
A_k
}
\]

其中 \(A_k \propto (r_k/a)^2\)。

代码位置：

- [`src/dynamics/nh_mod.F90:168`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/nh_mod.F90#L168)
- [`src/dynamics/nh_mod.F90:178`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/nh_mod.F90#L178)

这一段尤其重要，因为它说明你现在不是只改了“质量柱”的垂直通量，而是把半层量的垂直平流倾向也同步改到了同一套几何守恒框架下。

## 5. 与垂直动量方程的关系

垂直平流之外，`implicit_w_solver` 里还引入了另一组几何因子：

\[
factor\_r = (a/r)^2
\]

代码位置：

- [`src/dynamics/nh_mod.F90:294`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/nh_mod.F90#L294)

这里的用途和“垂直平流面积修正”不一样：

- 垂直平流里用的是 \((r/a)^2\)，对应界面/控制体面积
- `implicit_w_solver` 里用的是 \((a/r)^2\)，用于竖直动量方程中的重力项、压强梯度项等几何系数重写

不要把这两类平方因子混成一个来源。它们都和 \(r^2\) 有关，但分别来自：

- 通量散度中的球壳面积
- 竖直动量方程的深大气 metric form

## 6. 水平算子的配套修改

虽然这份文档重点是垂直平流，但当前实现还有几处必须配套理解的水平算子变化。

### 6.1 动能梯度

`calc_grad_ke` 在深大气路径下对水平梯度乘一次 `a/r`：

- [`src/dynamics/operators_mod.F90:1436`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L1436)
- [`src/dynamics/operators_mod.F90:1452`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L1452)

原因是 `mesh%de_lon` / `mesh%de_lat` 基于参考半径 `a` 定义，而真实水平长度在半径 `r` 上，因此要补一个 `a/r`。

### 6.2 深大气 Hamilton 型几何修正

`deep_hamiton_modify_3d` 对不同位置的场统一做：

\[
fd \leftarrow fd / (a/r) = fd \cdot (r/a)
\]

代码位置：

- [`src/dynamics/operators_mod.F90:1313`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L1313)

这类操作更接近“把原本按参考球定义的量换算到实际半径壳面上”。

### 6.3 旋度计算

当前修正里还把深大气旋度调用的纬向半径参数纠正成了 `rdp_lat`，而不是误传 `rdp_lon`：

- [`src/dynamics/operators_mod.F90:792`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L792)

这属于实现一致性修复。

## 7. 扰动 PGF 的当前处理

`pgf_ptb_mod` 当前把四项分成两类：

- 压强梯度族：`tmp1`, `tmp2`
- 位势梯度族：`tmp3`, `tmp4`

目前实现为：

- `tmp1`, `tmp2` 乘 `factor_r = a/r`
- `tmp3`, `tmp4` 除 `factor_r = a/r`，即净效果乘 `r/a`

代码位置：

- [`src/dynamics/pgf/pgf_ptb_mod.F90:170`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/pgf/pgf_ptb_mod.F90#L170)
- [`src/dynamics/pgf/pgf_ptb_mod.F90:175`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/pgf/pgf_ptb_mod.F90#L175)
- [`src/dynamics/pgf/pgf_ptb_mod.F90:179`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/pgf/pgf_ptb_mod.F90#L179)
- [`src/dynamics/pgf/pgf_ptb_mod.F90:197`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/pgf/pgf_ptb_mod.F90#L197)
- [`src/dynamics/pgf/pgf_ptb_mod.F90:202`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/pgf/pgf_ptb_mod.F90#L202)
- [`src/dynamics/pgf/pgf_ptb_mod.F90:206`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/pgf/pgf_ptb_mod.F90#L206)

当前判断是：

- 位势相关两项不需要再改
- 压强相关两项已经补到了和深大气水平梯度一致的 `a/r`

## 8. 当前实现可概括为一条主线

如果只看“几何尺度变化”这一条线，当前工作树的实现逻辑可以压缩成：

1. 用 `calc_rdp` 生成局地半径场 `rdp*`
2. 水平一阶梯度相关项按需要补 `a/r`
3. 位势梯度类项按当前 metric form 处理成净 `r/a`
4. 垂直通量散度统一改成带 \((r/a)^2\) 的面积加权形式
5. 竖直动量方程另外使用 \((a/r)^2\) 的 metric 系数

这五步合起来，才是当前代码里“深大气几何改动”的完整闭环。

## 9. 目前最值得继续盯的点

从实现一致性角度，后续最值得继续检查的是：

- `mfz` / `qmfz` 的物理量纲和几何归一化定义是否和面积权重完全一致
- `calc_adv_lev` 与 `adv_run_tracers` 是否在所有半层量上都应共享同一套几何形式
- `implicit_w_solver` 中 \((a/r)^2\) 的推导是否和你当前采用的 Wood & Staniforth 写法完全同号、同归一化
- 若后续把重力改成随半径变化的非常量 `g(r)`，垂直几何与动力项还需再做一次一致性复核

## 10. 一句话结论

当前这版修改里，最关键的变化不是某一个 `factor_r` 的指数，而是：

**垂直平流已经从浅大气的“常面积柱体通量差”，改成了深大气的“球壳面积随半径变化的守恒通量散度”。**

这一步是整套几何修正里最核心、也最不应回退的一步。

## 11. 水平平流链条的本次补齐

这次除了垂直平流，还把水平平流里原先残留的 shallow 几何入口补齐了。核心思想只有一条：

- 质量通量先按局地物理量定义
- 水平 CFL、跨格质量、水平散度再统一消费局地半径几何

这样做是为了避免同一个 `a/r` 同时在“通量构造”和“散度算子”里各乘一次。

### 11.1 质量通量 `mfx/mfy` 改回局地通量

`calc_mf` 现在使用：

\[
mfx = dmg_{lon} \cdot u,\qquad mfy = dmg_{lat} \cdot v
\]

不再在这里预乘 `a/r`。

代码位置：

- [`src/dynamics/operators_mod.F90:727`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L727)

这样后面的 deep 散度只需要吃一次几何修正，语义上也更接近 ICON 中“质量通量本身”和“divergence metric”分离的结构。

### 11.2 水平 CFL 改用局地格距

`adv_batch_calc_cflxy_mass` 现在不再直接用参考球上的 `mesh%de_lon/de_lat`，而是改用：

- `adv_batch_edge_dx`
- `adv_batch_edge_dy`

它们在 deep 路径下分别返回：

\[
\Delta x_{local} = \Delta x_a \cdot r_x/a,\qquad
\Delta y_{local} = \Delta y_a \cdot r_y/a
\]

代码位置：

- [`src/dynamics/adv/adv_batch_mod.F90:594`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_batch_mod.F90#L594)
- [`src/dynamics/adv/adv_batch_mod.F90:713`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_batch_mod.F90#L713)

这对应“每层水平格距随高度增大”的 deep 几何事实。

### 11.3 tracer CFL 与整格穿越质量改用局地面元/边长

`adv_batch_calc_cflxy_tracer` 现在统一使用：

- `adv_batch_face_length_x/y`
- `adv_batch_area`

其中：

\[
L_{face} \propto r/a,\qquad A_{cell} \propto (r/a)^2
\]

代码位置：

- [`src/dynamics/adv/adv_batch_mod.F90:748`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_batch_mod.F90#L748)

这保证了 semi-Lagrangian/FFSL 在计算“一次跨过几个整格、剩余分数是多少”时，用的是局地壳面几何，而不是参考球几何。

### 11.4 Upwind/FFSL 的跨整格输送改成局地壳面面积

`upwind_mod` 和 `ffsl_mod` 里，原来凡是形如：

\[
\sum m \cdot \Delta x / \Delta t
\]

的整格穿越部分，现在都改成了：

\[
\sum m \cdot A_{local} / L_{face,local} / \Delta t
\]

代码位置：

- [`src/dynamics/adv/upwind_mod.F90:64`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/upwind_mod.F90#L64)
- [`src/dynamics/adv/ffsl_mod.F90:420`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/ffsl_mod.F90#L420)

所以现在不仅 CFL 本身是 deep 的，真正参与 full-cell transport 的质量和示踪物累积量也已经是 deep 的。

### 11.5 水平散度改成 deep directional divergence

除了 `div_operator(fx, rx, fy, ry, div)` 这类整散度，FFSL 分裂里还需要 x/y 方向单独散度。

因此新增了：

- `divx_operator_deep`
- `divy_operator_deep`

代码位置：

- [`src/meshes/latlon/latlon_operators_mod.F90:196`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/meshes/latlon/latlon_operators_mod.F90#L196)

对应调用点已经接入：

- `adv_batch_prepare_semilag_h`
- `ffsl_calc_mass_hflx_swift`
- `swift_prepare`
- `ffsl_calc_tracer_hflx_swift`

代码位置：

- [`src/dynamics/adv/adv_batch_mod.F90:1012`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_batch_mod.F90#L1012)
- [`src/dynamics/adv/ffsl_mod.F90:168`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/ffsl_mod.F90#L168)
- [`src/dynamics/adv/ffsl_mod.F90:240`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/ffsl_mod.F90#L240)
- [`src/dynamics/adv/ffsl_mod.F90:317`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/ffsl_mod.F90#L317)

### 11.6 `mass` 与 `tracer` 的主更新都已接入 deep 水平散度

`adv_run_mass` 和 `adv_run_tracers` 在 deep 路径下都已调用带 `rdp_lon/rdp_lat` 的散度。

代码位置：

- [`src/dynamics/adv/adv_mod.F90:263`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_mod.F90#L263)
- [`src/dynamics/adv/adv_mod.F90:343`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/adv/adv_mod.F90#L343)

到这里，GMCORE 里的水平平流链条已经从“只有一部分算子 deep 化”，变成了：

- CFL 用局地长度
- full-cell transport 用局地面积和边长
- 主更新散度用 deep metric
- FFSL 分裂方向散度也用 deep metric

## 12. 与 ICON deepatmo 实现的对照判断

你给的 ICON 代码里，最关键的信息不是某一个公式本身，而是它把 deep 几何 metric 明确分成了几类，而不是只用一个统一的 `r/a`。

### 12.1 ICON 把 gradient、divergence、volume 分开定义

在 [`mo_vertical_grid.f90:2438`](/Users/cuiqitao/Desktop/深大气/垂直分层设计/icon_deepatmo/mo_vertical_grid.f90#L2438) 开始，ICON 定义了：

- `deepatmo_gradh_mc`
- `deepatmo_divh_mc`
- `deepatmo_divzL_mc`
- `deepatmo_divzU_mc`
- `deepatmo_vol_mc`

这意味着 ICON 的基本立场是：

- 水平梯度 metric 不等于水平散度 metric
- 垂直下界面散度 metric 不等于垂直上界面散度 metric
- 体积 metric 也单独保存

所以从几何精细度看，ICON 明显比“一个 `rdp` 场外加插值”的做法更完整。

### 12.2 ICON 的水平质量散度是在 divergence 端乘 metric

在 [`mo_solve_nonhydro.f90:2130`](/Users/cuiqitao/Desktop/深大气/垂直分层设计/icon_deepatmo/mo_solve_nonhydro.f90#L2130)，ICON 的水平质量散度写成：

\[
\mathrm{div}_h(\mathrm{mass\_flux}) \times deepatmo\_divh\_mc
\]

也就是说，它更接近：

- 通量保持物理通量定义
- metric 由散度端统一消费

这和 GMCORE 这次把 `calc_mf` 里的预乘 `a/r` 去掉，是同一个方向。

### 12.3 ICON 的速度平流项也区分不同 metric

在 [`mo_velocity_advection.f90:379`](/Users/cuiqitao/Desktop/深大气/垂直分层设计/icon_deepatmo/mo_velocity_advection.f90#L379) 以及 [`mo_velocity_advection.f90:750`](/Users/cuiqitao/Desktop/深大气/垂直分层设计/icon_deepatmo/mo_velocity_advection.f90#L750)，ICON 对水平梯度类项使用 `deepatmo_gradh_*`，而不是 `deepatmo_divh_*`。

这进一步说明：

- `a/r` 风格的梯度修正
- 水平通量散度修正

在严格 deep 几何下并不应强行合并成同一个标量因子。

### 12.4 对 GMCORE 当前实现的判断

参照 ICON，我对 GMCORE 当前实现的判断是：

1. 方向是对的。
   现在 GMCORE 已经把水平平流从 shallow 几何推到“局地长度 + 局地面积 + deep divergence”的一致框架里了。

2. 复杂度还低于 ICON。
   GMCORE 当前仍主要依赖 `rdp/rdp_lon/rdp_lat/rdp_lev` 这一族半径场来近似所有 deep metric，没有像 ICON 那样显式拆成 `gradh/divh/divzL/divzU/vol`。

3. 因而它更像一个“几何一致的一阶闭环”，还不是 ICON 那种“按控制体精确推导出来的全 metric 池”。

## 13. `rdp` 先数值计算、再经 `interp_mod` 插值，会不会掉精度

会有一定误差，但我判断它现在**不是主导误差源**。

### 13.1 为什么一定会有误差

`calc_rdp` 先在格心通过质量和密度积分构造 `rdp`：

- [`src/dynamics/operators_mod.F90:1235`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L1235)

然后再通过：

- `average_run`
- `interp_run`

得到 `rdp_lon`、`rdp_lat`、`rdp_lev` 等场：

- [`src/dynamics/operators_mod.F90:1267`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/operators_mod.F90#L1267)
- [`src/dynamics/interp_mod.F90:71`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/dynamics/interp_mod.F90#L71)

其中水平边插值很多就是简单平均，竖直插值则是按层厚加权的线性插值。  
所以从“几何真值”角度说，它当然不是严格精确半径。

### 13.2 但它带来的误差通常小于“算子前后不一致”

在这次修改之前，GMCORE 更大的问题其实是：

- CFL 用参考球长度
- full-cell transport 用参考球面积
- 散度又在别处部分乘了 `a/r`

这种误差不是插值误差，而是**算子不一致误差**。  
和这种误差相比，`rdp` 做一次规则插值带来的截断误差通常是次一级的。

### 13.3 真正需要警惕的不是“插值本身”，而是“一个 `rdp` 包打天下”

参照 ICON，真正的潜在问题是：

- `gradh` 和 `divh` 并不完全相同
- `divzL` 和 `divzU` 也不完全相同
- `vol` 也有独立定义

如果 GMCORE 长期只用一个 `rdp*` 族去近似所有这些 metric，那么随着：

- 垂直层变厚
- 模式顶更高
- 深大气效应更强

它和严格几何推导之间的偏差会逐步显现。

### 13.4 当前阶段的工程判断

如果目标是“先把 GMCORE 从 shallow 几何改成自洽 deep 几何”，那现在这种做法是合理的，值得保留。  
如果目标是“尽量逼近 ICON 那套严格 deepatmo metric”，那后续应考虑进一步拆分出：

- 水平梯度 metric
- 水平散度 metric
- 垂直上下界面散度 metric
- 体积 metric

而不是继续只依赖 `rdp` 插值。

## 14. 当前结论

对你这次关心的两件事，我的结论是：

1. GMCORE 的 `mass`、`tracer`、CFL、FFSL/Upwind full-cell transport、分裂方向散度，这次已经基本补成 deep 水平平流版本了。
2. `rdp + interp_mod` 确实会引入一些近似，但它现在更像“合理的工程简化”，不是最危险的地方。
3. 真正和 ICON 还有差距的地方，在于 ICON 把 deep metric 明确拆成了 `gradh/divh/divzL/divzU/vol`，而 GMCORE 目前还是单一 `rdp*` 家族近似。
