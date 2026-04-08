# Gerkema 图11 下 `N = 2Ω` 的临界纬度传播算例

这个算例的目标是做一个最小对比：

- `shallow / TA-like`：`deepwater = .false.`
- `deep / 含 2Ωcosφ`：`deepwater = .true.`，并让 `use_hor_nct = .true.`

我们关心的是同一个惯性重力内波包向高纬传播时，是否会在传统惯性纬度附近被阻挡，以及在包含非传统项后是否还能继续向极传播。

## 1. 选取的基准参数

在 [`critical_wave_test_mod.F90`](/Users/cuiqitao/work_dir/gmcore_deep_control/src/tests/critical_wave_test_mod.F90) 里，默认值已经改成下面这一组：

- `N = 2Ω`
- `omega_wave = 1.5Ω`
- `lat0 = 20 deg`
- `pert_type = 2`

这样选的原因是：

- 频率必须满足 `f(lat0) < ω < N`，这样源区能支持 IGW。
- `ω = 1.5Ω` 时，传统与非传统两条临界纬度之间有大约 `15 deg` 的分离，信号足够明显。
- `lat0 = 20 deg` 离两条临界纬度都比较远，方便波包先形成再向高纬传播。

## 2. 理论预期

传统近似下，转折发生在惯性纬度：

`sin(phi_i) = omega_wave / (2Ω)`

代入 `omega_wave = 1.5Ω` 得到：

- `phi_i ≈ 48.59 deg`

Gerkema 文中图 11 所对应的全科氏力边界，可写成：

`sin^2(phi_c) = omega_wave^2 * (N^2 + 4Ω^2 - omega_wave^2) / (4Ω^2 N^2)`

代入 `N = 2Ω, omega_wave = 1.5Ω` 得到：

- `phi_c ≈ 64.06 deg`

所以这个算例最核心的检验量就是：

- shallow 情况下，波包主能量应在 `49 deg` 左右附近衰减或反射
- deep 情况下，波包若保持可辨识，应能继续推进到 `64 deg` 左右附近

## 3. 当前更推荐的激发方式

现在更推荐把 `crt` 当成“波源型”算例，而不是“直接塞一个解析波包”的初值算例。

默认行为是：

- 初值只保留平静背景
- 用一个平滑、短时、定频的局地热源去辐射 IGW
- 热源默认在纬向和垂直方向高斯局地化
- 默认 `forcing_zonal_uniform = 1`，即源在经向上是均匀带状的，专门用来激发南北传播，避免全球经度拼接处出现“断开的双包”

这样做的原因很直接：

- 在 `N = 2Ω` 这种弱层结设定下，自由 IGW 的水平与垂直尺度往往是同量级的。
- 对 30 层、1° 全球格点来说，想同时构造“解析极化正确、尺度又充分分辨、还不受经度拼接影响”的初值波包，其实很别扭。
- 用短时平滑波源更容易得到干净的辐射波列，也更接近“源激发后向临界纬度传播”的物理图景。

## 4. 如果你仍然想用旧的初值波包

旧的解析初值波包还在，但现在不是默认值：

- `pert_type = 0`：默认，不加初值波包
- `pert_type = 1`：纯 `theta` 波包
- `pert_type = 2`：解析极化的 `u/v/w/pt` 波包

如果你发现初值一开始就有明显条纹、拼接或断续，优先建议回到 `pert_type = 0`，只保留强迫源。

## 5. 当前测试模块做了什么

`crt` 测试现在会：

- 构造常 `N` 的静止干大气背景
- 如果 `pert_type > 0`，可以额外加入旧版初值波包
- 如果 `forcing_mode = 1`，则在积分过程中用平滑热源持续 `forcing_ncycles` 个周期

## 6. 推荐的 namelist 方向

建议先做两组完全相同的积分，只改深大气开关：

```fortran
&gmcore_control
  test_case      = 'crt'
  deepwater      = .false.
/

&critical_wave_control
  pert_type      = 0
  forcing_mode   = 1
  forcing_zonal_uniform = 1
  forcing_ncycles = 3.0d0
  forcing_pt_amp = 2.0d-4
  N_ratio        = 2.0d0
  omega_wave_ratio = 1.5d0
  lat0_deg       = 20.0d0
/
```

注意：

- 不要在 Fortran namelist 里写 `2.0d0 * 7.292115d-5` 这种表达式，很多情况下它不会被当成乘法表达式求值。
- 如果你更喜欢直接填 SI 数值，就写成：
  `N_freq = 1.4584230d-4`
  `omega_wave = 1.09381725d-4`
  `lat0 = 3.4906585d-1`

第二组只改成：

```fortran
&gmcore_control
  test_case      = 'crt'
  deepwater      = .true.
  use_hor_nct    = .true.
/
```

如果你希望更干净地只看 `2Ωcosφ` 的影响，可以保留 `deepwater = .true.`，同时明确检查其他 deep 子开关是否也一起打开了。

## 7. 最推荐看的诊断量

- 固定经度上 `|w|` 或 `|v|` 的纬高剖面随时间演变
- 波包包络最大值对应的纬度 `phi_peak(t)`
- 纬向积分后的能量指标，例如 `E ~ u'^2 + v'^2 + w'^2 + (g/N)^2 (theta'/theta)^2`

最直接的图就是：

- 同一时刻 shallow/deep 两张 `latitude-height` 剖面
- 在图上叠加 `48.6 deg` 与 `64.1 deg` 两条竖线

## 8. 如果数值上没有看到明显差异

优先按下面顺序调：

- 降低 `forcing_pt_amp`，先保证线性传播
- 增大 `forcing_ncycles`，让频谱更窄
- 把 `lat0_deg` 再往低纬放一点，比如 `15 deg`
- 如果想让波源更局地，就减小 `Ly`；如果想让辐射更平滑，就适当增大 `Ly`
- 如果要回到旧波包模式，再单独调 `pert_type`

## 9. 这一步之后最自然的下一步

如果 `N = 2Ω` 这组能看到 clear signal，下一步就很自然：

- 保持 `ω / Ω` 不变，扫描 `N / 2Ω = 0.5, 1, 2`
- 或保持 `N = 2Ω`，扫描 `ω / Ω = 1.2, 1.4, 1.5, 1.6`

这样可以直接对应回 Gerkema 图 11 的“允许频带随纬度变化”。
