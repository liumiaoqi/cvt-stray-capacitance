(*  stray_capacitance.wls
    电容式电压互感器杂散电容误差的定量复算
    对应课程报告《电容式电压互感器中电场分布与杂散电容误差分析》

    运行：  wolframscript -file stray_capacitance.wls
    需要：  Wolfram Engine 12.0+（开发者许可免费）

    脚本做四件事：
      1. 用报告图 4 的参数复算七个数据点，逐个对照报告取值
      2. 验证高压臂/低压臂的灵敏度之比等于 K-1
      3. 给出 10 pF 并入两臂各自的比值差，量化不对称
      4. 核 GB/T 20840.5 表 501 的准确级限值，并算出 1.0% 相当于几倍

    任一断言超差即退出码 1。
*)

pass = 0; fail = 0;

check[label_, got_, want_, tol_, unit_:""] := Module[{ok, shown},
  ok = Abs[got - want] <= tol;
  If[ok, pass++, fail++];
  shown = ToString[NumberForm[N[got, 8], {10, 6}]];
  Print[If[ok, "  OK   ", "  FAIL "],
        StringPadRight[label, 34], " got ", StringPadRight[shown, 12],
        " expect ", want, " +/- ", tol, " ", unit]
];

(* ---------- 报告图 4 所用参数 ---------- *)
c1 = 1.*10^-9;      (* 高压臂 1 nF *)
c2 = 1.*10^-6;      (* 低压臂 1 uF *)
dcMax = 10.*10^-12; (* 杂散电容峰值 10 pF *)
K[a_, b_] := (a + b)/a;
kn = K[c1, c2];     (* 额定变比取设计变比 *)

Print["[params]  C1 = 1 nF, C2 = 1 uF, Cspur = 10 pF*cos(theta)"];
Print["[params]  K = Kn = ", kn];
Print[""];

(* ---------- 1. 图 4 七个数据点 ---------- *)
Print["[1] loss ratio vs lead angle"];
kprime[thDeg_] := K[c1 + dcMax*Cos[thDeg*Degree], c2];
eps[thDeg_] := (kn/kprime[thDeg] - 1)*100;

reported = {{0, 0.99899}, {15, 0.96495}, {30, 0.86515}, {45, 0.70640},
            {60, 0.49950}, {75, 0.25856}, {90, 0.00000}};

Do[
  th = reported[[i, 1]];
  check["eps(" <> ToString[th] <> " deg), %", eps[th], reported[[i, 2]], 1.*10^-5, "%"],
  {i, Length[reported]}
];

(* ---------- 2. 灵敏度之比 ---------- *)
Print[""];
Print["[2] sensitivity asymmetry"];
d1 = Abs[D[K[a, b], a] /. {a -> c1, b -> c2}];
d2 = Abs[D[K[a, b], b] /. {a -> c1, b -> c2}];
check["|dK/dC1| / |dK/dC2|", d1/d2, kn - 1, 1.*10^-6, ""];

(* ---------- 3. 10 pF 并入两臂的对比 ---------- *)
Print[""];
Print["[3] same 10 pF on either arm"];
eHigh = (kn/K[c1 + dcMax, c2] - 1)*100;
eLow  = (kn/K[c1, c2 + dcMax] - 1)*100;
check["10 pF on HV arm, %", eHigh, 0.99899, 1.*10^-4, "%"];
check["10 pF on LV arm, %", eLow, -0.00100, 1.*10^-4, "%"];
Print["         ratio of effects = ", Abs[eHigh/eLow]];

(* ---------- 4. 准确级限值与倍数 ---------- *)
(* 限值引用自 GB/T 20840.5-2013 表 501（测量用 CVT）：
   0.2 级 ±0.2% / ±10'；0.5 级 ±0.5% / ±20'；1.0 级 ±1.0% / ±40'。
   这三组是标准给定值，脚本只核对「1.0% 相当于几倍」，不重新推导限值本身。 *)
Print[""];
Print["[4] limit multiple (limits per GB/T 20840.5-2013 Table 501)"];
check["max eps / 0.2 class limit", 0.99899/0.2, 5.0, 0.01, "x"];
check["max eps / 0.5 class limit", 0.99899/0.5, 2.0, 0.01, "x"];

Print[""];
Print["---- ", pass, " passed, ", fail, " failed ----"];
If[fail > 0, Exit[1]];
