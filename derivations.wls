(*  derivations.wls
    把报告里借自文献的定量关系重新推导一遍，看是否成立

    报告引用的文献给出了若干结论，但结论本身要能自己算出来才算数。
    本脚本复现其中可解析验证的四条：

      A. 温度：两臂温度系数不同 → 分压比偏移（报告只作定性描述）
      B. 频率：两臂损耗角不同 → 比值差随频率变化（报告只作定性描述）
      C. 能量法求电容：验证「W = CV²/2 → 从场解提取电容」这一步
         （liugang2020 用有限元算空间电场能量作为方程右端项，此为其基础）
      D. 感应系数矩阵与电容矩阵互逆：[C] = [α]^(-1)

    运行：wolframscript -file derivations.wls
    任一断言超差即退出码 1。
*)

pass = 0; fail = 0;
check[label_, got_, want_, tol_] := Module[{ok, shown},
  ok = Abs[got - want] <= tol;
  If[ok, pass++, fail++];
  shown = ToString[NumberForm[N[got, 8], {12, 8}]];
  Print[If[ok, "  OK   ", "  FAIL "], StringPadRight[label, 40],
        " got ", StringPadRight[shown, 14], " expect ", want, " +/- ", tol]
];

c10 = 1.*10^-9; c20 = 1.*10^-6;
K[a_, b_] := (a + b)/a;
k0 = K[c10, c20];
Print["[setup]  C1 = 1 nF, C2 = 1 uF, K = ", k0];
Print[""];

(* ================= A. 温度：两臂温度系数不同 ================= *)
Print["[A] temperature: unequal coefficients on the two arms"];
(*  C1(T) = C10 (1 + b1 dT),  C2(T) = C20 (1 + b2 dT)
    K'(dT) = 1 + (C20/C10) (1 + b2 dT)/(1 + b1 dT)
    一阶近似： dK/K = ((K-1)/K) * (b2 - b1) * dT                                  *)
b1 = 2.*10^-4; b2 = 1.*10^-4; dT = 50.;   (* 典型介质温度系数量级，非实测 *)
dKKexact = (K[c10*(1 + b1*dT), c20*(1 + b2*dT)] - k0)/k0;
dKKapprox = ((k0 - 1)/k0)*(b2 - b1)*dT;
check["dK/K exact vs first-order", dKKexact, dKKapprox, Abs[dKKapprox]*0.02];
Print["         dK/K = ", N[dKKexact*100, 6], " %  (即比值差 ", N[-dKKexact*100, 4], " %)"];

(* 若两臂介质相同（b1 = b2），分压比不变 —— 这是报告结论的核心 *)
check["same medium -> no drift", (K[c10*(1 + b1*dT), c20*(1 + b1*dT)] - k0)/k0, 0., 1.*10^-15];

(* ================= B. 频率：分压器本身不敏感，失谐发生在电抗器 ================= *)
Print[""];
Print["[B] frequency: the capacitive divider itself is frequency-flat"];
(*  并联模型 Y = j w C + 1/R = w C (tan d + j)，串联模型 Z = (tan d - j)/(w C)，
    两种模型下 w 都是公因子，故
        U2/U1 = C1(td1+j) / [C1(td1+j) + C2(td2+j)]
    与频率无关。即使两臂损耗角不同，结论也一样。                              *)
ratio[w_, td1_, td2_] := Module[{y1, y2},
  y1 = w*c10*(td1 + I); y2 = w*c20*(td2 + I);
  y1/(y1 + y2)
];
r50 = ratio[2.*Pi*50., 0.001, 0.0005];
r502 = ratio[2.*Pi*50.2, 0.001, 0.0005];
check["|ratio(50Hz) - ratio(50.2Hz)|", Abs[r50 - r502], 0., 1.*10^-15];
Print["         ratio magnitude (any freq) = ", N[Abs[r50], 8]];
Print["         -> 分压器对频率不敏感，报告 4 章「损耗角不同则频率影响」这一说法不成立"];

(* 真正的频率敏感性来自补偿电抗器：额定频率下 w0 L = 1/(w0 Cs)，
    偏离 delta 时残抗约为 X0 的 2 delta *)
Print[""];
Print["     compensation reactor detuning"];
delta = 0.004;   (* 50 Hz -> 50.2 Hz *)
residExact = ((1 + delta) - 1/(1 + delta));
check["residual reactance / X0 (exact)", residExact, 2*delta, 0.0001];
Print["         delta = ", delta, "  -> 残抗 = ", N[residExact*100, 4], " % of X0"];

(* ================= C. 能量法求电容 ================= *)
Print[""];
Print["[C] capacitance from field energy (coaxial check)"];
(*  同轴：E(r) = U/(r ln(b/a))，W = 1/2 eps Integrate[E^2 dV]，C = 2W/U^2
    解析值 C = 2 pi eps/ln(b/a)。用精确有理数做积分才能求出闭式。              *)
aR = 1/100; bR = 1/40;
eField[r_] := 1/(r*Log[bR/aR]);
wEnergy = (1/2)*eps0*Integrate[eField[r]^2*2*Pi*r, {r, aR, bR},
                              Assumptions -> eps0 > 0];
cEnergy = 2*wEnergy;                       (* U = 1，单位长度 *)
cAnalytic = 2*Pi*eps0/Log[bR/aR];
check["C(energy) / C(analytic)", Simplify[cEnergy/cAnalytic], 1, 1.*10^-12];
Print["         b/a = 2.5, C/eps0 = ", N[cAnalytic/eps0, 8]];

(* ================= D. 感应系数矩阵与电容矩阵互逆 ================= *)
Print[""];
Print["[D] potential-coefficient matrix and capacitance matrix are inverses"];
(*  多导体：Q = C V，V = alpha Q  =>  C = alpha^(-1)
    取一个三导体对称例子验证                                              *)
alpha = {{1.20, 0.35, 0.18}, {0.35, 1.20, 0.35}, {0.18, 0.35, 1.20}}*10^9;
cMat = Inverse[alpha];
check["alpha . C = I (off-diagonal sum)", Total[Abs[Flatten[alpha.cMat - IdentityMatrix[3]]]], 0., 1.*10^-12];
(* 自电容 = 该行全部元素之和；互电容 = 非对角元取负 *)
cSelf = Total[cMat[[1]]];
cMutual = -cMat[[1, 2]];
Print["         C_self(1) = ", N[cSelf*10^12, 6], " pF"];
Print["         C_mutual(1,2) = ", N[cMutual*10^12, 6], " pF"];
check["row sum consistency (Q for V={1,0,0})", cMat[[1, 1]] + cMat[[1, 2]] + cMat[[1, 3]], cSelf, 1.*10^-18];

(* ================= E. 概率：超差概率 ================= *)
Print[""];
Print["[E] probability of exceeding class 0.2"];
(*  设比值差 eps ~ N(mu, sigma^2)，限值 a。
    mu = 0 时 P = 2 Phi(-a/sigma)。
    注意：sigma 未实测，本段只复算报告表 1 的数值，不主张某一 sigma 成立。     *)
a = 0.2;
prob[a_, s_] := 2*CDF[NormalDistribution[0, 1], -a/s];
check["sigma=0.05 -> P", prob[0.2, 0.05], 6.3*10^-5, 1.*10^-6];
check["sigma=0.10 -> P", prob[0.2, 0.10], 0.046, 5.*10^-4];
check["sigma=0.20 -> P", prob[0.2, 0.20], 0.317, 5.*10^-4];
(* mu 已在校正前处于 1% 量级时，超差近乎必然 *)
pMu = 1 - (CDF[NormalDistribution[0, 1], (a - 1.0)/0.1] -
           CDF[NormalDistribution[0, 1], (-a - 1.0)/0.1]);
check["mu=1.0%, sigma=0.1% -> P", pMu, 1.0, 1.*10^-9];
Print["         -> 校正前超差是必然事件，故「校正先于控制离散」"];

(* ================= F. 切比雪夫上界（报告表 2） ================= *)
Print[""];
Print["[F] Chebyshev bound, table 2"];
(*  P(|eps - mu| >= k sigma) <= 1/k^2，k = a/sigma。
    这一条是定理，不是待验证的猜想；此处核的是表 2 每一格与该式一致，
    以及反向解出的 sigma 上限。                                              *)
cheb[a_, s_] := (s/a)^2;      (* 即 1/k^2 *)
check["sigma=0.05 -> 1/k^2", cheb[0.2, 0.05], 0.0625, 1.*10^-12];
check["sigma=0.10 -> 1/k^2", cheb[0.2, 0.10], 0.2500, 1.*10^-12];
check["sigma=0.20 -> 1/k^2", cheb[0.2, 0.20], 1.0000, 1.*10^-12];

(* 反向：要保证 P <= 5%，需 k >= sqrt(20)，故 sigma <= a/sqrt(20) *)
sigMax = 0.2/Sqrt[20.];
check["sigma max for P<=5% (%)", sigMax, 0.0447214, 1.*10^-6];
check["sigma max / a", sigMax/0.2, 1/Sqrt[20.], 1.*10^-12];
Print["         正态在该 sigma 处给 P = ", N[prob[0.2, sigMax], 8],
      "  （切比雪夫只保证 <= 5%，相差约 ",
      N[0.05/prob[0.2, sigMax], 4], " 倍）"];

Print[""];
Print["---- ", pass, " passed, ", fail, " failed ----"];
If[fail > 0, Exit[1]];
