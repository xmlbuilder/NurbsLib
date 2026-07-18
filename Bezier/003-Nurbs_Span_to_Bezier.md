# NURBS Span to Bezier Guide

## 1. 개요

- NURBS 곡선은 여러 개의 **Non-zero Knot Span**으로 구성된다.

- 각 Span은 정확히 하나의 Bezier Curve로 변환될 수 있다.

- 즉,

$$
\boxed{
\text{One Non-zero Knot Span}
\Longleftrightarrow
\text{One Bezier Curve}
}
$$

---

## 2. NURBS 정의

- 차수가 \(p\)인 NURBS(B-spline) 곡선은

$$
C(u) =
\sum_{i=0}^{n}
N_{i,p}(u)P_i
$$

- 로 정의된다.

- 여기서

- $N_{i,p}(u)$ : B-spline Basis
- $P_i$ : Control Point

이다.

---

## 3. Non-zero Knot Span

- Knot Vector가

$$
U=\{u_0,u_1,\cdots,u_m\}
$$

- 일 때 실제 Span은

$$
u_{i+1}>u_i
$$

- 인 구간뿐이다.

- 따라서 Span 개수는

$$
\boxed{
SpanCount=
\sum_i
\mathbf 1
\left(
u_{i+1}>u_i
\right)
}
$$

이다.

- 예를 들어

```text
0 0 0 0 1 2 3 4 4 4 4
```

- 에서는

```text
[0,1]
[1,2]
[2,3]
[3,4]
```

- 의 네 개 Span만 존재한다.

---


## 4. Span과 Bezier의 관계

- Span 하나는 하나의 다항식 조각이다.

- 경계 Knot의 Multiplicity가

$$
p+1
$$

이 되면

- 해당 Span은 독립적인 Bernstein 다항식이 된다.

- 즉

$$
\boxed{
\text{Clamped B-spline Span} =
\text{Bezier Curve}
}
$$

이다.

---

## 5. Parameter 변환

- Bezier는 항상

$$
0\le t\le1
$$

에서 정의된다.

- 원래 NURBS Parameter는

$$
u_a\le u\le u_b
$$

이다.

- 두 Parameter의 관계는

$$
u =
u_a+t(u_b-u_a)
$$

이며

- 역변환은

$$
t=
\frac{u-u_a}
{u_b-u_a}
$$

이다.

---

## 6. API 사용 순서

```cpp
const int span_count = curve.NonZeroSpanCount();

for (int i = 0; i < span_count; ++i)
{
    BezierCurve bezier;
    Interval domain;

    curve.ConvertNonZeroSpanToBezier(
        i,
        bezier,
        &domain);

    // domain = 원래 NURBS Parameter
}
```

---

## 7. Span 번호 예

- Degree=3

```text
- Knots

0 0 0 0 1 2 3 4 4 4 4
```

|Bezier Span|Original Parameter|
|-----------:|------------------|
|0|[0,1]|
|1|[1,2]|
|2|[2,3]|
|3|[3,4]|

- 예를 들어

```cpp
curve.ConvertNonZeroSpanToBezier(1,...)
```

은

$$
u\in[1,2]
$$

- 구간을 Bezier Curve로 반환한다.

---

## 8. Bezier 평가

- Bezier는

$$
B(t) =
\sum_{i=0}^{p}
B_i^p(t)Q_i
$$

로 계산한다.

여기서

- $Q_i$ : Bezier Control Point
- $B_i^p$ : Bernstein Basis

이다.

---

## 9. 왜 Span 단위로 사용하는가

- Bezier Patch는

- Adaptive Tessellation
- GPU Evaluation
- Subdivision
- Trim 계산

에 매우 적합하다.

- 그래서 대부분의 CAD Kernel은

```text
NURBS
    ↓
Bezier Patch
    ↓
Rendering / Tessellation
```

순서로 처리한다.

---

## 10. 전체 과정

```text
NURBS Curve
      │
      ▼
Non-zero Span Count
      │
      ▼
Span 선택
      │
      ▼
ConvertNonZeroSpanToBezier()
      │
      ▼
Bezier Curve
      │
      ▼
Bezier Evaluation
```

---

## 11. 핵심 정리

- NURBS에서 Bezier Curve를 얻는 순서는

$$
\boxed{
\text{Span Count}
\rightarrow
\text{Span 선택}
\rightarrow
\text{Bezier 추출}
}
$$

이다.

- Bezier Curve는 NURBS 전체에서 하나만 존재하는 것이 아니라,

**각 Non-zero Knot Span마다 정확히 하나씩 존재한다.**

---

