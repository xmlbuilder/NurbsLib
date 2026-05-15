# Power Reparam Matrix

🧩 두 함수의 핵심 차이 (한 줄 요약)

| Function | Purpose | Reparam Formula | Meaning |
|---|---|---|---|
| on_reparameterization_matrix_unit | Unit-domain reparameterization | u = a + (b-a)t | old polynomial f(u)를 new parameter t ∈ [0,1] 기준으로 재전개 |
| on_reparameterization_matrix | General-domain reparameterization | u = α·up + β | old polynomial f(u)를 new parameter up ∈ [a,b] 기준으로 재전개 |

즉:

- unit 버전은:

```math
f(u)\rightarrow f(a+(b-a)t)
```

- 일반 버전은:

```math
f(u)\rightarrow f(\alpha up+\beta)
```

---

## 🧠 1) on_reparameterization_matrix_unit 
- unit-domain reparameterization

## ✔ 치환식

```math
u = a + (b-a)t
```

여기서:

```math
t \in [0,1]
```

이고, `u`는 기존 polynomial의 변수이다.

즉 기존 함수가:

```math
f(u)=\sum_{k=0}^{p} c_k u^k
```

일 때 새 함수는:

```math
g(t)=f(a+(b-a)t)
```

이다.

따라서:

```math
g(t)=\sum_{i=0}^{p} d_i t^i
```

이고:

```math
d_i=\sum_{k=i}^{p} c_k {k \choose i} a^{k-i}(b-a)^i
```

이다.

---

## ✔ 행렬 형태

```math
R[i,k]={k \choose i}a^{k-i}(b-a)^i
```

단:

```math
i \le k
```

일 때만 값이 있고 나머지는 0이다.

즉 upper-right triangular matrix 구조를 가진다.

---

## ✔ 의미

이 함수는:

```text
old polynomial f(u)
```

를:

```text
new local parameter t ∈ [0,1]
```

기준 polynomial로 다시 표현한다.

즉:

```text
u-domain의 [a,b] 구간을
local parameter [0,1] 위에서 재표현
```

하는 것이다.

---

## ✔ 대표 사용처

- Bézier subdivision
- Bézier segment extraction
- curve segment `[a,b]`를 local `[0,1]` polynomial로 재표현
- trimming curve 일부 구간을 local parameter로 재정의
- power basis polynomial의 local segment 생성
- Bernstein subdivision 내부 처리

---

## 🧠 2) on_reparameterization_matrix 
- 일반 domain reparameterization

## ✔ 치환식

```math
u=\alpha up+\beta
```

조건:

```math
up=a \Rightarrow u=u0
```

```math
up=b \Rightarrow u=u1
```

따라서:

```math
\alpha=\frac{u1-u0}{b-a}
```

```math
\beta=u0-\alpha a
```

---

## ✔ 기존 함수

기존 polynomial:

```math
f(u)=\sum_{i=0}^{p} A_i u^i
```

를 새 변수 `up` 기준으로:

```math
g(up)=f(\alpha up+\beta)
```

로 바꾼다.

즉:

```math
g(up)=\sum_{j=0}^{p} A'_j up^j
```

이며:

```math
A'_j=\sum_{i=j}^{p}A_i{i \choose j}\beta^{i-j}\alpha^j
```

이다.

---

## ✔ 행렬 형태

```math
M[j,i]={i \choose j}\beta^{i-j}\alpha^j
```

단:

```math
j \le i
```

일 때만 값이 있고 나머지는 0이다.

즉 역시 upper-right triangular matrix 구조를 가진다.

---

## ✔ 대표 사용처

- 임의 domain remapping
- surface reparameterization
- trimming domain 변경
- NURBS patch normalization
- arbitrary parameter interval remapping
- power basis surface domain 변경

---

# 🎯 결정적 차이

`on_reparameterization_matrix_unit(p, a, b)` 는 사실상 다음 일반 함수와 동일하다.

```rust
on_reparameterization_matrix(p, a, b, 0.0, 1.0)
```

왜냐하면:

```math
up \in [0,1]
```

일 때:

```math
up=0 \Rightarrow u=a
```

```math
up=1 \Rightarrow u=b
```

이므로:

```math
u=a+(b-a)up
```

가 되기 때문이다.

즉 unit 버전은 일반 버전의 특수한 경우이다.

---

# 🧱 언제 어떤 함수를 써야 하나?

## ✔ on_reparameterization_matrix_unit

이럴 때 사용:

```text
원래 polynomial의 [a,b] 구간을
새 local parameter [0,1] 위에서 표현하고 싶다.
```

예:

```text
Bezier curve의 u=[0.25,0.75] 구간을
새 Bezier [0,1]로 만들기
```

---

## ✔ on_reparameterization_matrix

이럴 때 사용:

```text
새 parameter domain이 [0,1]이 아니라
임의의 [a,b]일 때
```

예:

```text
old u-domain [u0,u1]을
new up-domain [a,b] 기준으로 표현
```

---

## 🎯 한 줄 요약

- unit 버전:

```math
f(u)\rightarrow f(a+(b-a)t),\quad t\in[0,1]
```

- 일반 버전:

```math
f(u)\rightarrow f(\alpha up+\beta),\quad up\in[a,b]
```

- 즉 둘 다:

```text
polynomial variable substitution matrix
```

이지만,

- unit 버전은 local `[0,1]`
- 일반 버전은 arbitrary domain

을 대상으로 한다.

---

```cpp
bool basis::on_reparameterization_matrix_unit(
    const int p,
    const Real a,
    const Real b,
    Matrix& r)
{
    if (p < 0)
        return false;

    constexpr Real param_tol = 1.0e-12;

    const Real s = b - a;

    if (std::abs(s) <= param_tol)
        return false;

    const int n = p + 1;

    if (!r.Create(n, n))
        return false;

    (void)r.Fill(0.0);

    for (int i = 0; i <= p; ++i)
    {
        const Real si = std::pow(s, i);

        for (int k = i; k <= p; ++k)
        {
            const Real cki = on_binomial_real(k, i);
            const Real apow = std::pow(a, k - i);

            r(i, k) = cki * apow * si;
        }
    }

    return true;
}
```

```cpp
bool basis::on_reparameterization_matrix(
    const int p,
    const Real a,
    const Real b,
    const Real ap,
    const Real bp,
    Matrix& rem)
{
    if (p < 0)
        return false;

    constexpr Real param_tol = 1.0e-12;

    if (std::abs(b - a) <= param_tol)
        return false;

    const int n = p + 1;

    if (!rem.Create(n, n))
        return false;

    (void)rem.Fill(0.0);

    const Real c = (bp - ap) / (b - a);
    const Real d = (ap * b - bp * a) / (b - a);

    rem(0, 0) = 1.0;

    for (int i = 1; i <= p; ++i)
    {
        rem(0, i) = d * rem(0, i - 1);
        rem(i, i) = c * rem(i - 1, i - 1);
    }

    for (int i = 1; i < p; ++i)
    {
        Real fact = rem(i, i);

        for (int j = i + 1; j <= p; ++j)
        {
            fact *= d;
            rem(i, j) = on_binomial_real(j, i) * fact;
        }
    }

    return true;
}
```
