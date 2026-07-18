
## Why Knot Insertion Produces a Bezier Curve

### 1. NURBS curve

- A NURBS/B-spline curve is

$$
C(u)=\sum_{i=0}^{n}N_{i,p}(u)P_i
$$

where

- $N_{i,p}(u)$ : B-spline basis
- $P_i$ : control points
- $p$ : degree

---

### 2. Knot insertion preserves geometry

- After knot insertion,

$$
C'(u)=C(u)
$$

- Only the knot vector and control points change.

---

### 3. Knot multiplicity

- If the multiplicity of knot $u_k$ is $r$,

- the number of additional insertions required is

$$
p-r
$$

- After insertion,

$$
r'=p+1
$$

---

### 4. Continuity

- For a knot of multiplicity $r$

$$
C^{\,p-r}
$$

- continuity exists.

- When

$$
r=p+1
$$

- the continuity becomes

$$
C^{-1}
$$

meaning the polynomial pieces become completely independent.

---

### 5. Why it becomes a Bezier curve

- A Bezier curve is

$$
B(t)=
\sum_{i=0}^{p}
B_i^p(t)Q_i
$$

- defined on

$$
0\le t\le1.
$$

- After both span boundaries have multiplicity $p+1$,

- the B-spline basis over that span is exactly identical to the Bernstein basis.

- Therefore

$$
\boxed{
\text{Clamped B-spline Span} =
\text{Bezier Curve}
}
$$

---

### 6. Parameter mapping

- If the original span is

$$
[u_a,u_b]
$$

- then

$$
u=u_a+t(u_b-u_a)
$$

- and

$$
t=\frac{u-u_a}{u_b-u_a}.
$$

---

### 7. Knot insertion formula

- Each inserted control point is

$$
Q_i=(1-\alpha_i)P_{i-1}+\alpha_iP_i
$$

- with

$$
\alpha_i=
\frac{u-u_i}
{u_{i+p}-u_i}.
$$

---

### 8. Algorithm

```text
For every non-zero span

    Check multiplicity
          ↓
Insert missing knots
          ↓
Boundary multiplicity = p+1
          ↓
Extract p+1 control points
          ↓
Bezier curve
```

---
