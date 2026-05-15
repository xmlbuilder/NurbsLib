# de Casteljau Bezier

- 컨트롤 포인트들 $P_i$ 에서 시작해서,
- de Casteljau 재귀식
```math
P_{k,i}=(1-u_0)P_{k-1,i}+u_0P_{k-1,i+1}
```
- 를 k=1부터 n까지 반복해서 내려가면,
- 마지막에 남는 $P_{n,0}(u_0)$ 가 바로 Bézier 곡선 위의 실제 점 좌표다.


## Bézier 곡선의 두 가지 표현
- Bernstein 합 표현
```math
C(u)=\sum _{i=0}^nB_{i,n}(u)\, P_i
```
- de Casteljau 재귀 알고리즘
```math
P_{k,i}(u)=(1-u)P_{k-1,i}(u)+uP_{k-1,i+1}(u),\quad C(u)=P_{n,0}(u)
```
- 이 둘이 **완전히 같은 점 C(u)** 를 준다.

- 직접 눈에 보이게 전개.

### 1. n=2 (2차)에서 먼저 확인
- 2차 Bézier:
```math
C(u)=(1-u)^2P_0+2u(1-u)P_1+u^2P_2
```
- de Casteljau:
- 1단계:
```math
P_{1,0}=(1-u)P_0+uP_1
```
```math
P_{1,1}=(1-u)P_1+uP_2
```
- 2단계:
```math
P_{2,0}=(1-u)P_{1,0}+uP_{1,1}
```
- 이제 $P_{2,0}$ 을 전개해보자:

- 즉,
```math
P_{2,0}=C(u)
```
- 2차에서는 두 표현이 완전히 동일함이 직접 보인다.

### 2. 일반 n에 대해 아이디어: “선형 보간의 반복 = Bernstein 조합”
- 핵심 아이디어는:
    - de Casteljau는 선형 보간을 n번 반복하는 알고리즘
    - 선형 보간을 반복하면, 결국 각 P_i에 붙는 계수가 Bernstein 다항식이 된다
- 이를 조금 더 구조적으로 보자.

### 3. de Casteljau를 계수 형태로 다시 쓰기
- 초기:
```math
P_{0,i}=P_i
```
- 1단계(k=1):
```math
P_{1,i}=(1-u)P_{0,i}+uP_{0,i+1}
```
- 이걸 “계수” 관점에서 보면:
```math
P_{1,i}=\sum _{j=0}^n\alpha _{i,j}^{(1)}(u)\, P_j
```
- 실제로는:
- $P_{1,0}=(1-u)P_0+uP_1$
- $P_{1,1}=(1-u)P_1+uP_2$
- …
- 즉, 각 단계에서 계수들이 $(1-u), u$ 로 섞여 나간다.
- k단계에서:
```math
P_{k,i}=\sum _{j=0}^n\alpha _{i,j}^{(k)}(u)\, P_j
```
- 그리고 재귀식:
```math
P_{k,i}=(1-u)P_{k-1,i}+uP_{k-1,i+1}
```
- 계수로 쓰면:
```math
\alpha _{i,j}^{(k)}(u)=(1-u)\, \alpha _{i,j}^{(k-1)}(u)+u\, \alpha _{i+1,j}^{(k-1)}(u)
```
- 초기 조건:
```math
\alpha _{i,j}^{(0)}(u)=\left\{ \, \begin{array}{ll}\textstyle 1&\textstyle i=j\\ \textstyle 0&\textstyle i\neq j\end{array}\right.
``` 

### 4. 이 재귀식은 사실 Bernstein 재귀식과 동일하다
- Bernstein 다항식 재귀식:
```math
B_{i,n}(u)=(1-u)B_{i,n-1}(u)+uB_{i-1,n-1}(u)
```
- 지금 우리가 얻은 계수 재귀식:
```math
\alpha _{i,j}^{(k)}(u)=(1-u)\, \alpha _{i,j}^{(k-1)}(u)+u\, \alpha _{i+1,j}^{(k-1)}(u)
```
- 인덱스를 적절히 바꾸면,
- **k단계에서 i번째 점에 붙는 j번째 컨트롤 포인트의 계수** 가
    결국 어떤 Bernstein 다항식이 된다는 걸 볼 수 있음.
- 특히 마지막 단계 k=n,i=0에서:
```math
P_{n,0}(u)=\sum _{j=0}^n\alpha _{0,j}^{(n)}(u)\, P_j
```
- 이때
```math
  \alpha _{0,j}^{(n)}(u)
```
- 가 바로
```math
$B_{j,n}(u)
```
- 가 된다.
- 즉:
```math
P_{n,0}(u)=\sum _{j=0}^nB_{j,n}(u)\, P_j=C(u)
```

### 5. 직관적으로 정리하면
- de Casteljau는 **컨트롤 포인트들을 u 비율로 계속 섞어가는 과정**
- 이 섞임의 모든 조합을 전개하면,  
    각 $P_i$ 에 붙는 계수는 **u와 (1-u)의 조합** 으로 이루어진 다항식이 된다
- 그 다항식이 바로 Bernstein 다항식 $B_{i,n}(u)$
- 그래서 마지막에 나오는 $P_{n,0}(u)$ 는
```math
P_{n,0}(u)=\sum _{i=0}^nB_{i,n}(u)\, P_i=C(u)
```
- 즉, “Bernstein 합 표현”과 **de Casteljau 알고리즘** 은
    완전히 같은 곡선, 같은 점을 다른 방식으로 계산하는 것일 뿐이다.

## 6. 한 문장 요약
- de Casteljau는 컨트롤 포인트를 선형 보간으로 반복해서 섞는 알고리즘이고,  
    그 반복 보간에서 자연스럽게 생기는 계수들이 바로 Bernstein 다항식이기 때문에,  
- 마지막 점 $P_{n,0}(u)$ 는
```math
P_{n,0}(u)=\sum _{i=0}^nB_{i,n}(u)\, P_i=C(u)
```
- 로, Bernstein 표현과 완전히 동일한 결과를 준다.

```rust
fn de_casteljau_scalar(ctrl: &[Real], t: Real) -> Real {
    let mut tmp = ctrl.to_vec();
    let n = tmp.len();
    if n == 0 {
        return 0.0;
    }
    for r in 1..n {
        for i in 0..(n - r) {
            tmp[i] = (1.0 - t) * tmp[i] + t * tmp[i + 1];
        }
    }
    tmp[0]
}
```
```cpp
Real on_de_casteljau_scalar(
    const std::vector<Real>& ctrl,
    const Real t)
{
    const std::size_t n = ctrl.size();

    if (n == 0)
        return 0.0;

    std::vector<Real> tmp = ctrl;

    for (std::size_t r = 1; r < n; ++r)
    {
        for (std::size_t i = 0; i < n - r; ++i)
        {
            tmp[i] =
                (1.0 - t) * tmp[i]
                + t * tmp[i + 1];
        }
    }

    return tmp[0];
}
```
```rust
fn de_casteljau_point3(ctrl: &[Point3D], t: Real) -> Point3D {
    let mut tmp = ctrl.to_vec();
    let n = tmp.len();
    if n == 0 {
        return Point3D::zero();
    }
    for r in 1..n {
        for i in 0..(n - r) {
            tmp[i] = Point3D {
                x: (1.0 - t) * tmp[i].x + t * tmp[i + 1].x,
                y: (1.0 - t) * tmp[i].y + t * tmp[i + 1].y,
                z: (1.0 - t) * tmp[i].z + t * tmp[i + 1].z,
            };
        }
    }
    tmp[0]
}
```
```cpp
Point3D on_de_casteljau_point3(
    const std::vector<Point3D>& ctrl,
    const Real t)
{
    const std::size_t n = ctrl.size();

    if (n == 0)
        return Point3D::Zero();

    std::vector<Point3D> tmp = ctrl;
    for (std::size_t r = 1; r < n; ++r)
    {
        for (std::size_t i = 0; i < n - r; ++i)
        {
            tmp[i].x =
                (1.0 - t) * tmp[i].x
                + t * tmp[i + 1].x;

            tmp[i].y =
                (1.0 - t) * tmp[i].y
                + t * tmp[i + 1].y;

            tmp[i].z =
                (1.0 - t) * tmp[i].z
                + t * tmp[i + 1].z;
        }
    }
    return tmp[0];
}
```
- 을 그대로 코드로 옮기면 바로 이 형태가 된다.
---

## 🟦 1. 입력: 컨트롤 포인트 배열
- 예를 들어 $ctrl = [P_0,P_1,\dots ,P_n]$
- 수식에서는:
```math
P_{0,i}=P_i
```
- 코드에서는:
```rust
let mut tmp = ctrl.to_vec();
```

- 즉, $tmp[i] = P_{0,i}$
    - 초기 삼각형의 맨 위 줄을 그대로 복사

## 🟦 2. 바깥 루프: k = 1..n (삼각형의 단계)
- 수식에서는:
```math
k=1,2,\dots ,n
```
- 코드에서는:
```rust
for r in 1..n {
```

- 여기서 r이 바로 수식의 k에 해당한다.
- 즉:
    - r = 1 → 첫 번째 보간
    - r = 2 → 두 번째 보간
    - …
    - r = n → 마지막 보간

## 🟦 3. 안쪽 루프: i = 0..(n-k)
- 수식에서는:
```math
i=0,1,\dots ,n-k
```
- 코드에서는:
```rust
for i in 0..(n - r) {
```

- 완전히 동일하다.

## 🟦 4. 핵심 보간식:
```math
tmp[i] = (1-t) * tmp[i] + t * tmp[i+1]
```
- 수식:
```math
P_{k,i}(t)=(1-t)P_{k-1,i}(t)+tP_{k-1,i+1}(t)
```
- 코드:
```math
tmp[i] = (1.0 - t) * tmp[i] + t * tmp[i + 1];
```

- 이건 수식을 1:1로 옮긴 것이다.
- $tmp[i] → P_{k-1,i}$
- $tmp[i+1] → P_{k-1,i+1}$
- $tmp[i] = … → P_{k,i}$ 로 업데이트
- 즉, 삼각형의 다음 줄을 계산하는 과정을 그대로 구현한 것.

## 🟦 5. 반복이 끝나면 tmp[0]만 남는다
- 수식에서는:
```math
C(t)=P_{n,0}(t)
```
- 코드에서는:
```rust
tmp[0]
```

- 즉:
- 첫 번째 줄: $P0, P1, P2, …, Pn$
- 두 번째 줄: $P_{1,0}, P_{1,1}, …$
- 세 번째 줄: $P_{2,0}, P_{2,1}, …$
- …
- 마지막 줄: $P_{n,0}$
- 이 마지막 하나가 곡선 위의 점.

## 🟦 6. Point3D 버전도 완전히 동일
- 스칼라 버전:
```rust
tmp[i] = (1-t)*tmp[i] + t*tmp[i+1];
```

- Point3D 버전:
```rust
tmp[i] = Point3D {
    x: (1-t)*tmp[i].x + t*tmp[i+1].x,
    y: (1-t)*tmp[i].y + t*tmp[i+1].y,
    z: (1-t)*tmp[i].z + t*tmp[i+1].z,
};
```

- 즉, 벡터 보간을 성분별로 수행한 것일 뿐 수식은 완전히 동일하다.

## 🟦 7. 전체 대응표
| 수식 표현                                   | 코드 표현                                   |
|---------------------------------------------|----------------------------------------------|
| $P_{0,i} = P_i$                             | tmp = ctrl.to_vec()                          |
| k = 1 .. n                                  | for r in 1..n                                |
| i = 0 .. (n - k)                            | for i in 0..(n - r)                          |
| $P_{k,i} = (1 - t) * P_{k-1,i} + t * P_{k-1,i+1}$ | tmp[i] = (1.0 - t) * tmp[i] + t * tmp[i+1] |
| $C(t) = P_{n,0}$                              | tmp[0]                                       |


- 즉, 코드가 de Casteljau 알고리즘을 100% 정확하게 구현하고 있다.

## 🟦 8. 한 문장 요약
- 이 코드는 de Casteljau 수식
```math
P_{k,i}=(1-t)P_{k-1,i}+tP_{k-1,i+1}
```
- 을 그대로 두 중첩 for문으로 구현한 것이며, 마지막 tmp[0]이 곡선 위의 점 C(t)이 된다.

---
