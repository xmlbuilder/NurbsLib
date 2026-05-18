# LU 기반 선형 시스템 해법 분석

다음 함수는 부분 피벗(Partial Pivoting)을 사용하는 가우스 소거 기반의 선형 시스템 해법입니다.

```cpp
bool maths::on_solve_linear_system(
    Matrix& A,
    std::vector<double>& b)
```

이 함수는 다음 선형 시스템을 풉니다.

$$
A x = b
$$

---

# 전체 알고리즘 구조

이 함수는 내부적으로 다음 과정을 수행합니다.

1. Pivot 선택
2. 행 교환(Row Swap)
3. Elimination (LU 분해 진행)
4. Forward Substitution 효과를 b에 직접 반영
5. Back Substitution 수행
6. 최종 해를 b에 저장

---

# 1. Pivot 선택

현재 열에서 가장 큰 절댓값을 가진 행을 찾습니다.

```cpp
int max_i = k;
double max_v = std::fabs(A[k][k]);

for (int i = k + 1; i < N; ++i)
{
    if (const double v = std::fabs(A[i][k]); v > max_v)
    {
        max_v = v;
        max_i = i;
    }
}
```

목적은 다음 식의 분모가 너무 작아지는 것을 방지하는 것입니다.

$$
m=\frac{A_{ik}}{A_{kk}}
$$

작은 pivot은 수치 오차를 크게 증가시킵니다.

---

# 2. Row Swap (행 교환)

Pivot 행이 현재 행과 다르면 행을 교환합니다.

```cpp
A.SwapRows(k, max_i);
std::swap(b[k], b[max_i]);
```

행렬만 바꾸면 안 되고 우변 벡터도 반드시 같이 교환해야 합니다.

원래 시스템:

$$
A x = b
$$

행 교환 후:

$$
P A x = P b
$$

여기서 \(P\)는 permutation matrix 입니다.

---

# 3. Elimination (전진 소거)

핵심 코드:

```cpp
const double akk = A[k][k];

for (int i = k + 1; i < N; ++i)
{
    const double m = A[i][k] / akk;

    A(i, k) = m;

    for (int j = k + 1; j < N; ++j)
        A[i][j] -= m * A[k][j];

    b[i] -= m * b[k];
}
```

---

# 3-1. Multiplier 계산

$$
m=\frac{A_{ik}}{A_{kk}}
$$

이 값은 LU 분해에서 Lower matrix \(L\) 의 원소입니다.

코드에서는 아래쪽 삼각 영역에 직접 저장합니다.

```cpp
A(i, k) = m;
```

즉:

$$
L_{ik}=m
$$

---

# 3-2. Upper Matrix(U) 생성

다음 식으로 아래쪽 항들을 제거합니다.

$$
A_{ij}=A_{ij}-mA_{kj}
$$

즉:

```cpp
A[i][j] -= m * A[k][j];
```

결과적으로 행렬은 위삼각 행렬 형태로 변합니다.

```math
U=
\begin{bmatrix}
* & * & * \\
0 & * & * \\
0 & 0 & *
\end{bmatrix}
```

---

# 4. 실제 Forward Substitution은 어디서 수행되는가?

이 함수의 가장 중요한 부분입니다.

일반적인 LU 시스템은 보통 다음 두 단계로 풉니다.

## Step 1

$$
L y=b
$$

(Forward Substitution)

## Step 2

$$
U x=y
$$

(Back Substitution)

---

하지만 현재 구현은 별도의 \(y\) 벡터를 만들지 않습니다.

대신 elimination 과정 중에 직접 b를 갱신합니다.

```cpp
b[i] -= m * b[k];
```

수식으로 쓰면:

$$
b_i=b_i-mb_k
$$

이 연산이 반복되면서 b는 점점 변형됩니다.

최종적으로 elimination이 끝나면:

$$
b \rightarrow y
$$

상태가 됩니다.

즉:

```text
Forward Substitution이 elimination 과정 안에 포함되어 있다.
```

라고 볼 수 있습니다.

---

# 5. Back Substitution

이제 위삼각 행렬 시스템을 풉니다.

```cpp
for (int i = N - 1; i >= 0; --i)
{
    double sum = b[i];

    for (int j = i + 1; j < N; ++j)
        sum -= A[i][j] * b[j];

    b[i] = sum / A[i][i];
}
```

---

# 5-1. 수식 형태

현재 시점에서:

$$
U x = y
$$

상태입니다.

위삼각 행렬은 아래에서부터 역순으로 풉니다.

---

## 마지막 행

$$
U_{nn}x_n=y_n
$$

따라서:

$$
x_n=\frac{y_n}{U_{nn}}
$$

---

## 일반적인 행

$$
U_{ii}x_i+\sum_{j=i+1}^{n-1}U_{ij}x_j=y_i
$$

정리하면:

$$
x_i=
\frac{
y_i-\sum_{j=i+1}^{n-1}U_{ij}x_j
}{
U_{ii}
}
$$

---

# 코드와 수식 대응

## 코드

```cpp
sum -= A[i][j] * b[j];
```

## 수식

$$
y_i-\sum_{j=i+1}^{n-1}U_{ij}x_j
$$

---

## 코드

```cpp
b[i] = sum / A[i][i];
```

## 수식

$$
x_i=
\frac{
y_i-\sum_{j=i+1}^{n-1}U_{ij}x_j
}{
U_{ii}
}
$$

---

# 최종 메모리 상태

함수 종료 후:

## 행렬 A 내부

- Upper Triangle → \(U\)
- Lower Triangle → \(L\) 의 multiplier

즉:

$$
A \rightarrow LU
$$

형태로 저장됩니다.

---

## 벡터 b 내부

초기에는 RHS 벡터였지만:

$$
b \rightarrow y \rightarrow x
$$

로 변환됩니다.

최종적으로:

```cpp
b == x
```

입니다.

---

# 전체 알고리즘 요약

## 입력

$$
A x=b
$$

---

## LU 분해 수행

$$
P A=L U
$$

---

## Elimination 중

$$
b \rightarrow y
$$

즉:

$$
L y=P b
$$

를 암묵적으로 수행

---

## 마지막 단계

$$
U x=y
$$

를 Back Substitution으로 해결

---

# 핵심 포인트

이 구현에서:

```cpp
b[i] -= m * b[k];
```

이 부분이 사실상 Forward Substitution 역할입니다.

그리고:

```cpp
for (int i = N - 1; i >= 0; --i)
```

이 부분이 Back Substitution 입니다.

따라서 이 함수는:

```text
LU 분해 + implicit forward substitution + back substitution
```

구조라고 볼 수 있습니다.

---
## 소스 코드
```cpp
bool maths::on_solve_linear_system(
    Matrix& A,
    std::vector<double>& b,
    double zero_tol)
{
    const int N = A.RowCount();

    if (!A.IsValid())
        return false;

    if (N == 0 || A.ColCount() != N)
        return false;

    if (static_cast<int>(b.size()) != N)
        return false;

    for (int k = 0; k < N; ++k)
    {
        int pivot_row = k;
        double pivot_abs = std::fabs(A(k, k));

        for (int i = k + 1; i < N; ++i)
        {
            const double v = std::fabs(A(i, k));
            if (v > pivot_abs)
            {
                pivot_abs = v;
                pivot_row = i;
            }
        }

        if (pivot_abs <= zero_tol)
            return false;

        if (pivot_row != k)
        {
            if (!A.SwapRows(k, pivot_row))
                return false;

            std::swap(b[k], b[pivot_row]);
        }

        const double pivot = A(k, k);

        for (int i = k + 1; i < N; ++i)
        {
            const double m = A(i, k) / pivot;

            A(i, k) = m; // L 저장

            for (int j = k + 1; j < N; ++j)
                A(i, j) -= m * A(k, j);

            b[i] -= m * b[k];
        }
    }

    for (int i = N - 1; i >= 0; --i)
    {
        const double diag = A(i, i);

        if (std::fabs(diag) <= zero_tol)
            return false;

        double sum = b[i];

        for (int j = i + 1; j < N; ++j)
            sum -= A(i, j) * b[j];

        b[i] = sum / diag;
    }

    return true;
}
```
---
