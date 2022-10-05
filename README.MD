Verification circuit for Lagrange interpolation. Two arrays as input and an array of polynomial coefficients as output. The size of arrays is 2^n. It is recommended to build components, which work for any n, but you can fix n = n' (thus you can fix the number of rows).

${\displaystyle {\begin{aligned}\ell _{j}(x)&={\frac {(x-x_{0})}{(x_{j}-x_{0})}}\cdots {\frac {(x-x_{j-1})}{(x_{j}-x_{j-1})}}{\frac {(x-x_{j+1})}{(x_{j}-x_{j+1})}}\cdots {\frac {(x-x_{k})}{(x_{j}-x_{k})}}\\[10mu]&=\prod _{\begin{smallmatrix}0\leq m\leq k\\m\neq j\end{smallmatrix}}{\frac {x-x_{m}}{x_{j}-x_{m}}}.\end{aligned}}}$

```
//|------|------|------ - ---- | ------ ------ -| ------ ------ | ------ -  ----|
//|  a   | b    |        W2    |         W3     |         W4    |         W5    |
//|------|------|------ c ---- | ------ sx [0] -| ------ sr [0] | ------ r  ----|
//|  2   | 4    |-0.114286     | -189           |21.6           |13.5           |
//|  3   | 8    |0.333333      | 111            |-12.6857       |-12.75         |
//|  9   | 7    |0.0833333     | -19            |2.17143        |3.5            |
//|  7   | 10   |-0.25         | 1              |-0.114286      | -0.25         |
//|------|------|------ ------ | ------ sx [1] -|------ sr [1] -|------ ------ -|
//|      |      |              | -126           |-42            |
//|      |      |              | 95             |31.6667        |
//|      |      |              | -18            |-6             |
//|      |      |              | 1              |0.333333       |
//|------|------|------ ------ | ------ sx [2] -|------ sr [2] -|
//|      |      |              | -42            |-3.5           |
//|      |      |              | 41             |3.41667        |
//|      |      |              | -12            |-1             |
//|      |      |              | 1              |0.0833333      |
//|------|------|------ ------ | ------ sx [3] -|------ sr [3] -|
//|      |      |              | -54            |13.5           |
//|      |      |              | 51             |-12.75         |
//|      |      |              | -14            |3.5            |
//|      |      |              | 1              |-0.25          |
//|------|------|------ ------ | ------ ------ -|------ ------ -|
```