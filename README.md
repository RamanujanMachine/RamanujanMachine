# MasseyRamanujan

The Ramanujan Machine is an algorithmic approach to discover new mathematical conjectures. For the time being, the project is focused in number theory, specifically on finding formulas relating fundamental constants like pi, e, and the Riemann zeta function values to various continued fractions.

For more information, please go to [RamanujanMachine.com](https://www.RamanujanMachine.com).

## Installation

Clone the repo and install the requirements in the source/requirements.txt file. If you have pip, it can be done by running
```
pip install -r requirements.txt
```
under the source/ folder. That's it, you are now ready to discover new conjectures.

## Running the code

The source code exists in the source/ folder and should be run from there. Results that are generated in the examples below will be printed on the screen as well as to a LaTeX and a PDF under the source/results/ folder for your convenience.

### Cool examples

#### e
To find a few formulas for e, run (don't forget to be under the source/ folder when you run it)
```python
python main.py MITM_RF -lhs_constant e -num_of_cores 1 -lhs_search_limit 5 -poly_a_order 2 -poly_a_coefficient_max 5 -poly_b_order 2 -poly_b_coefficient_max 5
```

#### pi
To find a few formulas for pi, run (don't forget to be under the source/ folder when you run it)
```python
python main.py MITM_RF -lhs_constant pi -num_of_cores 1 -lhs_search_limit 20 -poly_a_order 2 -poly_a_coefficient_max 13 -poly_b_order 3 -poly_b_coefficient_max 11 --polynomial_shift1_bn
```

#### Riemann Zeta function at 3 (Aépry's constant)
To find a few formulas related to the Riemann zeta function at 3 (Zeta of 3 is called [Apéry's constant](https://www.wikiwand.com/en/Ap%C3%A9ry%27s_constant) and has a role in the electron's gyromagnetic ratio), run:
```python
python main.py MITM_RF -lhs_constant zeta -function_value 3 -num_of_cores 2 -lhs_search_limit 14 -poly_a_order 3 -poly_a_coefficient_max 20 -poly_b_order 3 -poly_b_coefficient_max 20 --zeta3_an --zeta_bn
```

#### Catalan constant
To find a few formulas related to the [catalan constant](https://www.wikiwand.com/en/Catalan%27s_constant), you can run the following code. This one takes a bit longer generate the hash table for and make take a few minutes.
```python
python main.py MITM_RF -lhs_constant catalan pi-acosh_2 -num_of_cores 1 -lhs_search_limit 8 -poly_a_order 3 -poly_a_coefficient_max 15 -poly_b_order 2 -poly_b_coefficient_max 5 --catalan_bn
```

Now that you've seen how to run the basic code, you can tweak the search parameters and find new conjectures of your own. To do so, please read the next section.

## Testing for performance

To test the speed of various modules or functions, use the inherent cProfile module in python, e.g.
```python
python -m cProfile -s tottime main.py MITM_RF -lhs_constant e -num_of_cores 1 -lhs_search_limit 5 -poly_a_order 3 -poly_a_coefficient_max 10 -poly_b_order 2 -poly_b_coefficient_max 10
```

In addition, we added a decorator called measure_performance under utils that allows you to test individual functions by adding @measure_performance before the function definition.

### Tweaking the search parameters

Under the source/ folder,
```python
python main.py
```
runs the code. The infrastructure supports various algorithms for discovery of constants and for now can only be run using the MITM_RF toggle.

##### MITM_RF module: 
this is our new MITM implementation. The program will "mine" new Continued Fraction conjectures of the type:
![LHS_RHS](images/LHS_RHS.png)
The code let's you control the equation space scanned by the algorithm. To get more information about what you can control and tweak, run
```python
python main.py MITM_RF -h
```

Parameters that you can currently control without changing the code itself include:

* -lhs_constant {zeta,e,pi,catalan,golden_ratio,khinchin,euler-mascheroni,pi-acosh_2} [{zeta,e,pi,catalan,golden_ratio,khinchin,euler-mascheroni,pi-acosh_2} ...] constants to search for - initializing the left-hand-side hash table
* -function_value FUNCTION_VALUE Which value of the function are we assessing (assuming LHS constant takes an arguments)
* -lhs_search_limit LHS_SEARCH_LIMIT The limit for the LHS coefficients
* -num_of_cores NUM_OF_CORES The number of cores to run on
* -poly_a_order POLY_A_ORDER the number of free coefficients for {a_n} series
* -poly_a_coefficient_max POLY_A_COEFFICIENT_MAX The maximum value for the coefficients of the {a_n} polynomial
* -poly_b_order POLY_B_ORDER the number of free coefficients for {b_n} series
* -poly_b_coefficient_max POLY_B_COEFFICIENT_MAX The maximum value for the coefficients of the {b_n} polynomial
* custom {a_n} series generator:
  if defined, poly_a_order is ignored.if not defined the default polynomial
  will be used

  - --zeta3_an            Generator3[x3_, x0_] := {x0, 2 *x0 + x3, 3*x3, 2*x3}.
                        this was found to be useful for zeta3 searches
  - --zeta5_an            Generator5[x5_, x3_, x0_] := {x0, 2 *x0 + x3 - 2 *x5,
                        3*x3 - 5*x5, 2*x3, 5*x5, 2*x5}
  - --polynomial_shift1_an
                        a[n] = m(m(...(x[1]*m + x[0]) + x[2]) + ...) + x[k],
                        where m=n+1
  - --polynomial_an       a[n] = n(n(...(x[1]*n + x[0]) + x[2]) + ...) + x[k].
                        this is the default generator

* custom {b_n} series generator:
  if defined, poly_b_order is ignored. if not defined the default polynomial
  will be used

  - --zeta_bn             b[n] = x[0]*(n+1)^d - x[1]*(n+1)^(d-1). where
                        d=function_value this was found to be useful for zeta
                        values searches.
  - --catalan_bn          x[0]*(2*n+1)^4 + x[1]*(2*n+1)^3
  - --polynomial_shift1_bn
                        b[n] = m(m(...(x[1]*m + x[0]) + x[2]) + ...) + x[k],
                        where m=n+1
  - --polynomial_shift2n1_bn
                        b[n] = m(m(...(x[1]*m + x[0]) + x[2]) + ...) + x[k],
                        where m=2*n+1
  - --integer_factorization_bn
                        b[n] = x[0]*(term1)^d1*(term2)^d2*..., where
                        sum(d)=x[1]this is a unique generator using
                        combination permutations instead of cartesian product.
                        current terms are: (n+1), (n+2), (n+3), (2n), (2n-1),
                        (2n+3), (2n+5)
  - --polynomial_bn       b[n] = n(n(...(x[1]*n + x[0]) + x[2]) + ...) + x[k].
                        this is the default generator


more detailed information regarding this module can be found in documentation/enumerate_over_gcf.pdf
