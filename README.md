# MasseyRamanujan

The Ramanujan Machine is an algorithmic approach to discover new mathematical conjectures. For the time being, the project is focused in number theory, specifically on finding formulas relating fundamental constants like pi, e, and the Riemann zeta function values to various continued fractions.

For more information, please go to [RamanujanMachine.com](https://www.RamanujanMachine.com).

## Installation

Clone the repo and install the package. If you have pip, it can be done by running
```
pip install -e .
```
under the main folder. That's it, you are now ready to discover new conjectures.

## Running the code

To use the RamanujanMachine, you'll have to write a short script that combines three elements:

- Left Hand Side Hash Table (`LHSHashTable`) - A data structure that holds expressions made from the required constant.
  for example, the following code will generate all Mobius transforms of `e`, keeping the parameters between -5 and 5.
  Also, it will save the generated domain under `saved_hash`, and will load the data from it on the next execution.
```python
 from ramanujan.LHSHashTable import LHSHashTable
 from ramanujan.constants import g_const_dict

 saved_hash = 'e_lhs_dept5_db'
 lhs_search_limit = 5
 lhs = LHSHashTable(
    saved_hash,
    lhs_search_limit,
    [g_const_dict['e']])
```
- A Polynomial Domain (Any class under `poly_domains`) - A definition for a family of `an` and `bn` polynomials. 
Based on those families GCFs will be generated. 
  
  The following code will generate the simplest domain, where an and bn are polynomials of degree 2, and the 
  coefficients have no connection, so every permutation of integers between -5 and 5 for the coefficients will be 
  generated
 ```python
  from ramanujan.poly_domains.CartesianProductPolyDomain import CartesianProductPolyDomain
  
  poly_search_domain = CartesianProductPolyDomain(
    2, [-5, 5],
    2, [-5, 5])
```

- An enumerator (any class under `enumerators`) - The glue that holds the two. This is where MITM algorithm resides.

This class will calculate the gcf and decide which conjectures pass the first calculation, and which pass the second,
as described in our paper.
For example, creating `EfficientGCFEnumerator` using the `LHSHashTable` and `poly_domain` defined above:
```python
from ramanujan.enumerators.EfficentGCFEnumerator import EfficentGCFEnumerator

enumerator = EfficentGCFEnumerator(
    lhs,
    poly_search_domain,
    [g_const_dict['e']],
    lhs_search_limit
    )
```


And thats it! start your execution by running:
```python
results = enumerator.full_execution()
```


### Cool examples
Examples given here can be found under `scripts/paper_results`
#### e
To find a few formulas for e, use the following poly domain:
```python
poly_search_domain = CartesianProductPolyDomain(
    2, [-5, 5],
    2, [-5, 5])
```

#### pi
To find a few formulas for pi, use:
```python
poly_search_domain = CartesianProductPolyDomain(
    1, [-13, 13],
    2, [-11, 11])
```

#### Riemann Zeta function at 3 (Aépry's constant)
To find a few formulas related to the Riemann zeta function at 3 (Zeta of 3 is called 
[Apéry's constant](https://www.wikiwand.com/en/Ap%C3%A9ry%27s_constant) and has a role in the electron's gyromagnetic
ratio), use:
```python
poly_search_domain = Zeta3Domain1(
    [(2, 2), (1, 1), (-20, 20), (-20, 20)],
    (-20, -1))
```

Now that you've seen how to run the basic code, you can tweak the search parameters and find new conjectures of your own.
To do so, please read the next section.

### Tweaking the search parameters

If you wish to tweak the searched series, you can create a new class that extends `AbstractPolyDomains`, and defines
your new polynomial families. If there is no complex connection between an and bn, you can extend 
`CartesianProductPolyDomain`.

I'll add an extended example soon