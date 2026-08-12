# Apply the PREP2 upper-limb prognosis algorithm

Implements the day-3 PREP2 decision tree from SAFE score, age,
motor-evoked potential status, and NIH Stroke Scale score.

## Usage

``` r
prep2(safe_score, age, mep_status = NULL, nihss = NULL)
```

## Arguments

- safe_score:

  Shoulder Abduction and Finger Extension MRC sum, an integer from 0 to
  10.

- age:

  Age in years.

- mep_status:

  Optional MEP status: logical, `"present"`/`"MEP+"`, or
  `"absent"`/`"MEP-"`.

- nihss:

  Optional day-3 NIHSS total, an integer from 0 to 42.

## Value

A `prep2_result` containing category, decision pathway, inputs, and a
clinical description.

## References

Stinear CM, Byblow WD, Ackerley SJ, Smith MC, Stinear JW, Barber PA
(2017). PREP2: A biomarker-based algorithm for predicting upper limb
function after stroke. *Annals of Clinical and Translational Neurology*,
4:811-820. [doi:10.1002/acn3.488](https://doi.org/10.1002/acn3.488)

## Examples

``` r
prep2(safe_score = 8, age = 65)
#> PREP2 prognosis: Excellent 
#>   Pathway: SAFE>=5 -> age<80 
#>   Full or near-full upper-limb recovery is expected. 
prep2(2, 70, mep_status = "absent", nihss = 5)
#> PREP2 prognosis: Limited 
#>   Pathway: SAFE<5 -> MEP- -> NIHSS<7 
#>   Modest recovery is expected; emphasize task-specific training. 
```
