### Input data

- Please create the file in CSV format.
- Enter the substrate's InChI or SMILES in the first column.
- Enter the reaction conditions in the second column and subsequent columns. Positive examples are 1, negative examples are -1, and unlabeled examples are 0.
- Enter descriptors in the columns following the reaction conditions. Descriptors are not required when calculating ECFP OnBits, ECFP Tanimoto, or ECFP ConvexHull.

- If targeting phenols, you can download and use the 28-dimensional dataset from [here](https://github.com/poclab-web/streamlit-substrate-scope-metrics/blob/main/data/phenols_28d.csv).


#### Example 1 (2 reaction conditions, 3 dimensions)

| InChI                                                                   | condition_1 | condition_2 |         HOMO | MolLogP |   MolWt |
| :---------------------------------------------------------------------- | ----------: | ----------: | -----------: | ------: | ------: |
| InChI=1S/C14H22O/c1-13(2,3)10-7-8-12(15)11(9-10)14(4,5)6/h7-9,15H,1-6H3 |           1 |           1 | -7.285308122 |  3.9872 | 206.329 |
| InChI=1S/C10H14O/c1-4-9-5-7(2)8(3)6-10(9)11/h5-6,11H,4H2,1-3H3          |           0 |           1 | -7.107889794 | 2.57144 | 150.221 |
| InChI=1S/C10H14O/c1-10(2,3)8-4-6-9(11)7-5-8/h4-7,11H,1-3H3              |          -1 |           0 | -7.359867358 |  2.6897 | 150.221 |
| InChI=1S/C8H10O/c1-6-3-4-8(9)7(2)5-6/h3-5,9H,1-2H3                      |           1 |           1 |  -7.21918442 | 2.00904 | 122.167 |
| InChI=1S/C11H16O/c1-8-5-6-10(12)9(7-8)11(2,3)4/h5-7,12H,1-4H3           |           1 |           1 |   -7.1838096 | 2.99812 | 164.248 |


#### Example 2 (4 reaction conditions, no descriptors)

| SMILES                                           | K3FeCN6 | MesAcr+BF4- | CuCl | Diacetyl |
| :----------------------------------------------- | ------: | ----------: | ---: | -------: |
| CC(C)(C)c1ccc(O)c(C(C)(C)C)c1                    |       1 |           1 |    1 |        1 |
| COc1cc(C=O)ccc1O                                 |      -1 |           0 |    0 |        0 |
| COC(=O)C(Cc1ccc(O)c(C(C)(C)C)c1)NC(=O)OCc1ccccc1 |       0 |           1 |    0 |        0 |
| CCCCc1ccc(O)cc1                                  |       1 |           0 |    0 |        0 |
| Cc1ccc(O)c(C)c1                                  |       1 |           1 |    0 |        1 |
| COc1cc(C#N)ccc1O                                 |      -1 |           0 |    0 |        0 |
| Cc1cccc(-c2ccc(O)c(C(C)(C)C)c2)c1                |       0 |           0 |    1 |        0 |
| Cc1cc(O)c(C(C)C)cc1C                             |       0 |           1 |    0 |        1 |
| CC(C)(c1ccccc1)c1ccc(O)c(C(C)(C)c2ccccc2)c1      |       0 |           0 |    1 |        0 |


---
### Operation Guide

#### Data entry
- Enter data on the :material/home: Home page.
- Standardization is required for Range, Distance, and ConvexHull calculations. When the “Standardize” toggle is enabled, the input data will be standardized.
- Enter the number of reaction conditions included in the input data and upload the input data in CSV format.
- Once the input data is loaded, the number of substrates, descriptors, and names of reaction conditions will be displayed.


#### ECFP Settings
- You can change the radius and the number of bits in :material/settings: Advanced settings at the bottom of the :material/home: Home page. The default values are radius 2 and 2048 bits.

#### Evaluation of the Substrate Scope
- Select a metric from the sidebar on the left side of the screen.
- All graphs can be zoomed in or out.

> :material/warning: Warning  
> - ECFP Tanimoto calculations require at least two positive examples.  
> - ECFP ConvexHull and ConvexHull calculations require at least three positive examples.  
> - ConvexHull calculations require descriptors to be two-dimensional or higher.  

#### Reset
- Clicking the “Reset” button at the bottom of the :material/home: Home will clear the input data and all cached information.
