## 1. Perform Scan calculation
### Input file 
Route section example
```
#p opt=modredundant ucam-b3lyp/aug-cc-pvdz
```
Contraints (in this case the bond length will decrease)
```
B 4 3 S 16 -0.05 
```

Input file example
```
#p opt=modredundant ucam-b3lyp/aug-cc-pvdz

t1-opt-from-traj-30

0 3
 C                 -1.01731900    0.76494100    0.09660500
 C                 -1.01705900   -0.76522100   -0.09664500
 C                  0.40644300   -1.16730000    0.08019200
 C                  0.40610100    1.16744200   -0.07991200
B 4 3 S 16 -0.05 

```

### Analysis
1. Select the point with high energy to perform TS calculation
   
