## Branch: Run2_UL_NanoAOD
### CMSSW version: CMSSW_10_6_30
**Instructions**

1. Enable CentOS7 on singularity by running the following command in your home or work area.

```
cmssw_el7
```

2. Initialize CMSSW_10_6_30.
```
Singularity > cmsrel CMSSW_10_6_30
Singularity > cd CMSSW_10_6_30/src
Singularity > cmsenv
```

3. Clone this repository from GitHub and switch to this branch.

4. Run ```init.sh```.
```
. init.sh
```