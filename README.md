# PopHWE
PopHWE is a statistical method and software tool to test for Hardy–Weinberg Equilibrium (HWE) in the presence of population stratification. It separates HWE deviations caused by population structure from those caused by other reasons, especially genotype data errors. See [PopHWE documentation](https://github.com/jimmy-penn/PopHWE/blob/main/PopHweDocumentation.md) for software documentation.

Developer: [Yumi (Jimmy) Jin, PhD](mailto:yumijin@hotmail.com).

### Required libraries and files

PopHWE is currently implemented in python. The PopHWE repository can be downloaded using the following command:

``` sh
$ git clone https://github.com/jimmy-penn/pophwe.git
```

### Running medium tests

Test scripts and test cases are placed under `medium_testing` directory. Test cases are saved in `test_manifest.txt`. Perl script `test_pophwe.pl` is used for manually running these test cases. 
1. Make sure the Perl and Python interpreters `perl` and `python` can be called from current directory `.`. 
1. Under `medium_testing` directory, execute: `perl test_pophwe.pl test_manifest.txt`. 
1. If a another python interpreter (not the default one) should be used, pass it to the script by executing: `perl test_pophwe.pl test_manifest.txt --python <python intepreter>`. 
1. If source code is updated, update `test_manifest.txt` to add new test cases or modify existing cases, and execute `perl test_pophwe.pl test_manifest.txt --ub` to update the baseline. 
1. Check the baseline files and make sure they are all correct, then execute `test_pophwe.pl test_manifest.txt` (without option --ub) again. 
1. Make sure all test cases pass.
