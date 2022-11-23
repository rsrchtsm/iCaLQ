# iCaLQ

#### _LHC Dileption Limits Calculator_
Alpha version

## Introduction
This is a dilepton limits calculator. LHC data is taken from [here](https://www.hepdata.net/record/ins1782650). The FeynRules models of the leptoquarks considered can be found in the repository. The alpha version includes only the U1 leptoquark. Physics details and implementation details can be found in this [paper].

The calculator can be used in two modes: [interactive](#interactive-mode) and [non-interactive](#non-interactive-mode).


## Setting Up

The calculator is written in python3 and only needs the following packages: numpy, pandas, scipy and sympy. The steps involving virtual environment might not be necessary for most people.

Making sure packages requirements are met:
```sh
pip install numpy sympy scipy pandas
```

If you wish to use a virtual environment, do the following:
```sh
cd <icalq-direcory>
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Interactive Mode

To use the calculator in interactive mode,
```sh
python3 icalq.py
```
You will be greeted with iCaLQ art and a list of available commands. The prompt will be `icalq > `. Example commands include:
```
icalq > mass = 1234
icalq > couplings = LM23L LM33R
icalq > significance = 1
icalq > ignore_single_pair = no
icalq > systematic_error = 0.1
icalq > import_model = U1
icalq > status
icalq > help
icalq > initiate
```
The first 4 commands set the input parameters as `<command> = <value>`. The rest do not take any argument.

- `mass` should be an integer between 1000 and 3000 (inclusive).
- `couplings` should list couplings in the format _LMxyL_ or _LMxyR_ where _x_ is the quark number and _y_ is the lepton number and the last _L_ or _R_ denote left or right handedness respectively.
- `significance` takes values 1 or 2.
- `ignore_single_pair` takes values _yes_ or _no_. Input _yes_ means that the single and pair productions will be ignored and this will speed up calculations.
- `systematic-error` denotes the systematic error margin. Default is 10%.
- `import_model` to specify which leptoquark model to use.
- `status` displayes the current values of input parameters.
- `help` displays the list of commands available.
- `initiate` will compute the chi-square polynomial and its minima corresponding to the current values of input parameters.

Once initiated, the calculator goes to input values accepting mode and the prompt changes to ` > `. This prompt will only accept queries of the `<f1> <f2> ... <fn>`, where _\<f1\>_ to _\<fn\>_ are floating point numbers, _n_ is the number of couplings in the parameters and the values are space separated. The expected order of couplings (which is the same as the input order) is mentioned after initiating.

Corresponding to every query, the delta chi-square value will be displayed. Whether this is allowed withing the {1,2} sigma limit is also displayed.

Type `done` to exit query mode. An example query after executing the above commands would be:
```
 > 0.1 0
 > 0.37 0.0001
 > 0.5 0.7
 > done
```

The prompt returns to the input mode and input parameters hold the previous values which can be updated. Finally, to exit the calculator,
```
icalq > exit
```

## Non-interactive Mode

To use the calculator in non-interactive mode, use the tag `-ni` or `--non-interactive`. Note that input cards and query values must be provided in this mode. An example of usage in this mode:
```sh
python3 icalq.py --non-interactive --no-banner --input-card=card.txt --input-values=values.txt --output-yes=yes.csv --output-no=no.csv
```

```sh
python3 icalq.py [options]
```
Options:
- `--help`: Display this help message.
- `--input-card=[filename]`: Input card file. Line 1: mass, line 2: couplings, line 3: ignore_single_pair, line 4: sigma. These are same as input parameter values mentioned in the interactive version.
- `--input-values=[filename]`: Input values to check from the given file. Each line contains a query value. If there are _n_ couplings, then each line would be `<f1> <f2> ... <fn>`, where _\<fi>_ are float values.
- `--non-interactive` or `-ni`: Run in non-interactive mode. This requires input-card and input-values to be specified
- `--no-banner` or `-nb`: icalq banner is not printed.
- `--output-yes=[filename]`: Specify the name of output file (allowed values) (overwrites the existing file). Default: icalq_yes.csv
- `--output-no=[filename]`: Specify the name of output file (disallowed values) (overwrites the existing file). Default: icalq_no.csv

## Caveats

The calculator is in its alpha version and will undergo extensive testing. Currently only U1 leptoquark is used, more leptoquark models will be added later. Any comments regarding the calculator can be emailed to [us](mailto:yash.chaurasia@research.iiit.ac.in).