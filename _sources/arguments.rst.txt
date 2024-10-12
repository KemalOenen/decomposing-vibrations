Arguments of Nomodeco
======================

The output an usage of Nomodeco can be modified using different arguments. One can return a list of all possible arguments with::
    nomodeco --help

**Arguments-Table**

+------------------------+---------------------------+--------------------------------------------------------+
| Argument               | Attributes                | Function                                               |
|                        |                           |                                                        |
+========================+===========================+========================================================+
| \-\-help               | None                      | shows a list of possible arugments                     |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-log                | None                      | set an additonal log file                              |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-matrix_opt         | ved, diag, contr (default)| choose the matrix used for optimization                |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-penalty1           |[INFTFREQ-PENALTY]         | penalty value for assymetric intrinsic frequencys      |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-penality2          |[INTFC-PENALTY]            | penality value for unphysical contributions            |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-heatmap            | VED, PED, contr           | return a heatmap for the specified matrix              |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-csv                | VED, PED, contr           | return a csv for the specified matrix                  |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-latex_tab          | None                      | generate additional latex repr. of contribution table  |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-molpro             | molpro.out file           | molpro.out files can be specified and used as input    |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-gv                 | gaussian.log file         | gaussian.log files can be specified and used as input  |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-pymolpro           | None                      | Use pymolpro integration to run molpro calculation     |
+------------------------+---------------------------+--------------------------------------------------------+
| \-\-comb               | {1,2,3}                   | Adjust which coordinates get added for set generation  |
+------------------------+---------------------------+--------------------------------------------------------+

**Explanation of comb**

With the parameter comb one can adjust which IC types get included in the generation of IC sets for hydrogen-bonded clusters:

* comb = 1, here all covalent as well as ICs which containing a hydrogen bond are included
* comb = 2, here all covalent and hydrogen bonded as well as acceptor donor bonds and angles are included
* comb = 3, here all covalent asnd hydrogen and all acceptor donor bonds are included (default)

