# Certifying the Anosov property

This Python code verifies that a specific surface subgroup of SL(3,R) is Anosov, based on the criteria in [R].

The strategy is based on the algorithm originally described by Kapovich-Leeb-Porti in [KLP14], see also [KLP23].

Contact Max Riestenberg at riestenberg@mis.mpg.de with any questions.


# Prerequisites

You need a working version of Python.

Before you can run the computation, you first need to generate the list of words of length 8 e.g. using KBMAG [HT23], save them as a file somewhere, and change a line in the code to access that file. The list of words is unfortunately too large to store here. You can also contact me for it.

Note that the presentation is given by <a,b,c,d|acDbACdB=1>.

The code takes my about 35 minutes to run on my tablet.

# References
[HT23] Derek Holt and The GAP Team, kbmag, Knuth-Bendix on Monoids and Automatic Groups, Version 1.5.11 (2023)
(Refereed GAP package), https://gap-packages.github.io/kbmag.

[KLP14] Michael Kapovich, Bernhard Leeb, and Joan Porti, Morse actions of discrete groups on symmetric space, (2014).

[KLP23] Michael Kapovich, Bernhard Leeb, and Joan Porti, Morse actions of discrete groups on symmetric spaces: Local-to-global principle, (2023).

[R] J. Maxwell Riestenberg, Certifying Anosov representations, (2024).


