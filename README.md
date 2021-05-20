# The Cash Paradox (Replication of Jiang & Shao in Julia)

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pedrovergaramerino.github.io/CashParadox.jl/stable) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pedrovergaramerino.github.io/CashParadox.jl/dev)

GitHub Actions : [![Build Status](https://github.com/pedrovergaramerino/CashParadox.jl/workflows/CI/badge.svg)](https://github.com/pedrovergaramerino/CashParadox.jl/actions?query=workflow%3ACI+branch%3Amaster)


[![Coverage Status](https://coveralls.io/repos/pedrovergaramerino/CashParadox.jl/badge.svg?branch=master)](https://coveralls.io/r/pedrovergaramerino/CashParadox.jl?branch=master)
[![codecov.io](http://codecov.io/github/pedrovergaramerino/CashParadox.jl/coverage.svg?branch=master)](http://codecov.io/github/pedrovergaramerino/CashParadox.jl?branch=master)

This package provides a `Julia` infrastructure to replicate the paper *[The Cash Paradox](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3367850)* by Jiang & Shao (2019).

## Installation

In your julia REPL type

```julia
] add https://github.com/pedrovergaramerino/CashParadox.jl.git
```

## Showcase

This package provides several functions to replicate main figures and tables displayed in the paper and its *[Online Appendix](https://ideas.repec.org/p/red/append/18-268.html)*

For instance, in order to get the results from table 1 it suffices to type

```julia
using CashParadox

Table1()

```

The documentation shows all replicated figures and tables, as the functions used to create them.
