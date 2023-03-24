[![Build documentation](https://github.com/rouskinlab/dreem/actions/workflows/documentation.yaml/badge.svg?branch=main)](https://github.com/rouskinlab/dreem/actions/workflows/documentation.yaml)
[![Run tests](https://github.com/rouskinlab/dreem/actions/workflows/tests.yaml/badge.svg?branch=main)](https://github.com/rouskinlab/dreem/actions/workflows/tests.yaml)

# DREEM

Prof. Silvi Rouskin's [DREEM algorithm](https://www.nature.com/articles/s41586-020-2253-5).

## Documentation

The documentation is available on [Github Pages](https://rouskinlab.github.io/dreem).

## Contributors

Yves Martin, Scott Grote, Matthew Allan, Albéric de Lajarte.

## For developers

Before pushing to main, ensure that the code passes the tests:

```
pytest test/test_pipeline.py
```

and that the docs compile:

```
cd docs
make html
```
